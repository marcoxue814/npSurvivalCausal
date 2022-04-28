Q = 200 #simulation runs
simu_data = replicate(Q, list())
for(q in 1:Q){
  n = 200 #sample size
  set.seed(12345 + 2 * q)
  
  Z = matrix(rnorm(4*n),ncol=4,nrow=n)
  
  # overlap parameter
  alpha = 2 # Corresponding overlap is weak
  
  prop = 1 / (1 + exp(alpha * (1.2 * Z[,1] - 1.5 * Z[,2] - Z[,4])))
  treat = rbinom(n, 1, prop)
  
  # event time
  lambda = 0.5 + 0.1 * treat + Z[,1]^2 + 0.1 * exp(Z[,1] * Z[,2]) + 0.2 * (Z[, 1] + Z[, 3] + 6)^2
  tt = rexp(n, rate = lambda)
  
  # covariates
  X = cbind(exp(Z[,1]/2),  Z[,2]/(1+exp(Z[,1]))+10,
            (Z[,1]*Z[,3]/25+0.6)^3, (Z[,2]+Z[,4]+20)^2)
  
  # censoring time
  C = rexp(n, rate = 0.5) # Corresponding censoring rate is 10%
  
  # observed time
  Y = pmin(tt, C)
  Delta = (tt <= C) * 1
  
  
  simu_data[[q]]$Y = Y
  simu_data[[q]]$X = X
  simu_data[[q]]$treat = treat
  simu_data[[q]]$prop = prop
  simu_data[[q]]$censor = C
  simu_data[[q]]$event = tt
  simu_data[[q]]$hazard = lambda
  simu_data[[q]]$failures = Delta
  
}

save(simu_data, file = "ste_simu_data_weak.rda")

Q = 200
load(file = "ste_simu_data_weak.rda")
results = replicate(Q, list())
for(q in 1:Q){
  print(paste0("Working on iter=",q))
  
  X = simu_data[[q]]$X
  treat = simu_data[[q]]$treat

  # Sobolev kernel
  Xstd = transform.sob(X)$Xstd # standardize X to [0,1]^p
  K = getGram(Xstd) # get Gram matrix using Sobolev kernel

  # design a grid for the tuning parameter
  nlam1 <- 25
  # lambda for DTE
  lams1 <- exp(seq(log(1e-8), log(1e-3), len=nlam1))

  # # kernel based methods for T = 1
  fit1_prop = DTE.ncb.SN(treat, K, lam1s=lams1)

  # # kernel based methods for T = 0
  fit0_prop = DTE.ncb.SN(1-treat, K, lam1s=lams1)


  # estimation of propensity using logistic regression
  fit_log = glm(treat ~ X, family = 'binomial')
  
  nlam2 <- 50
  # lambda for DTE
  lams2 <- exp(seq(log(1e-8), log(1), len=nlam2))
  
  # # kernel based methods for T = 1
  fit1_init = ATE.ncb.SN(treat, K, lam1s=lams2)
  
  # # kernel based methods for T = 0
  fit0_init = ATE.ncb.SN(1-treat, K, lam1s=lams2)


  results[[q]]$fit1_prop = fit1_prop
  results[[q]]$fit0_prop = fit0_prop
  results[[q]]$fit1_init = fit1_init
  results[[q]]$fit0_init = fit0_init
  results[[q]]$fit_log = fit_log
}

save(results, file = "results_weak_weights.rda")

# get estimation of the true cdf of Y1 and Y0
n = 10000
surv_Y1 = matrix(0, nrow = n, ncol = 100)
surv_Y0 = matrix(0, nrow = n, ncol = 100)
for(i in 1:n){
  set.seed(12345 + 2 * i)
  print(paste0("Working on iter=",i))
  N = 1e5
  Z <- matrix(rnorm(4*N),ncol=4,nrow=N)
  alpha = 2
  prop <- 1 / (1 + exp(alpha * (1.2 * Z[,1] - 1.5 * Z[,2] - Z[,4])))
  treat <- rbinom(N, 1, prop)
  
  # event time
  lambda = 0.5 + 0.1 * treat + Z[,1]^2 + 0.1 * exp(Z[,1] * Z[,2] ) + 0.2 * (Z[, 1] + Z[, 3] + 6)^2
  tt = rexp(N, rate = lambda)
  
  Y1 = tt[as.logical(treat)]  
  Y0 = tt[as.logical(1-treat)]
  
  cdf_Y1 = ecdf(Y1)
  cdf_Y0 = ecdf(Y0)
  
  surv_Y1[i,] = 1 - cdf_Y1(seq(0, 1, length.out = 100))
  surv_Y0[i,] = 1 - cdf_Y0(seq(0, 1, length.out = 100)) 
}
# Get Monte Carlo simulation estimation of the true functions F1 and F0
surv_Y1_true = apply(surv_Y1, 2, mean)
surv_Y0_true = apply(surv_Y0, 2, mean)
save(surv_Y1_true, file = "true_Y1_times_weak.rda")
save(surv_Y0_true, file = "true_Y0_times_weak.rda")

load("~/true_Y1_times_weak.rda")
load("~/true_Y0_times_weak.rda")


Q = 200
sup_prop0 = rep(0, Q)
l2_prop0 = rep(0, Q)
sup_prop1 = rep(0, Q)
l2_prop1 = rep(0, Q)

sup_AKME0 = rep(0, Q)
l2_AKME0 = rep(0, Q)
sup_AKME1 = rep(0, Q)
l2_AKME1 = rep(0, Q)

sup_init0 = rep(0, Q)
l2_init0 = rep(0, Q)
sup_init1 = rep(0, Q)
l2_init1 = rep(0, Q)

prop_fn1_value = matrix(0, nrow = Q, ncol = 100)
prop_fn0_value = matrix(0, nrow = Q, ncol = 100)
AKME_fn1_value = matrix(0, nrow = Q, ncol = 100)
AKME_fn0_value = matrix(0, nrow = Q, ncol = 100)
init_fn1_value = matrix(0, nrow = Q, ncol = 100)
init_fn0_value = matrix(0, nrow = Q, ncol = 100)


for (q in 1:Q){
  print(paste0("Working on iter=",q))
  
  X = simu_data[[q]]$X
  Y = simu_data[[q]]$Y
  treat = simu_data[[q]]$treat
  prop = simu_data[[q]]$prop
  Delta = simu_data[[q]]$failures
  
  fit1_prop = results[[q]]$fit1_prop
  fit0_prop = results[[q]]$fit0_prop
  fit_log = results[[q]]$fit_log
  
  # proposed weight
  w1_prop = fit1_prop$w
  w0_prop = fit0_prop$w
  w_prop = w1_prop + w0_prop
  
  # initial weight
  w1_init = fit1_init$w
  w0_init = fit0_init$w
  w_init = w1_init + w0_init
  
  # AKME weight
  phat = fit_log$fitted.values
  w_AKME = (treat == 1) * (1/phat) + (treat == 0) * (1) / (1 - phat)
  
  # Observed event time in treatment group
  observed_event_time1 = Y[Delta*treat == 1]
  n_observed_event_time1 = length(observed_event_time1)
  
  # Observed event time in control group
  observed_event_time0 = Y[Delta*(1-treat) == 1]
  n_observed_event_time0 = length(observed_event_time0)
  
  # AKME survival probabilities
  AKME_surv = ipw.survival(times = Y, failures = Delta, variable = treat, weights = w_AKME)
  
  #Proposed survival probabilities
  prop_surv = ipw.survival(times = Y, failures = Delta, variable = treat, weights = w_prop)
  
  # AKME S_Z(1)
  AKME_surv_prob1 = subset(as.data.frame(AKME_surv$table.surv), variable == 1)$survival
  AKME_surv_prob1 = AKME_surv_prob1[1:(n_observed_event_time1+1)]
  
  # Estimated survival function of S_Z(1) from AKME
  AKME_fn1 = stepfun(sort(observed_event_time1), 
                     AKME_surv_prob1, right = FALSE)
  
  # AKME S_Z(0)
  AKME_surv_prob0 = subset(as.data.frame(AKME_surv$table.surv), variable == 0)$survival
  AKME_surv_prob0 = AKME_surv_prob0[1:(n_observed_event_time0+1)]
  
  # Estimated survival function of S_Z(0) from AKME
  AKME_fn0 = stepfun(sort(observed_event_time0), 
                     AKME_surv_prob0, right = FALSE)
  
  # Proposed S_Z(1)
  prop_surv_prob1 = subset(as.data.frame(prop_surv$table.surv), variable == 1)$survival
  prop_surv_prob1 = prop_surv_prob1[1:(n_observed_event_time1+1)]
  
  # Estimated survival function of S_Z(1) from proposed method
  prop_fn1 = stepfun(sort(observed_event_time1), 
                     prop_surv_prob1, right = FALSE)
  
  # Proposed S_Z(0)
  prop_surv_prob0 = subset(as.data.frame(prop_surv$table.surv), variable == 0)$survival
  prop_surv_prob0 = prop_surv_prob0[1:(n_observed_event_time0+1)]
  
  # Estimated survival function of S_Z(0) from proposed method
  prop_fn0 = stepfun(sort(observed_event_time0), 
                     prop_surv_prob0, right = FALSE)
  
  #Proposed survival probabilities with initial weights
  init_surv = ipw.survival(times = Y, failures = Delta, variable = treat, weights = w_init)
  
  
  # init S_Z(1)
  init_surv_prob1 = subset(as.data.frame(init_surv$table.surv), variable == 1)$survival
  init_surv_prob1 = init_surv_prob1[1:(n_observed_event_time1+1)]
  
  # Estimated survival function of S_Z(1) from init
  init_fn1 = stepfun(sort(observed_event_time1),
                     init_surv_prob1, right = FALSE)
  
  # init S_Z(0)
  init_surv_prob0 = subset(as.data.frame(init_surv$table.surv), variable == 0)$survival
  init_surv_prob0 = init_surv_prob0[1:(n_observed_event_time0+1)]
  
  # Estimated survival function of S_Z(0) from init
  init_fn0 = stepfun(sort(observed_event_time0),
                     init_surv_prob0, right = FALSE)
  
  
  prop_fn1_value[q, ] = prop_fn1(seq(0, 1, length.out = 100))
  prop_fn0_value[q, ] = prop_fn0(seq(0, 1, length.out = 100))
  
  AKME_fn1_value[q, ] = AKME_fn1(seq(0, 1, length.out = 100))
  AKME_fn0_value[q, ] = AKME_fn0(seq(0, 1, length.out = 100))
  
  init_fn1_value[q, ] = init_fn1(seq(0, 1, length.out = 100))
  init_fn0_value[q, ] = init_fn0(seq(0, 1, length.out = 100))
  
  delta_Y1_true = rep(1/100, 100)
  delta_Y0_true = rep(1/100, 100)
  
  
  
  sup_prop1[q] = max(abs(prop_fn1(seq(0, 1, length.out = 100)) - surv_Y1_true))
  l2_prop1[q] = sqrt( sum((prop_fn1(seq(0, 1, length.out = 100)) - surv_Y1_true)^2 * delta_Y1_true ))
  # 
  sup_AKME1[q] = max(abs(AKME_fn1(seq(0, 1, length.out = 100)) - surv_Y1_true))
  l2_AKME1[q] = sqrt( sum((AKME_fn1(seq(0, 1, length.out = 100)) - surv_Y1_true)^2 * delta_Y1_true ))
  
  sup_init1[q] = max(abs(init_fn1(seq(0, 1, length.out = 100)) - surv_Y1_true))
  l2_init1[q] = sqrt( sum((init_fn1(seq(0, 1, length.out = 100)) - surv_Y1_true)^2 * delta_Y1_true ))
  
  sup_prop0[q] = max(abs(prop_fn0(seq(0, 1, length.out = 100)) - surv_Y0_true))
  l2_prop0[q] = sqrt( sum((prop_fn0(seq(0, 1, length.out = 100)) - surv_Y0_true)^2 * delta_Y0_true))
  # 
  sup_AKME0[q] = max(abs(AKME_fn0(seq(0, 1, length.out = 100)) - surv_Y0_true))
  l2_AKME0[q] = sqrt( sum((AKME_fn0(seq(0, 1, length.out = 100)) - surv_Y0_true)^2 * delta_Y0_true))
  
  sup_init0[q] = max(abs(init_fn0(seq(0, 1, length.out = 100)) - surv_Y0_true))
  l2_init0[q] = sqrt( sum((init_fn0(seq(0, 1, length.out = 100)) - surv_Y0_true)^2 * delta_Y0_true))
}


# Boxplots of supremum distance 
df1 = data.frame(SupNorm = c(sup_prop1, sup_AKME1, sup_init1), 
                 Estimator =  c(rep("Proposed", 200), 
                                rep("AKME", 200),
                                rep("Initial", 200)))
df1$Estimator = factor(df1$Estimator,     # Reorder factor levels
                         c("Initial", "AKME", "Proposed"))
p1 = ggplot(df1, aes(x=Estimator, y=SupNorm)) + 
  geom_boxplot() +
  labs(x = "", y = "") + 
  theme(axis.text=element_text(size=20),
         axis.title=element_text(size=20))
p1 

df2 = data.frame(SupNorm = c(sup_prop0, sup_AKME0, sup_init1), 
                 Estimator =  c(rep("Proposed", 200), 
                                rep("AKME", 200),
                                rep("Initial", 200)))
df2$Estimator = factor(df2$Estimator,     # Reorder factor levels
                       c("Initial", "AKME", "Proposed"))
p2 = ggplot(df2, aes(x=Estimator, y=SupNorm)) + 
  geom_boxplot() +
  labs(x = "", y = "") + 
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20))
p2 




delta_Y1_true = rep(1/100, 100)

delta_Y0_true = rep(1/100, 100)

# Calculate ise

ise_prop1 = apply(prop_fn1_value, 1, function(x) sum(delta_Y1_true * (x - surv_Y1_true)^2))

ise_prop0 = apply(prop_fn0_value, 1, function(x) sum(delta_Y0_true * (x - surv_Y0_true)^2))

ise_akme1 = apply(AKME_fn1_value, 1, function(x) sum(delta_Y1_true * (x - surv_Y1_true)^2))

ise_akme0 = apply(AKME_fn0_value, 1, function(x) sum(delta_Y0_true * (x - surv_Y0_true)^2))

ise_init1 = apply(init_fn1_value, 1, function(x) sum(delta_Y1_true * (x - surv_Y1_true)^2))

ise_init0 = apply(init_fn0_value, 1, function(x) sum(delta_Y0_true * (x - surv_Y0_true)^2))

# Calculate bias
bias_prop1 = apply(prop_fn1_value, 2, mean) - surv_Y1_true

bias_prop0 = apply(prop_fn0_value, 2, mean) - surv_Y0_true

bias_AKME1 = apply(AKME_fn1_value, 2, mean) - surv_Y1_true

bias_AKME0 = apply(AKME_fn0_value, 2, mean) - surv_Y0_true

bias_init1 = apply(init_fn1_value, 2, mean) - surv_Y1_true

bias_init0 = apply(init_fn0_value, 2, mean) - surv_Y0_true


# Calculate ASAE
asae_prop1 = max(abs(bias_prop1))

asae_prop0 = max(abs(bias_prop0))

asae_akme1 = max(abs(bias_AKME1))

asae_akme0 = max(abs(bias_AKME0))

asae_init1 = max(abs(bias_init1))

asae_init0 = max(abs(bias_init0))


#Calculate ISB
isb_prop1 = sum(delta_Y1_true * bias_prop1^2)

isb_prop0 = sum(delta_Y0_true * bias_prop0^2)

isb_akme1 = sum(delta_Y1_true * bias_AKME1^2)

isb_akme0 = sum(delta_Y0_true * bias_AKME0^2)

isb_init1 = sum(delta_Y1_true * bias_init1^2)

isb_init0 = sum(delta_Y0_true * bias_init0^2)


#Calculate MISE
aise_prop1 = mean(ise_prop1)

aise_prop0 = mean(ise_prop0)

aise_akme1 = mean(ise_akme1)

aise_akme0 = mean(ise_akme0)

aise_init1 = mean(ise_init1)

aise_init0 = mean(ise_init0)



# Calculate SD of ISE
sd_ise_prop1 = sd(ise_prop1)

sd_ise_prop0 = sd(ise_prop0)

sd_ise_akme1 = sd(ise_akme1)

sd_ise_akme0 = sd(ise_akme0)

sd_ise_init1 = sd(ise_init1)

sd_ise_init0 = sd(ise_init0)


# Boxplot of ISE

df5 = data.frame(ISE = c(ise_prop1, ise_akme1, ise_init1), 
                 Estimator =  c(rep("Proposed", 200), 
                                rep("AKME", 200),
                                rep("Initial", 200)))
df5$Estimator = factor(df5$Estimator,     # Reorder factor levels
                       c("Initial", "AKME", "Proposed"))
p5 = ggplot(df5, aes(x=Estimator, y=ISE)) + 
  geom_boxplot() +
  labs(x = "", y = "") + 
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20))
p5


df6 = data.frame(ISE = c(ise_prop0, ise_akme0, ise_init0), 
                 Estimator =  c(rep("Proposed", 200), 
                                rep("AKME", 200),
                                rep("Initial", 200)))
df6$Estimator = factor(df6$Estimator,     # Reorder factor levels
                       c("Initial", "AKME", "Proposed"))
p6 = ggplot(df6, aes(x=Estimator, y=ISE)) + 
  geom_boxplot() +
  labs(x = "", y = "") + 
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20))
p6 



table1 = cbind(rbind(mean(sup_prop1) * 100, mean(sup_AKME1) * 100, mean(sup_init1) * 100),
               rbind(sd(sup_prop1) * 100, sd(sup_AKME1) * 100, sd(sup_init1) * 100),
               rbind(asae_prop1 * 100, asae_akme1 * 100, asae_init1 * 100),
               rbind(aise_prop1 * 1e4, aise_akme1 * 1e4, aise_init1 * 1e4),
               rbind(sd_ise_prop1 * 1e4, sd_ise_akme1 * 1e4, sd_ise_init1 * 1e4),
               rbind(isb_prop1 * 1e4, isb_akme1 * 1e4, isb_init1 * 1e4),
               rbind(mean(sup_prop0) * 100, mean(sup_AKME0) * 100, mean(sup_init0) * 100),
               rbind(sd(sup_prop0) * 100, sd(sup_AKME0) * 100, sd(sup_init0) * 100),
               rbind(asae_prop0 * 100, asae_akme0 * 100, asae_init0 * 100),
               rbind(aise_prop0 * 1e4, aise_akme0 * 1e4, aise_init0 * 1e4),
               rbind(sd_ise_prop0 * 1e4, sd_ise_akme0 * 1e4, sd_ise_init0 * 1e4),
               rbind(isb_prop0 * 1e4, isb_akme0 * 1e4, isb_init0 * 1e4)
              
)
rownames(table1) = c("Proposed", "AKME", "Initial")
colnames(table1) = c("sup_mean", "sup_sd", "asae", 
                     "ise_mean", "ise_sd", "isb",
                     "sup_mean", "sup_sd", "asae", 
                     "ise_mean", "ise_sd", "isb")
library(xtable)
xtable(table1, digits = 1)






