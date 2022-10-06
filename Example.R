set.seed(123)

n = 200 #sample size
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


#------------------------------------
# RKHS-based covariate balancing
#------------------------------------
# Sobolev kernel
Xstd = transform.sob(X)$Xstd # standardize X to [0,1]^p
K = getGram(Xstd) # get Gram matrix using Sobolev kernel

# design a grid for the tuning parameter
nlam1 = 25
# lambda for DTE
lams1 = exp(seq(log(1e-8), log(1e-3), len=nlam1))

# kernel based methods for D = 1
fit1 = STE.ncb.SN(treat, K, lam1s=lams1)
# kernel based methods for D = 0
fit0 = STE.ncb.SN(1-treat, K, lam1s=lams1)

# RKHS-based covariate balancing weights
w = fit1$w + fit0$w

library(RISCA)
surv = ipw.survival(times = Y, failures = Delta, variable = treat, weights = w)