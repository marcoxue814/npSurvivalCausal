STE.ncb.core <- function(ind, K, lam1s, lam2s=1e4*lam1s, lower=1, upper=Inf,
                         thresh.ratio=1e-8, traceit=TRUE, w0=NULL, maxit=100,
                         maxit2=25, xtol_rel=1e-4, xtol_rel2=1e-2,
                         check=FALSE, full=FALSE){ 
  
  if (check){
    if (!requireNamespace("Rglpk", quietly = TRUE)){
      stop("Package \"Rglpk\" needed for SLP algorithm (check=TRUE) to work. Please install it.", call. = FALSE)
    }
  }
  
  N <- length(ind)
  
  # construct a vector
  tind <- as.logical(ind)
  n1 <- sum(tind)
  n2 <- sum(!tind)
  
  # lower and upper bound
  if (length(lower)==1) lower <- rep(lower, n1)
  if (length(upper)==1) upper <- rep(upper, n1)
  
  # construct P1 and q1
  e <- eigen(K)
  thresh <- e$values[1]*thresh.ratio
  eind <- (e$values>=thresh)
  r <- sum(eind)
  
  if (traceit)
    cat("Number of components used: r = ", r, "with threshod.ratio =", thresh.ratio, "\n")
  
  if (is.null(w0)){
    w0 <- lower+1
  }
  w1 = w0
  
  nlams <- length(lam1s)
  lam1s <- sort(lam1s) # lams is sorted
  outlist <- list()
  ws <- matrix(0, nr=nlams, nc=N)
  SNs <- array(dim=nlams)
  unorms <- array(dim=nlams)
  fittedus <- array(dim=c(nlams, nc=N))
  rmseTo1s <- array(dim=nlams)
  outlist <- NULL
  outlist2 <- NULL
  alg <- rep(1,nlams)
  obj1 <- array(dim=nlams)
  obj2 <- array(dim=nlams)
  fittedus_trim = array(dim=c(nlams, nc=N))
  ws_trim = matrix(0, nr=nlams, nc=N)
  obj_trim = array(dim=nlams)
  SNs_trim = array(dim=nlams)
  fittedfs_trim = array(dim=c(nlams, nc=N))
  fnorms = array(dim=nlams)
  
  if (length(lam2s)==1){
    lam2s <- rep(lam2s, nlams)
  }
  
  for (i in (1:nlams)){
    nP1 <- e$vectors[,eind]/sqrt(N) # P1/sqrt(N)
    nq1 <- -lam1s[i]*N/e$values[eind] # -lam1s[i]*N1*q1
    # reorder
    oo <- order(nq1)
    nq1 <- nq1[oo]
    nP1 <- nP1[,oo]
    
    if (traceit) cat("####", i, ":", "\tlam1 =", lam1s[i], "\tlam2 =", lam2s[i], "\n")
    
    res <- nloptr::nloptr(x0=w0, eval_f=eval_obj_grad, lb=lower, ub=upper,
                          N=N, r=r, tind=tind, nP1=nP1, nq1=nq1, lam2=lam2s[i],
                          opts=list(algorithm="NLOPT_LD_LBFGS", xtol_rel=xtol_rel,
                                    print_level=0, maxeval=maxit, check_derivatives=F))
    
    if (full)  outlist[[i]] <- res
    obj1[i] <- res$objective
    ws[i,tind] <- res$sol
    temp <- t(nP1) %*% (ws[i,]-1)
    ee <- eigen(temp %*% t(temp) + diag(nq1))
    
    ## checking if the largest eigenvalue has multiplicity 1
    if (check){
      if (abs(ee$values[1]-ee$values[2])/abs(ee$values[1]) < 1e-6){
        cat("multiplicity issue of largest eigenvalues\n")
        res2 <- mineig(w0=res$sol, N=N, r=r, tind=tind, nP1=nP1, nq1=nq1, lam2=lam2s[i],
                       lower=lower, upper=upper, rho0=1e-3, tau0=1e-6,
                       tol=xtol_rel2, max.iter=maxit2, verbose=T)
        if (full) outlist2[[i]] <- res2
        obj2[i] <- eval_obj_grad(res2$w, N=N, r=r, tind=tind, nP1=nP1, nq1=nq1, lam2=lam2s[i])$objective
        
        # udpate w
        ws[i,tind] <- res2$w
        temp <- t(nP1) %*% (ws[i,]-1)
        ee <- eigen(temp %*% t(temp) + diag(nq1))
        alg[i] <- 2
      }
    }
    
    # compute SN
    SNs[i] <- (sum(temp*ee$vectors[,1]))^2
    unorms[i] <- sqrt(sum(nq1/(-lam1s[i]) * ee$vectors[,1]^2))
    fittedus[i,] <- as.vector(nP1 %*% ee$vectors[,1] * (N))
    rmseTo1s[i] <- sqrt(mean((ws[i,]-1)^2))
    
    # trim fitted u
    fittedus_trim[i,] = fittedus[i,]
    ulind = as.logical(fittedus[i,]<0)
    fittedus_trim[i, ulind] = 1e-6
    uuind = as.logical(fittedus[i,]>1)
    fittedus_trim[i, uuind] = 1 - 1e-6
    
    # calculate ws based on truncated u
    K_inv = eigen(K)$vectors %*% diag( 1 / eigen(K)$values) %*% solve(eigen(K)$vectors)
    alpha_trim = as.vector(K_inv %*% log(fittedus_trim[i,] / (1-fittedus_trim[i,])))
    
    res2 = nloptr::nloptr(x0 = res$sol, eval_f = outer_obj_grad, lb=lower, ub=upper,
                          N=N, tind=tind, lam1s = lam1s[i], lam2s=lam2s[i], alpha = alpha_trim, K = K,
                          opts=list(algorithm="NLOPT_LD_LBFGS", xtol_rel = xtol_rel,
                                    print_level=0, maxeval=maxit, check_derivatives=F)
    
    ws_trim[i,tind] = res2$solution
    obj_trim[i] = res2$objective
    
    res3 = nloptr::nloptr(x0 = alpha_trim, eval_f = inner_obj_grad, w = ws_trim[i, tind], tind = tind, K=K, N=N, lam1s=lam1s[i],
                          opts=list(algorithm="NLOPT_LD_LBFGS", xtol_rel = 1e-4, maxeval = 200))
    
    alpha_trim_h = res3$solution
    fittedfs_trim[i,] = c(logitinv(K %*% alpha_trim_h))
    fnorms[i] = sum(fittedfs_trim[i,]^2)/N
    
    z <- rep(-1, N)
    z[tind] <- ws_trim[i, tind]-1
    SNs_trim[i] = (z %*% logitinv(K %*% alpha_trim_h)/N)^2
  }
  
  return(list(outlist=outlist, outlist2=outlist2, ws=ws, SNs=SNs, unorms=unorms,fnorms = fnorms, ws_trim = ws_trim, obj_trim = obj_trim, SNs_trim = SNs_trim,
              fittedus=fittedus, fittedfs_trim = fittedfs_trim, rmseTo1s=rmseTo1s, lam1s=lam1s, lam2s=lam2s, alg=alg, obj1=obj1, obj2=obj2,
              check=check))
}


inner_obj_grad = function(alpha, w, K, N, tind, lam1s) {
  z <- rep(-1, N)
  z[tind] <- w-1
  # S1 = S2 = 0
  grad = NULL
  for(j in 1:N) {
    S1 = sum(z * logitinv(K %*% alpha))
    S2 = sum(z * logitinv(K%*%alpha) * (1-logitinv(K%*%alpha)) * K[, j])
    grad = c(grad, -2/N^2*S1*S2)
  }
  grad = grad + lam1s * 2 * K %*% alpha
  obj = -c(1/N^2*t(logitinv(K %*% alpha))%*%(z %*% t(z)) %*% logitinv(K %*% alpha) - lam1s*t(alpha) %*% K %*% alpha)
  
  return(list("objective" = obj,
              "gradient" = as.vector(grad)))
}

outer_obj_grad = function(w, N, K, alpha, tind, lam1s, lam2s) {
  inner_max = nloptr(x0 = alpha, eval_f = inner_obj_grad, w = w, tind = tind, K=K, N=N, lam1s=lam1s,
                     opts=list(algorithm="NLOPT_LD_LBFGS", xtol_rel = 1e-2, maxeval = 25))
  
  alpha_h = inner_max$solution
  z <- rep(-1, N)
  z[tind] <- w-1
  
  S = sum(z * logitinv(K %*% alpha_h))
  
  grad = 2/N^2 * S * logitinv(K %*% alpha_h)[tind] + lam2s * 2/N * w
  return(list("objective" = -inner_max$objective + lam2s*sum(w^2)/N,
              "gradient" = as.vector(grad)))
}

STE.ncb.SN <- function(ind, K, lam1s, lam2s=1e4*lam1s, lower=1, upper=Inf,
                       thresh.ratio=1e-8, traceit=TRUE, w0=NULL, maxit=100,
                       maxit2=25, xtol_rel=1e-4, xtol_rel2=1e-2, method=1,
                       check=FALSE, full=FALSE){
  ores <- STE.ncb.core(ind=ind, K=K, lam1s, lam2s=lam2s, lower=lower, upper=upper,
                       thresh.ratio=thresh.ratio, traceit=traceit, w0=w0,
                       maxit=maxit, maxit2=maxit2, xtol_rel=xtol_rel, xtol_rel2=xtol_rel2,
                       check=check, full=full)
  
  # find which SN is the smallest
  if (method==1){
    ind <- which.min(ores$SNs_trim)
  } else if (method==2){
    ind <- which(diff(ores$SNs_trim)/diff(ores$lam1s) > (-1e-6))[1]
  }
  return(list(w=ores$ws_trim[ind,], ind=ind, warns=c(ind==1, ind==length(lam1s)), ores=ores))
}
