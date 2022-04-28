#################
#### kernels ####
#################

#########################################
#########################################

K.gauss <- function(x,y,sigma=1){
  return(exp(-sum((x-y)^2) / sigma^2))
}

#########################################
#########################################

# standardize L2-norm=1
K.gauss2 <- function(x,y,sigma=1){
  return(dnorm( sqrt(sum((x-y)^2)), mean=0, sd=sigma))
}

#########################################
#########################################

K.laplace <- function(x,y,alpha){
  return(exp(-alpha*sum(abs(x-y))))
}


################################
##### Sobolev RK function ######
################################
# not optimized

#########################################
#########################################

k1 <- function(t){
  return(t-.5)
}

#########################################
#########################################

k2 <- function(t){
  return( (k1(t)^2-1/12)/2 )
}

#########################################
#########################################

k4 <- function(t){
  return( (k1(t)^4-k1(t)^2/2+7/240)/24 )
}


#########################################
#########################################

K.sob <- function(s,t){
  ans <- 1 + k1(s)*k1(t) + k2(s)*k2(t) - k4(abs(s-t))
  return(ans)
}

#########################################
#########################################

K.sob.prod <- function(x, y){
  p <- length(x)
  out <- 1
  for (j in (1:p)){
    out <- out*K.sob(x[j], y[j])
  }
  return(out)
}

#########################################
#########################################

standard <- function(x){
  return((x-mean(x))/sd(x))
}

#########################################
#########################################

transform.sob <- function(X){
  Xlim <- apply(X, 2, range)
  Xstd <- matrix(nr=nrow(X), nc=ncol(X))
  for (i in (1:ncol(X))){
    Xstd[,i] <- (X[,i]-Xlim[1,i])/diff(Xlim[,i])
  }
  return(list(Xstd=Xstd, Xlim=Xlim))
}


##################################################
#### general function to produce Gram matrix #####
##################################################

getGram <- function(X, standardize=T, ker="sob", sigma=1, alpha=1){
  # note that sobolev kernel is not compatible with standardize=T, and it requires X to be
  # transformed to [0,1]^p
  if (ker=="sob") standardize <- F
  if (standardize) X <- apply(X, 2, standard)
  N <- nrow(X)
  K <- matrix(nr=N, nc=N)
  for (i in (1:N)){
    for (j in (1:N)){
      if (ker=="gauss"){
        K[i,j] <- K.gauss(X[i,], X[j,], sigma=sigma)
      } else if (ker=="gauss2"){
        K[i,j] <- K.gauss2(X[i,], X[j,], sigma=sigma)
      } else if (ker=="laplace"){
        K[i,j] <- K.laplace(X[i,], X[j,], alpha=alpha)
      } else if (ker=="sob"){
        K[i,j] <- K.sob.prod(X[i,], X[j,])
      }
    }
  }
  return(K)
}

logitinv = function(x) {
  return(1 / (1 + exp(-x)))
}


dpr1_eigval <- function(i, d, z, need.order=T, vec=F){
  if (need.order){
    o <- order(d)
    d <- d[o]
    z <- z[o]
  }
  rho <- sum(z^2)
  z <- z/sqrt(rho)
  res <- .Call("dpr1_eigval", as.integer(i), d, z, as.double(rho), PACKAGE="ATE.ncb")
  if (vec){
    vv <- 1/res[[1]]*z
    res[[4]] <- vv/sqrt(sum(vv^2))
  }
  return(res)
}


# compute the full eigendecomposition for matrix diag(d) + outer(z, z)
dpr1_eigen <- function(d, z, Kstart=1, Kstop=length(d), need.order=T, order.eigen=T){
  if (need.order){
    o <- order(d)
    d <- d[o]
    z <- z[o]
  }
  rho <- sum(z^2)
  z <- z/sqrt(rho)
  res <- .Call("dpr1_eigen", as.integer(Kstart), as.integer(Kstop), d, z, as.double(rho), PACKAGE="ATE.ncb")
  if (order.eigen){
    nn <- ncol(res[[1]])
    res[[1]] <- res[[1]][,nn:1]
    res[[2]] <- res[[2]][nn:1]
  }
  names(res) <- c("vectors", "values", "info")
  return(res)
}

#################
#### ATE.ncb ####
#################

#########################################
#########################################

eval_obj_grad <- function(w, N, r, tind, nP1, nq1, lam2){
  z <- rep(-1, N)
  z[tind] <- w-1
  V <- sum(w^2)/N
  z <- t(nP1) %*% z
  dpr1 <- dpr1_eigval(r, nq1, z, need.order=F, vec=T)
  v <- dpr1[[4]]
  #v <- eigen(z %*% t(z) + diag(nq1))$vectors[,1]
  return(list("objective"=dpr1[[2]] + lam2*V,
              "gradient"=as.vector(2 * as.vector(t(v) %*% z) * nP1[tind,] %*%
                                     v) + 2*lam2*(w)/N)) # modified
}
