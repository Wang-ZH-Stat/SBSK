
library(MASS)
library(GeneralizedHyperbolic)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# RKHS kernel
k <- function(x, y, sigma, kernel = "Gaussian"){
  if(kernel == "Gaussian"){
    return(exp(-sum((x-y)^2)/(2*sigma^2)))
  }
  if(kernel == "Laplacian"){
    return(exp(-sum(abs(x-y))/sigma))
  }
}

# Kernel for U-statistic
h <- function(x1, x2, y1, y2, sigma, kernel = "Gaussian"){
  return(k(x1,x2,sigma,kernel=kernel) + k(y1,y2,sigma,kernel=kernel) - 
           k(x1,y2,sigma,kernel=kernel) - k(x2,y1,sigma,kernel=kernel))
}

# An unbiased estimator of MMD^2
MMD2 <- function(x, y, sigam, kernel = "Gaussian"){
  res <- 0
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  for(i in c(1:n)){
    for(j in c(1:n)){
      if(i != j){
        res <- res + h(x[i,],x[j,],y[i,],y[j,],sigam,kernel=kernel)
      }
      else{
        res <- res
      }
    }
  }
  return(res/(n*(n-1)))
}

# The special function
v <- function(mu){
  num <- (2/mu)*(pnorm(mu/2)-0.5)
  den <- (mu/2)*pnorm(mu/2)+dnorm(mu/2)
  return(num/den)
}

# Average run length
get_ARL <- function(B0, b){
  mu <- b*sqrt(2*(2*B0-1)/(B0*(B0-1)))
  ET <- exp(b^2/2)/b*((2*B0-1)/(B0*(B0-1)*sqrt(2*pi))*v(mu))^-1
  return(ET)
}

# Threshold
get_b_scan <- function(ARL, B0){
  f <- function(b, ARL, B0) {return(ARL - get_ARL(B0, b))}
  root <- uniroot(f, c(2, 5), ARL = ARL, B = B0, tol = 1e-4)
  return(root$root)
}

# Variance of ZB
get_VarZB <- function(X, N, B0, sigma, iter = 10000, kernel = "Gaussian", improve = FALSE){
  X <- as.matrix(X)
  n <- nrow(X)
  res <- 0
  prob_list <- rep(1, n)
  for(i in c(1:iter)){
    if(improve == FALSE){
      id <- sample(n, 8)
    }
    if(improve == TRUE){
      id <- sample(n, 8, prob=prob_list)
      prob_list[id] <- prob_list[id] - 1/10000
      prob_list <- prob_list / sum(prob_list)
    }
    sam_X <- as.matrix(X[id,])
    h1 <- h(sam_X[1,],sam_X[2,],sam_X[3,],sam_X[4,], sigma, kernel=kernel)
    h2 <- h(sam_X[5,],sam_X[6,],sam_X[7,],sam_X[8,], sigma, kernel=kernel)
    res <- res + (h1^2+h2^2)/(2*N) + h1*h2*(N-1)/N
  }
  res <- res/iter * 2/(B0*(B0-1))
  return(res)
}

# Detect the change point
dete_scan <- function(X, N, B0, b, sigma, kernel = "Gaussian", improve = FALSE){
  X <- as.matrix(X)
  n <- nrow(X)
  VarZB <- get_VarZB(X[c(1:(0.5*n)), ], N, B0, sigma, kernel=kernel, improve=improve)
  for(t in c((0.5*n):n)){
    ZB <- 0
    Y <- as.matrix(X[c((t-B0+1):t), ])
    for(l in c(1:N)){
      idl <- sample(t-B0, B0)
      Xl <- as.matrix(X[idl, ])
      ZB <- ZB + MMD2(Xl, Y, sigma, kernel=kernel)/N
    }
    if(ZB/sqrt(VarZB) > b){
      return(c(t,1))
    }
    if(t == n){
      # fail to detect
      return(c(n,0))
    }
  }
}

# Hotelling T2 statistic
dete_T2 <- function(X, B0, b){
  X <- as.matrix(X)
  n <- nrow(X)
  for(t in c((0.5*n):n)){
    Y <- as.matrix(X[c((t-B0+1):t), ])
    Xr <- as.matrix(X[c(1:(t-B0)), ])
    Sigma0 <- cov(Xr)
    T2 <- as.numeric(B0*t((apply(Y,2,mean)-apply(Xr,2,mean)))%*%
                       solve(Sigma0)%*%(apply(Y,2,mean)-apply(Xr,2,mean)))
    if(T2 > b){
      return(c(t,1))
    }
    if(t == n){
      # fail to detect
      return(c(n,0))
    }
  }
}

get_b_T2 <- function(ARL, B0, range){
  for(b in range){
    num <- length(range)
    print(paste0("b=",b))
    sum_t <- 0
    sum_fal <- 0
    iter <- 500
    set.seed(999)
    for(it in c(1:iter)){
      X <- mvrnorm(200, mu = rep(0,10), Sigma = diag(10))
      t <- dete_T2(X, B0, b)
      sum_t <- sum_t + (t[1]-100)
      sum_fal <- sum_fal + t[2]
    }
    print(paste0("sum_t=",sum_t))
    print(paste0("sum_fal=",sum_fal))
    if(sum_fal != 0){
      avg_t <- sum_t/sum_fal
      print(paste0("avg_t=",avg_t))
      if(num == 1){
        return(c(b, avg_t))
      }
      else{
        if(avg_t > ARL){
          return(c(b, avg_t))
        }
        if(b == range[num]){
          return(c(b_max, avg_t))
        }
      }
    }
  }
}


# Generalized likelihood ratio
dete_GLR <- function(X, B0, b){
  X <- as.matrix(X)
  n <- nrow(X)
  for(t in c((0.5*n):n)){
    Y <- as.matrix(X[c((t-B0+1):t), ])
    Xr <- as.matrix(X[c(1:(t-B0)), ])
    Xt <- as.matrix(X[c(1:t), ])
    T2 <- t*log(det(cov(Xt))) - B0*log(det(cov(Y))) - (t-B0)*log(det(cov(Xr)))
    if(T2 > b){
      return(c(t,1))
    }
    if(t == n){
      # fail to detect
      return(c(n,0))
    }
  }
}

get_b_GLR <- function(ARL, B0, range){
  for(b in range){
    num <- length(range)
    print(paste0("b=",b))
    sum_t <- 0
    sum_fal <- 0
    iter <- 500
    set.seed(999)
    for(it in c(1:iter)){
      X <- mvrnorm(200, mu = rep(0,10), Sigma = diag(10))
      t <- dete_GLR(X, B0, b)
      sum_t <- sum_t + (t[1]-100)
      sum_fal <- sum_fal + t[2]
    }
    print(paste0("sum_t=",sum_t))
    print(paste0("sum_fal=",sum_fal))
    if(sum_fal != 0){
      avg_t <- sum_t/sum_fal
      print(paste0("avg_t=",avg_t))
      if(num == 1){
        return(c(b, avg_t))
      }
      else{
        if(avg_t > ARL){
          return(c(b, avg_t))
        }
        if(b == range[num]){
          return(c(b_max, avg_t))
        }
      }
    }
  }
}


