###--------------------------------------------------------------------------###
###-- Gaussian process vector autoregressions and macroeconomic uncertainty -###
###---------------- Hauzenberger, Huber, Marcellino, & Petz -----------------###
###--------------- Journal of Business & Economic Statistics-----------------###
###--------------------------------------------------------------------------###
###-------------------------- GP-VAR code files -----------------------------###
###------------------------- Eq-by-eq estimation ----------------------------###
###--------------------------------------------------------------------------###
rm(list = ls())

w.dir <- ""
est.dir   <-  paste0(w.dir, "eqbyeq/"); dir.create(est.dir, showWarnings = FALSE)
girf.dir  <-  paste0(w.dir, "girfs/") ; dir.create(girf.dir, showWarnings = FALSE)
func.dir  <-  paste0(w.dir, "gpvar_funcs/")

library(MASS)
library(snow)
library(snowfall)
library(RcppArmadillo)
library(Rcpp)
library(coda)
library(stochvol)
library(Matrix)
library(mvtnorm)
library(truncnorm)
library(zoo)
library(readr)

###--------------------------------------------------------------------------###
###---------------------- Auxiliary functions -------------------------------###
###--------------------------------------------------------------------------###
source(paste0(func.dir, "gp_eqbyeq_mcmc.R"))  
sourceCpp(paste0(func.dir, "sqexp_kernel.cpp"))

###--------------------------------------------------------------------------###
###--------------------------------- Set-up ---------------------------------###
###--------------------------------------------------------------------------###
model.setup <-  list(
  data.set   = "BB",  # Basu & Bundick dataset
  stdz       = TRUE,  # Standardize data
  unc.ind    = "VXO", # Uncertainty indicator
  p          =  4  ,  # No. of lags
  sv         = "homo",  # Stochastic volatility
  hyper.grid =  list("h"= seq(0.1 , 10, length.out =  50) , "sig" = seq(0.5,5, length.out =  20)),
  h.sc       =  "med-heur",               # Median heuristic
  c          =  c("own" = 2, "other" = 2) # Shrinkage
)

mcmc.setup <-  list( # MCMC preliminaries
  nsave = 2500,
  nburn = 2500,
  nthin = 4
) 

list2env(model.setup, .GlobalEnv)
list2env(mcmc.setup, .GlobalEnv)

###--------------------------------------------------------------------------###
###---------------------- Basu & Bundick data set-up ------------------------###
###--------------------------------------------------------------------------###
basu_bundick_data_i <- read.csv(paste0(w.dir, "data/BB_realization.csv"))
matplot(basu_bundick_data_i, type = "l")
var.set <- colnames(basu_bundick_data_i)

basu_bundick_data_i[,1] <- 100*log(basu_bundick_data_i[,1])
basu_bundick_data_i[,2] <- 100*log(basu_bundick_data_i[,2])
basu_bundick_data_i[,3] <- 100*log(basu_bundick_data_i[,3])
basu_bundick_data_i[,4] <- 100*log(basu_bundick_data_i[,4])
basu_bundick_data_i[,5] <- 100*log(basu_bundick_data_i[,5])
basu_bundick_data_i[,6] <- 100*(basu_bundick_data_i[,6] - 1)
basu_bundick_data_i[,7] <- 100*log(basu_bundick_data_i[,7])
yraw <- basu_bundick_data_i
MM <- ncol(yraw)

if(stdz){
  sc.mean <- apply(yraw,2, mean)
  sc.sd <- apply(yraw,2, sd)
  yraw <- apply(yraw,2,function(x) (x-mean(x))/sd(x))
}else{
  sc.sd <- rep(1, MM)
  sc.mean <- apply(yraw,2, mean)
  yraw <- apply(yraw,2,function(x) x-mean(x))
}

Xraw <- mlag(yraw,p)
Y <- as.matrix(yraw[(p+1):nrow(yraw),])
X <- as.matrix(Xraw[(p+1):nrow(Xraw),])

N <- nrow(Y)
K <- ncol(X)

## OLS variance scalings
sigma_sq <- sigma_sq_org <- matrix(0,MM,1) #vector which stores the residual variance
rownames(sigma_sq) <- rownames(sigma_sq_org) <- colnames(Y)

for (i in 1:MM){
  Ylag_i <- embed(Y[,i], p+1)
  Y_i <- Ylag_i[,1]
  Ylag_i <- Ylag_i[,2:(p+1)]
  alpha_i <- solve(t(Ylag_i)%*%Ylag_i)%*%t(Ylag_i)%*%Y_i
  sigma_sq[i,1] <- sigma_sq_org[i,1] <- (1/(nrow(Y)-p))*t(Y_i-Ylag_i%*%alpha_i)%*%(Y_i-Ylag_i%*%alpha_i) 
}

###--------------------------------------------------------------------------###
###------------------------------- Directories ------------------------------###
###--------------------------------------------------------------------------###
dir <- paste0(w.dir, data.set, "-", unc.ind, "-p", p)
dir.create(dir, showWarnings = FALSE)
dir <- paste0(dir, "/GP_eqbyeq/")
dir.create(dir, showWarnings = FALSE) 

###--------------------------------------------------------------------------###
###-------------------------- Eq.-by-eq estimation --------------------------###
###--------------------------------------------------------------------------###
nr <- 1
for(nr in 1:MM){

foldername <- paste0(dir, nr, ".rda")
###--------------------------------------------------------------------------###
###---------------------- Median heuristics scalings ------------------------###
###--------------------------------------------------------------------------###
if(h.sc == "med-heur"){
  sigma_sq <- sigma_sq_org
  sigma_sq[,1]  <- 5e-1*sigma_sq[,1]

  X.own   <- X[,seq(nr,K, by = MM)]
  X.other <- X[,-seq(nr,K, by = MM)]
    
  sc.own <- sc.other <- matrix(NA, N, N)
  for(nn in 1:N){
    for(kk in 1:N){
      sc.own[nn,kk] <- sqrt(sum((X.own[nn,] - X.own[kk,])^2))
      sc.other[nn,kk] <- sqrt(sum((X.other[nn,] - X.other[kk,])^2))
    }  
  }
  sc.own   <- median(sc.own[lower.tri(sc.own)])/sigma_sq[nr,1]
  sc.other <- median(sc.other[lower.tri(sc.other)])/sigma_sq[nr,1]
}else{
  sc.own <- sqrt(p)
  sc.other <- sqrt(K-p)
}
hyper.grid$sc.own <- sc.own
hyper.grid$sc.other <- sc.other
  
###--------------------------------------------------------------------------###
###----------------------- Contemporaneous variables ------------------------###
###--------------------------------------------------------------------------###
if(nr == 1) zraw <- matrix(0, N, 0) else zraw <- as.matrix(yraw[,1:(nr-1)])
Z <- as.matrix(zraw[(p+1):nrow(zraw),])

###--------------------------------------------------------------------------###
###------------------------ Estimate GP regression --------------------------###
###--------------------------------------------------------------------------###
return.obj <- gp_estim(nr = nr, Y = Y, X = X, Z = Z, nsave = nsave, nburn = nburn, nthin = nthin, sv = sv, c = c, grid = hyper.grid)
  
F.i.1 <- return.obj$F[,,1]
F.i.2 <- return.obj$F[,,2]
A0.i <- return.obj$gamma[,] * (-1)
sigma2.i <- return.obj$sigma2

h.i.own     <- return.obj$lambda[,"h.own"]
h.i.other   <- return.obj$lambda[,"h.other"]
sig.i.own   <- return.obj$lambda[,"sig.own"]
sig.i.other <- return.obj$lambda[,"sig.other"]

# Storing only relevant inverse Kernels, requires less memory
Kn.inv.store <- list()
temp.grid.own <- list("h"=unique(h.i.own),"sig"=unique(sig.i.own))
temp.grid.other <- list("h"=unique(h.i.other),"sig"=unique(sig.i.other))
Kn.i.own.inv <- array(NA, c(N, N, length(temp.grid.own$h), length(temp.grid.own$sig))) #own
Kn.i.other.inv <- array(NA, c(N, N, length(temp.grid.other$h), length(temp.grid.other$sig))) #other
  
# Pre-calculate the log Kernel so we do not have to re-calculate the kernel for every set of grid values
LogKnn.own <- LogGaussKernel(Xin = t(X[,seq(nr,K, by = MM)]),Xstar = t(X[,seq(nr,K, by = MM)]), sc = sc.own)
LogKnn.other <- LogGaussKernel(Xin = t(X[,-seq(nr,K, by = MM)]),Xstar = t(X[,-seq(nr,K, by = MM)]), sc = sc.other)

for(ii in 1:(length(temp.grid.own$sig))){
  for(jj in 1: (length(temp.grid.own$h))){
    Knn.own <- temp.grid.own$sig[[ii]] * exp(temp.grid.own$h[[jj]]*LogKnn.own) 
    Knn.own.inv <- solve(Knn.own + diag(N))
    Kn.i.own.inv[,,jj,ii] <- Knn.own.inv 
  }
}
  
for(kk in 1:(length(temp.grid.other$sig))){
  for(ll in 1:(length(temp.grid.other$h))){
    Knn.other <- temp.grid.other$sig[[kk]]*exp(temp.grid.other$h[[ll]]*LogKnn.other)
    Knn.other.inv <- solve(Knn.other + diag(N))
    Kn.i.other.inv[,,ll,kk] <- Knn.other.inv 
  }
}
  
  
str.list = list(F.i.1 = F.i.1, F.i.2 = F.i.2, sigma2.i = sigma2.i, A0.i = A0.i, h.i.own = h.i.own, h.i.other = h.i.other, sig.i.own = sig.i.own, sig.i.other = sig.i.other,Kn.i.own.inv = Kn.i.own.inv, Kn.i.other.inv = Kn.i.other.inv, sc.own = sc.own, sc.other = sc.other)
save(file = foldername, "str.list")
}

stop("Eq-by-eq estimation done!")
