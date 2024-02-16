###--------------------------------------------------------------------------###
###-- Gaussian process vector autoregressions and macroeconomic uncertainty -###
###---------------- Hauzenberger, Huber, Marcellino, & Petz -----------------###
###--------------- Journal of Business & Economic Statistics-----------------###
###--------------------------------------------------------------------------###
###-------------------------- GP-VAR code files -----------------------------###
###------------------------------- GIRFs ------------------------------------###
###--------------------------------------------------------------------------###
rm(list = ls())

w.dir <- ""
est.dir   <-  paste0(w.dir, "eqbyeq/"); dir.create(est.dir, showWarnings = FALSE)
girf.dir  <-  paste0(w.dir, "girfs/") ; dir.create(girf.dir, showWarnings = FALSE)
func.dir  <-  paste0(w.dir, "gpvar_funcs/")

library(MASS)
library(Matrix)
library(mvtnorm)
library(truncnorm)
library(crayon)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
library(readr)

###--------------------------------------------------------------------------###
###---------------------- Auxiliary functions -------------------------------###
###--------------------------------------------------------------------------###
source(paste0(func.dir, "gp_eqbyeq_func_svimh_main.R"))  
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

shk.sc    <- c("pos" = 1)
no.shks   <- length(shk.sc)
girf.prds <- seq(5, 120, 1) # relevant for GIRFs
N.girf    <- length(girf.prds)

girf.setup <- list(
  nsave     = 5000,
  shk.var   = "vxo",
  shk.sc    = shk.sc,
  no.shks   = no.shks,
  shk.stdz  = c("prd-wise"),  # choose one: "none", "avg", "prd-wise" 
  nhor      = 17,   
  girf.prds = girf.prds, # relevant for GIRFs
  N.girf    = N.girf
  
)

list2env(model.setup, .GlobalEnv)
list2env(girf.setup, .GlobalEnv)

###--------------------------------------------------------------------------###
###---------------------- Basu & Bundick data set-up ------------------------###
###--------------------------------------------------------------------------###
basu_bundick_data_i <- read_csv(paste0(w.dir, "data/BB_realization.csv"))
var.names <- var.set <- colnames(basu_bundick_data_i)

basu_bundick_data_i[,1] <- 100*log(basu_bundick_data_i[,1])
basu_bundick_data_i[,2] <- 100*log(basu_bundick_data_i[,2])
basu_bundick_data_i[,3] <- 100*log(basu_bundick_data_i[,3])
basu_bundick_data_i[,4] <- 100*log(basu_bundick_data_i[,4])
basu_bundick_data_i[,5] <- 100*log(basu_bundick_data_i[,5])
basu_bundick_data_i[,6] <- 100*(basu_bundick_data_i[,6]-1)
basu_bundick_data_i[,7] <- 100*log(basu_bundick_data_i[,7])

yraw <- basu_bundick_data_i
yraw <- ts(as.matrix(yraw)) # data matrix (ranges from t = 1 to t = end.smp)
MM <- ncol(yraw)
var.order <- 1:MM

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

###--------------------------------------------------------------------------###
###------------------------------- Directories ------------------------------###
###--------------------------------------------------------------------------###
dir <- paste0(w.dir, data.set, "-", unc.ind, "-p", p)
eqbyeq.dir <- paste0(dir, "/GP_eqbyeq/")
eq.files   <- list.files(eqbyeq.dir, pattern = paste0(1:MM, ".rda$", collapse = "|"))

print(eq.files)
print(eqbyeq.dir)
  
if(length(eq.files) != MM){
    miss.eqs <- setdiff(paste0(1:MM, ".rda"), eq.files)
    stop("Equation(s) are missing.")
    save(file = paste0(eqbyeq.dir, "/!missingEqs.rda"), list = c("miss.eqs"))
}
  

###--------------------------------------------------------------------------###
###--------------------------- Storage objects ------------------------------###
###--------------------------------------------------------------------------###
A0.store <- A0.inv.store <- array(0, c(nsave,MM,MM))
dimnames(A0.store) <- dimnames(A0.inv.store) <- list(c(1:nsave), var.names, var.names)
A.str.store <- A.store <- array(0,c(nsave,K, MM))
dimnames(A.str.store) <- dimnames(A.store) <- list(c(1:nsave), 1:K, var.names)
F1.str.store <- F2.str.store <- F1.store <- F2.store <- F.store <- array(0,c(nsave,N,MM))
dimnames(F.store) <- dimnames(F1.str.store) <- dimnames(F2.str.store) <- dimnames(F1.store) <- dimnames(F2.store) <- list(c(1:nsave), 1:N, var.names)
SIGMA.store <- array(0,c(nsave,N,MM,MM))
dimnames(SIGMA.store) <- list(c(1:nsave), 1:N, var.names, var.names)
sigma2.store <- array(0,c(nsave,N,MM))
dimnames(sigma2.store) <- list(c(1:nsave), 1:N, var.names)
h.store <- sig.store <- array(0,c(nsave,MM,2))
sc.store <- matrix(0,MM,2)
Kn.inv.store <- list()
temp.grid.own <- temp.grid.other <- list()
# GIRFs 
girf.store <- array(0, c(nsave, MM, nhor,no.shks))
dimnames(girf.store) <- list("draws" = 1:nsave,"variable"= var.names, "horizon"=0:(nhor-1),  "sign"= names(shk.sc))
girf_tbyt.store <- array(NA, c(nsave, MM, nhor, no.shks, N.girf))
dimnames(girf_tbyt.store) <- list("draws" = 1:nsave,"variable"= colnames(Y), "horizon"=0:(nhor-1),  "sign"= names(shk.sc), "date" = girf.prds)

###--------------------------------------------------------------------------###
###----- Load estimates and store them according to the structural VAR ------###
###--------------------------------------------------------------------------###
for(ii in 1:MM){
  load(paste0(eqbyeq.dir,"/", ii, ".rda"))
  list2env(str.list,globalenv())
  Xginv <- MASS::ginv(X)  
  
  A0.store[,ii,ii] <- 1
  if(ii > 1) A0.store[,ii,1:(ii-1)] <- A0.i
  
  F1.str.store[,,ii] <- F.i.1
  F2.str.store[,,ii] <- F.i.2
  A.str.store[,,ii] <- t(Xginv%*%t(F.i.1 + F.i.2))
  
  sigma2.store[,,ii] <- sigma2.i
  h.store[,ii,1] <- h.i.own
  h.store[,ii,2] <- h.i.other
  sig.store[,ii,1] <- sig.i.own
  sig.store[,ii,2] <- sig.i.other
  sc.store[ii,1] <- sc.own
  sc.store[ii,2] <- sc.other
  
  # Creating a new Grid with only the h and sig values that are used in estimation
  kn.h.ident.own <- unique(h.store[,ii,1])
  kn.sig.ident.own <- unique(sig.store[,ii,1])
  kn.h.ident.other <- unique(h.store[,ii,2])
  kn.sig.ident.other <- unique(sig.store[,ii,2])
  temp.grid.own[[ii]] <- list("h"=kn.h.ident.own,"sig"=kn.sig.ident.own)
  temp.grid.other[[ii]] <- list("h"=kn.h.ident.other,"sig"=kn.sig.ident.other)
    
  Kn.inv.store[[ii]] <- list()
  Kn.inv.store[[ii]][[1]]<- Kn.i.own.inv
  Kn.inv.store[[ii]][[2]]<- Kn.i.other.inv
}
 
###--------------------------------------------------------------------------###
###----------------------------- Reduced form VAR ---------------------------###
###--------------------------------------------------------------------------###
jj <- 1
for(jj in 1:nsave){
  A0.inv.temp <- try(solve(A0.store[jj,,]), silent = FALSE); if(is(A0.inv.temp, "try-error")) A0.inv.temp <- (MASS::ginv(A0.store[jj,,]))
  A0.inv.store[jj,,] <- A0.inv.temp
  A.store[jj,,] <- A.str.store[jj,,]%*%t(A0.inv.temp)
  for(tt in 1:N){
      F1.store[jj,tt,] <- F1.str.store[jj,tt,] %*% t(A0.inv.temp)
      F2.store[jj,tt,] <- F2.str.store[jj,tt,] %*% t(A0.inv.temp)
      F.store[jj,tt,] <- (F1.str.store[jj,tt,] + F2.str.store[jj,tt,])%*% t(A0.inv.temp)
      
      SIGMA.store[jj,tt,,] <- A0.inv.store[jj,,]%*%diag(sigma2.store[jj,tt,])%*%t(A0.inv.store[jj,,])
  }
}
 
###--------------------------------------------------------------------------###
###------------------  Structural VAR analysis: GIRFs -----------------------###
###--------------------------------------------------------------------------###
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

tt.pos <- which(time(yraw)[-c(1:p)] %in% (girf.prds+p)) # Get positions of periods! 
tt.pos <- girf.prds[1]
count <- 1 

for(tt.pos in girf.prds){
  
gi <- 1
for (gi in 1:nsave){
girf <- X1.girf <- list() # Storage for GIRFs
A.draw <- A.store[gi,,]

# Shock at h = 0: 
shk.tt.1 <- rnorm(MM, 0, sqrt(sigma2.store[gi,tt.pos,]))
shk.tt.0 <- rnorm(MM, 0, sqrt(sigma2.store[gi,tt.pos,]))
names(shk.tt.1) <- names(shk.tt.1) <- var.names
shk.tt.1["vxo"] <- 0
shk.tt.0[] <- shk.tt.1[] <- 0

# Unconditional forecast
if(p == 1){
    X0.fcst <- X0.girf <- matrix(Y[tt.pos,] + shk.tt.0,1,MM) 
}else{
    X0.fcst <- X0.girf <- matrix(c(Y[tt.pos,] + shk.tt.0,X[tt.pos,1:(MM*(p-1))]),1,MM*p) 
}

# Conditional forecast on Cholesky shock 
SIGMA.t <-  SIGMA.store[gi,tt.pos,,] 
shock.chol <- t(chol(SIGMA.t)) 

for(ss in 1:no.shks){
  girf[[ss]] <- matrix(0,MM,nhor)
  if(shk.stdz == "none"){
    shk.stdz.sc <- 1
  }else if(shk.stdz == "avg"){
    shk.stdz.sc <- mean(sqrt(sigma2.store[gi,,shk.var]))
  }else if(shk.stdz == "prd-wise"){
    shk.stdz.sc <- shock.chol[shk.var,shk.var]
  }
  
  girf[[ss]][,1] <- shock.chol[,shk.var]*shk.sc[[ss]]/shk.stdz.sc + shk.tt.1 
  if(p == 1){
      X1.girf[[ss]] <- matrix(girf[[ss]][,1] + Y[tt.pos,],1,MM)  
  }else{
      X1.girf[[ss]] <- matrix(c(girf[[ss]][,1], rep(0, MM*(p-1))) + c(Y[tt.pos,],X[tt.pos,1:(MM*(p-1))]),1,MM*p)
  }
}
names(girf) <- names(X1.girf) <- names(shk.sc)
  
for (hh in 2:nhor){
  fhat0.out <- matrix(0,MM,1) 
  fhat1.out <- matrix(0,MM,no.shks)
  for (mm in 1:MM){
    own.ident <- seq(mm, K, by = MM)       # identify own lags
    other.ident <- setdiff(1:K, own.ident) # identify other lags
    if(mm == 1) {
        cfe.1 <- Y[,mm] - F1.store[gi,,mm] # conditional forecast error for other lags 
        cfe.2 <- Y[,mm] - F2.store[gi,,mm] # conditional forecast error for own lags 
    }else{
        cfe.1 <- Y[, mm] - F1.store[gi,,mm] + as.numeric(Y[,1:(mm - 1)] %*% as.matrix(A0.store[gi,mm,1:(mm - 1)])) # conditional forecast error for other lags
        cfe.2 <- Y[, mm] - F2.store[gi,,mm] + as.numeric(Y[,1:(mm - 1)] %*% as.matrix(A0.store[gi,mm,1:(mm - 1)])) # conditional forecast error for own lags
    }
      
    # Identification of equation specific h and sig values
    h.own.ident <- which(h.store[gi, mm, 1] == temp.grid.own[[mm]]$h)
    sig.own.ident <- which(sig.store[gi, mm, 1] == temp.grid.own[[mm]]$sig)
    h.other.ident <- which(h.store[gi, mm, 2] == temp.grid.other[[mm]]$h)
    sig.other.ident <- which(sig.store[gi, mm, 2] == temp.grid.other[[mm]]$sig)
      
    # Median heuristics
    sc.own <- sc.store[mm,1]
    sc.other <- sc.store[mm,2]

      ###------------- GIRFs forecasts (just expectation) -------------------###
      # Prediction-step with own lags
      KKn0.out <- GKcpp(Xin = t(X[,own.ident]),Xstar = t(X0.girf[,own.ident, drop = F]),h = h.store[gi,mm,1],sigf = sig.store[gi,mm,1], sc = sc.own)
      # KKn0.out <- GaussKernel(Xin = t(X[,own.ident]),Xstar = t(X0.girf[,own.ident, drop = F]),h = h.store[gi,mm,1],sigf = sig.store[gi,mm,1], sc = sc.own)
      own0.out <- get.GPpred(Kstar = sig.store[gi,mm,1], KKn.out = KKn0.out, Kn.inv = Kn.inv.store[[mm]][[1]][,,h.own.ident, sig.own.ident], cfe = cfe.2, sigma2 = sigma2.store[gi,,mm], tt.pos = tt.pos)
      
      # Prediction-step with other lags 
      KKn0.out <- GKcpp(Xin = t(X[,other.ident]),Xstar = t(X0.girf[,other.ident, drop = F]),h = h.store[gi,mm,2],sigf = sig.store[gi,mm,2],sc = sc.other)
      # KKn0.out   <- GaussKernel(Xin = t(X[,other.ident]),Xstar = t(X0.girf[,other.ident, drop = F]),h = h.store[gi,mm,2],sigf = sig.store[gi,mm,2],sc = sc.other)
      other0.out <- get.GPpred(Kstar = sig.store[gi,mm,2], KKn.out = KKn0.out, Kn.inv = Kn.inv.store[[mm]][[2]][,,h.other.ident,sig.other.ident], cfe = cfe.1, sigma2 = sigma2.store[gi,, mm], tt.pos = tt.pos)
      
      # Add both
      fhat0.out[mm] <- other0.out + own0.out #  + rnorm(1, 0, sqrt(sigma2.store[gi,tt.pos,mm]))
    
      for(ss in 1:no.shks){
        # Prediction-step with own lags
        KKn1.out <- GKcpp(Xin = t(X[,own.ident]),Xstar = t(X1.girf[[ss]][,own.ident, drop = F]),h = h.store[gi,mm,1],sigf = sig.store[gi,mm,1], sc = sc.own)
        # KKn1.out <- GaussKernel(Xin = t(X[,own.ident]),Xstar = t(X1.girf[[ss]][,own.ident, drop = F]),h = h.store[gi,mm,1],sigf = sig.store[gi,mm,1], sc = sc.own)
        own1.out <- get.GPpred(Kstar = sig.store[gi,mm,1], KKn.out = KKn1.out, Kn.inv = Kn.inv.store[[mm]][[1]][,,h.own.ident,sig.own.ident], cfe = cfe.2, sigma2 = sigma2.store[gi,,mm], tt.pos = tt.pos)
        # Prediction-step with other lags
        KKn1.out <- GKcpp(Xin = t(X[,other.ident]),Xstar = t(X1.girf[[ss]][,other.ident, drop = F]),h = h.store[gi,mm,2],sigf = sig.store[gi,mm,2], sc = sc.other)
        # KKn1.out <- GaussKernel(Xin = t(X[,other.ident]),Xstar = t(X1.girf[[ss]][,other.ident, drop = F]),h = h.store[gi,mm,2],sigf = sig.store[gi,mm,2], sc = sc.other)
        other1.out <- get.GPpred(Kstar = sig.store[gi, mm, 2], KKn.out = KKn1.out, Kn.inv = Kn.inv.store[[mm]][[2]][,,h.other.ident,sig.other.ident], cfe = cfe.1, sigma2 = sigma2.store[gi,, mm], tt.pos = tt.pos)
        # Add both 
        fhat1.out[mm,ss] <- other1.out + own1.out # + rnorm(1, 0, sqrt(sigma2.store[gi,tt.pos,mm]))
      }
  }
  
  ###-------------------- Map structural into reduced form ------------------###
  fhat0.out <- A0.inv.store[gi, , ] %*% fhat0.out
  X0.girf <- matrix(c(fhat0.out, X0.girf[1:(MM * (p - 1))]), 1, K) # Create new X.out for (h+1)-step-ahead
  for(ss in 1:no.shks){
    fhat1.out[,ss] <- A0.inv.store[gi, , ] %*% fhat1.out[,ss]
    X1.girf[[ss]] <- matrix(c(fhat1.out[,ss], X1.girf[[ss]][1:(MM * (p - 1))]), 1, K) # Create new X.out for (h+1)-step-ahead
    girf[[ss]][, hh] <- (fhat1.out[,ss] - fhat0.out) # store GIRFs
  }
  
}
for(ss in 1:no.shks) girf.store[gi, , , ss] <- girf[[ss]]*matrix(sc.sd, MM, nhor) 
print(gi)
}
print(paste0("Period: ",tt.pos," completed."))
girf_tbyt.store[,,,,count]  <- girf.store
count <- count + 1

} 
 
###--------------------------------------------------------------------------###
###------------------- Posterior summary quantities -------------------------###
###--------------------------------------------------------------------------###
girf_avg.store <- apply(girf_tbyt.store, c(1,2,3,4), mean)
girf_tbyt_post <- apply(girf_tbyt.store, c(2,3,4,5), quantile, c(0.16,0.25, 0.5, 0.75, 0.84), na.rm = T)
girf_avg_post <- apply(girf_avg.store, c(2,3,4), quantile, c(0.16,0.25, 0.5, 0.75, 0.84), na.rm = T)
dimnames(girf_tbyt_post) <- list("moment" = c("llow", "low","med","high", "hhigh"),"variable"= colnames(Y), "horizon"=0:(nhor-1),  "sign"= names(shk.sc), "date" = girf.prds)
dimnames(girf_avg_post) <- list("moment" = c("llow", "low","med","high", "hhigh"),"variable"= colnames(Y), "horizon"=0:(nhor-1),  "sign"= names(shk.sc))
girf.list <- list(girf_tbyt_post = girf_tbyt_post, girf_avg_post = girf_avg_post)

girf.dir <- paste0(dir, "/GIRFs/")
dir.create(girf.dir, showWarnings = FALSE) 
save(file = paste0(girf.dir, "/girf_", shk.var, "_", shk.stdz, "collected.rda"), list = c("girf.list"))