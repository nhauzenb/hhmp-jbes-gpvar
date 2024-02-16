###--------------------------------------------------------------------------###
###-------------------- Conjugate GP regression function --------------------###
###---------------- Hauzenberger, Huber, Marcellino, & Petz -----------------###
###--------------- Journal of Business & Economic Statistics-----------------###
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
###----------------------------- MCMC sampler -------------------------------### 
###----------------------------- for GP-VAR ---------------------------------### 
###--------------------------------------------------------------------------###
gp_estim <- function(nr, Y, X, Z, nsave, nburn, nthin, sv, mu, c, grid){

###--------------------------------------------------------------------------###
###-------------------- Gibbs sampler preliminaries -------------------------###
###--------------------------------------------------------------------------###
ntot <- nburn+nsave*nthin
save.set <- seq(nthin, nsave*nthin, nthin) + nburn
save.ind <- 0
  
MM <- ncol(Y)
N <- nrow(Y)
K <- ncol(X)
Q <- ncol(Z)

y <- Y[,nr, drop = F] #This selects the corresponding equation 

# GP set-up
sc.own <- grid$sc.own
sc.other <- grid$sc.other

h.grid <- grid$h  
sig.grid <- grid$sig
c.own <- c[["own"]]
c.other <- c[["other"]]

###--------------------------------------------------------------------------###
###-------------------------- OLS estimates ---------------------------------###
###--------------------------------------------------------------------------###
XXinv <- solve(crossprod(X) + diag(K)*1e-1)
beta.ols <- beta.draw <-  XXinv%*%crossprod(X,y)
Em <- Em.str <-  y-X%*%beta.ols
sigma2.ols <- crossprod(y-X%*%beta.ols)/N
sigma2.ols <- as.numeric(sigma2.ols)  

#-OLS based on also adding Z 
if (Q>0){
    gamma.ols <- gamma.draw <- solve(crossprod(Z))%*%crossprod(Z, y-X%*%beta.ols)
    gamma.store <- matrix(NA, nsave, Q)
}else{
    gamma.draw <- 0
    gamma.store <- NULL
}

###--------------------------------------------------------------------------###
###--------------------- Prior and grid set-up ------------------------------###
###--------------------------------------------------------------------------###
a.prior <- N/2
b.prior.own <- a.prior/c.own
b.prior.other <- a.prior/c.other

# Kernel storage for own and other lag 
Kn.full.store <- Kn.inv.store <- Kn.chol.store  <- Kn.mult.store <- Kn.store <- array(NA, c(N, N, length(h.grid), length(sig.grid), 2))
det.store <- array(NA, c(length(h.grid), length(sig.grid), 2))
dimnames(Kn.store) <- list(seq(1, N), seq(1,N), h.grid, sig.grid, c("own", "other"))
km <- nr #Set no. of equations
  
# Pre-calculate the log Kernel so we do not have to re-calculate the kernel for every set of grid values
LogKnn.own <-   LogGKcpp(Xin = t(X[,seq(km,K, by = MM)]),  Xstar = t(X[,seq(km,K, by = MM)]),  sc = sc.own)
# LogKnn.own <-   LogGaussKernel(Xin = t(X[,seq(km,K, by = MM)]),  Xstar = t(X[,seq(km,K, by = MM)]),  sc = sc.own)
LogKnn.other <- LogGKcpp(Xin = t(X[,-seq(km,K, by = MM)]), Xstar = t(X[,-seq(km,K, by = MM)]), sc = sc.other)
# LogKnn.other <- LogGaussKernel(Xin = t(X[,-seq(km,K, by = MM)]), Xstar = t(X[,-seq(km,K, by = MM)]), sc = sc.other)
print(paste0("Number of grid operations: ",length(h.grid)*length(sig.grid)))

# Precomputes several important quantities over a grid (we can even pre-compute some of the matrix mults)
count <- 0
for (jj in 1:length(h.grid)){
  for (ii in 1:length(sig.grid)){
      
    Knn.own <- sig.grid[[ii]]*exp(h.grid[[jj]]*LogKnn.own)
    Knn.other <- sig.grid[[ii]]*exp(h.grid[[jj]]*LogKnn.other)
      
    Dnn.own.inv <- solve(Knn.own + diag(N))
    Dnn.other.inv <- solve(Knn.other + diag(N))
    
    det.store[jj,ii, 1] <-  determinant(Knn.own)$modulus   #We need the determinant of Knn.own + diag(N)
    det.store[jj,ii, 2] <-  determinant(Knn.other)$modulus #We need the determinant of Knn.own + diag(N)
      
    Kn.store[,,jj,ii,1] <- Knn.own
    Kn.store[,,jj,ii,2] <- Knn.other
    
    #Kn.inv.store[,,jj,ii,1] <-  (MASS::ginv(Knn.own))
    #Kn.inv.store[,,jj,ii,2] <- (MASS::ginv(Knn.other))
    Knn.own.inv <- try(solve(Knn.own + diag(N)*1e-10), silent = FALSE); if(is(Knn.own.inv, "try-error")) Knn.own.inv <- (MASS::ginv(Knn.own))
    Knn.other.inv <- try(solve(Knn.other+ diag(N)*1e-10), silent = FALSE); if(is(Knn.other.inv, "try-error")) Knn.other.inv <- (MASS::ginv(Knn.other))
  
    Kn.inv.store[,,jj,ii,1] <- Knn.own.inv
    Kn.inv.store[,,jj,ii,2] <- Knn.other.inv
      
    Kn.mult.store[,,jj,ii,1] <- Knn.own %*% Dnn.own.inv
    Kn.mult.store[,,jj,ii,2] <- Knn.other %*% Dnn.other.inv
      
    Kn.chol.store[,,jj,ii,1] <- t(chol(Knn.own - Knn.own %*% Dnn.own.inv %*% Knn.own + diag(N)*1e-10))
    Kn.chol.store[,,jj,ii,2] <- t(chol(Knn.other - Knn.other %*% Dnn.other.inv %*% Knn.other + diag(N)*1e-10))
      
    # Include 
    Kn.full.store[,,jj,ii,1] <- Knn.own - Knn.own %*% Dnn.own.inv %*% Knn.own 
    Kn.full.store[,,jj,ii,2] <- Knn.other - Knn.other %*% Dnn.other.inv %*% Knn.other 
    
    count<- count+1
    if(count %% 50 == 0 || count == length(h.grid)*length(sig.grid)) print(paste0("Finished kernel operations for ",count," grid combinations."))
      
    }
}

grid.full <- expand.grid(h.grid, sig.grid)
Kinv.list.own <- Kinv.list.other <- list()

for (jj in 1:nrow(grid.full)){
  grid.sl <- grid.full[jj,]
  sl.h <- which(grid.sl[[1]] == h.grid)
  sl.sig  <- which(grid.sl[[2]] == sig.grid)
    
  Kinv.list.own[[jj]] <- Kn.inv.store[,,sl.h, sl.sig, 1]
  Kinv.list.other[[jj]] <- Kn.inv.store[,,sl.h, sl.sig, 2]
}
  
# Prior value for homoskedastic 
t0 <- 0.1
S0 <- 0.1
## Priors for SV
if(sv == "SV"){
  # Priors
  svprior <- list(
    "steq"  = "AR1",                                # RW vs. AR(1) state equation
    "rho"  = c("rho_m"  = 0.83, "rho_V"  = 0.0045), # Beta prior on AR(1) parameter, specified in terms of mean and variance B(a0, a1)
    "sig2" = c("sig2_m" = 0.1, "sig2_V"  = 1e-2),   # Inverse Gamma prior on the state innovation variance, specified in terms of mean and variance of IG(a0,a1)
    "h0"   = c("h0_m"   = log(sigma2.ols),   "h0_V"   = 1) # Initial state
  )
  # Newton-Raphson and MH set-up
  KIinv <- diag(1,N)
  
  svmhnr <- list("KIinv"  = KIinv, # Correlation structure between shocks defined by Kernel(s)
                 "maxitr" = 100,   # Maximum no. of iterations for Newton-Raphson
                 "tol"    = 1e-4,  # Tolerance level of Newton-Raphson  
                 "fast"   = FALSE) # Sparsify the Hessian of the likelihood to feature the same structure as the Hessian of the prior matrix
}else if(sv == "SVapprx"){
  sv_priors <- specify_priors(
    mu = sv_normal(mean = 0, sd = 1), # prior on unconditional mean in the state equation
    phi = sv_beta(shape1 = 25, shape2 = 1.5), #informative prior to push the model towards a random walk in the state equation (for comparability)
    sigma2 = sv_inverse_gamma(shape = 22, scale = 2.1), # Gamma prior on the state innovation variance
    nu = sv_infinity(),
    rho = sv_constant(0))
    # Initialization of SV processes
    svdraw <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
}

# Initialization of volatility processes
sigma2.draw <- rep(sigma2.ols, N)*1
ht <- log(sigma2.draw)
svpara <- c("h_rho"=0.95,"h_sig2"=0.01, "h0" = log(sigma2.ols))
I_T <- diag(N)

f.own <- rnorm(N,0,1)
f.other <- rnorm(N,0,1)

###--------------------------------------------------------------------------###
###---------------------------- Store objects -------------------------------###
###--------------------------------------------------------------------------###
f.store <- array(NA, c(nsave,N, 2)); dimnames(f.store) <- list(NULL, NULL, c("own", "other"))
ht.store <- sigma2.store <- array(NA, c(nsave, N))
lambda.store <- array(NA, c(nsave, 4)); dimnames(lambda.store) <- list(NULL, c("h.own", "sig.own", "h.other","sig.other"))
  
#Starting values for the Kernels
sl.h.own <- 1 
sl.sig.own <- 1
sl.h.other<- 1 
sl.sig.other <- 1
sl.sig.mu <- 10
  
if (Q>0){
    tau.hs <- 1
    zeta.hs <- 1
    nu.hs <- rep(1, Q)
    lam.hs <- rep(1, Q)
    psi.hs <- rep(1, Q)
}
  
nr_count <- svmh_acc <- 0
pb <- txtProgressBar(min = 0, max = ntot, style = 3) #start progress bar
print("Starting MCMC:")
irep <- 1
  
for (irep in seq_len(ntot)){

###------------------------------ Step 1: -----------------------------------###
###------------------- Sample contemp. relationships ------------------------###
###--------------------- and prior hyperparameters --------------------------###
###--------------------------------------------------------------------------###

###---------------------- Step 1.1: Sample covariances ----------------------###
if (Q>0){
  Z.norm <- Z * 1/sqrt(sigma2.draw)
  y.norm <- (y - f.own - f.other) * 1/sqrt(sigma2.draw)
  
if (Q>1){
    V.Z <- solve(crossprod(Z.norm) + diag(1/psi.hs))
}else if (Q==1){
      V.Z <- solve(crossprod(Z.norm) + 1/psi.hs)
}

m.Z <- V.Z %*% crossprod(Z.norm, y.norm)
gamma.draw <- m.Z + t(chol(V.Z))%*%rnorm(Q)

###------------------- Step 1.2: Sample HS hyperparameters ------------------###
hs_draw <- get.hs(bdraw=gamma.draw,lambda.hs=lam.hs,nu.hs=nu.hs,tau.hs=tau.hs,zeta.hs=zeta.hs)
lam.hs <- hs_draw$lambda
nu.hs <- hs_draw$nu
tau.hs <- hs_draw$tau
tau.hs[tau.hs > 10] <- 10
zeta.hs <- hs_draw$zeta
psi.hs <- hs_draw$psi
psi.hs[psi.hs > 100] <- 100
}
 
###------------------------------ Step 2: -----------------------------------###
###--------------- Sample stochastic volatility parameters ------------------###
###--------------------------------------------------------------------------###
if (Q > 0) y.tilde <- y - Z%*%gamma.draw else y.tilde <- y

if(sv == "SV"){
# Define KIinv for SV process
  KIn <- Kn.store[,,sl.h.own, sl.sig.own, 1]+Kn.store[,,sl.h.other, sl.sig.other, 2] + I_T
  KIn <- as(KIn, "TsparseMatrix")
  svmhnr$KIinv <- solve(KIn)
  svdraw <-  try(get.logvolas(eps=y.tilde, ht = ht, svpara = svpara,svmhnr = svmhnr), silent = FALSE)
    
  if(is(svdraw, "try-error")){
    Kn.chol.temp <- t(chol(Kn.store[,,sl.h.own, sl.sig.own, 1]+Kn.store[,,sl.h.other, sl.sig.other, 2]+diag(N)))
    Kn.vec.in <- as.numeric(sqrt(rowSums(Kn.chol.temp^2))) 
    y.dash <- y.tilde/Kn.vec.in 
      
    t1 <- t0 + N/2 # Posterior degree of freedoms
    S1 <- S0 + as.numeric(crossprod(y.dash))/2 # Posterior scaling
    sigma2.cons <- 1/rgamma(1, t1, S1) # Get variance
    sigma2.draw <- rep(sigma2.cons, N) 
    ht <- log(sigma2.draw) # Define log-variance
  }else{
    ht <- as.numeric(svdraw$ht)
    sigma2.draw <- exp(as.numeric(ht))  
    svmh_acc <- svmh_acc + svdraw$accept
    nr_count <- svdraw$nr_count
  }
  svpara <-  get.svpara(ht = ht, svpara = svpara, svprior = svprior)

}else if(sv == "SVapprx"){
  Kn.chol.temp <- t(chol(Kn.store[,,sl.h.own, sl.sig.own, 1]+Kn.store[,,sl.h.other, sl.sig.other, 2]+diag(N)))
  Kn.vec.in <- as.numeric(sqrt(rowSums(Kn.chol.temp^2)))
  y.dash <- y.tilde/Kn.vec.in 
  
  svdraw <- svsample_fast_cpp(y.dash, startpara = svdraw, startlatent = ht, priorspec = sv_priors)
  svdraw[c("mu", "phi", "sigma", "nu", "rho")] <- as.list(svdraw$para[, c("mu", "phi", "sigma", "nu", "rho")])
  ht <- t(svdraw$latent)
  ht[ht < -12] <- -12 
  sv.para <- as.numeric(svdraw$para[,c("mu", "phi", "sigma", "nu")])
      
  sigma2.draw <- exp(as.numeric(ht))  
    
}else if(sv == "homo"){
  Kn.chol.temp <- t(chol(Kn.store[,,sl.h.own, sl.sig.own, 1]+Kn.store[,,sl.h.other, sl.sig.other, 2]+diag(N)))
  Kn.vec.in <- as.numeric(sqrt(rowSums(Kn.chol.temp^2))) 
  y.dash <- y.tilde/Kn.vec.in 
    
  t1 <- t0 + N/2 # Posterior degree of freedoms
  S1 <- S0 + as.numeric(crossprod(y.dash))/2 # Posterior scaling
  sigma2.cons <- 1/rgamma(1, t1, S1) # Get variance
  sigma2.draw <- rep(sigma2.cons, N) 
  ht <- log(sigma2.draw) # Define log-variance
}
  
###------------------------------ Step 3: -----------------------------------###
###--- Sample the conditional fit related to the own lags and other lags ----###
###--------------------------------------------------------------------------###
#Step 3.1: Sample the  factors for the own lags
#Fast step which avoids matrix multis
Kn.sl.non <- Kn.store[,,sl.h.own, sl.sig.own, 1]
fhat <- diag(sqrt(sigma2.draw))%*%Kn.mult.store[,,sl.h.own,sl.sig.own,1]%*%diag(sqrt(1/sigma2.draw))%*%(y.tilde - f.other)
f.own <- fhat + (Kn.chol.store[,,sl.h.own, sl.sig.own,1]*sqrt(sigma2.draw))%*%rnorm(N)

#Step 3.2 Sample the  factors for the other lags
Kn.sl.non <- Kn.store[,,sl.h.other, sl.sig.other, 2]
fhat <- diag(sqrt(sigma2.draw))%*%Kn.mult.store[,,sl.h.other,sl.sig.other,2]%*%diag(sqrt(1/sigma2.draw))%*%(y.tilde - f.own)
f.other <- fhat + (Kn.chol.store[,,sl.h.other, sl.sig.other,2]*sqrt(sigma2.draw))%*%rnorm(N)

# Correction step for f.other: 
Kn.corr <- diag(sqrt(sigma2.draw))%*%Kn.full.store[,,sl.h.other, sl.sig.other,2]%*%diag(sqrt(sigma2.draw))
f.other.corr <- rowSums(Kn.corr)*sum(f.other)/sum(Kn.corr)
f.other <- f.other - f.other.corr
  
###------------------------------ Step 4: -----------------------------------###
###--------Sample the hyperparameters associated with the kernels -----------###
###--------------------------------------------------------------------------###

post.full <- matrix(NA, nrow(grid.full),2)
fhat.own <- f.own/sqrt(sigma2.draw)
fhat.other <- f.other/sqrt(sigma2.draw)
    
for (jj in 1:nrow(grid.full)){
  grid.sl <- grid.full[jj,]
  sl.h <- which(grid.sl[[1]] == h.grid)
  sl.sig  <- which(grid.sl[[2]] == sig.grid)
      
  K.own.inv <- Kinv.list.own[[jj]]
  K.other.inv <- Kinv.list.other[[jj]]
      
  post.full[jj,1] <- get.prior(fhat.own,   det.store[sl.h,sl.sig,1], K.own.inv,   N) + dgamma(grid.sl[[1]], a.prior, b.prior.own,   log=TRUE) + dgamma(grid.sl[[2]], a.prior, b.prior.own,   log=TRUE)
  post.full[jj,2] <- get.prior(fhat.other, det.store[sl.h,sl.sig,2], K.other.inv, N) + dgamma(grid.sl[[1]], a.prior, b.prior.other, log=TRUE) + dgamma(grid.sl[[2]], a.prior, b.prior.other, log=TRUE)
}

weights <- matrix(NA, nrow(grid.full),2)
for (ii in 1:2) weights[,ii] <- exp(post.full[,ii] - max(post.full[,ii]))/sum(exp(post.full[,ii] - max(post.full[,ii])))
  
    
sl.own <- sample(seq(1, nrow(grid.full)), 1, prob = weights[,1])
grid.own <- grid.full[sl.own,]
sl.h.own <- which(grid.own[[1]] == h.grid)
sl.sig.own  <- which(grid.own[[2]] == sig.grid)
    
sl.other <- sample(seq(1, nrow(grid.full)), 1, prob = weights[,2])
grid.other <- grid.full[sl.other,]
sl.h.other <- which(grid.other[[1]] == h.grid)
sl.sig.other  <- which(grid.other[[2]] == sig.grid)
 
###------------------------------ Step 5: -----------------------------------###
###----------------------------- Storage ------------------------------------###
###--------------------------------------------------------------------------###
if(irep %in% save.set){
  save.ind <- save.ind + 1
  f.store[save.ind,,1] <- f.own
  f.store[save.ind,,2] <- f.other
  sigma2.store[save.ind,] <- sigma2.draw
  ht.store[save.ind,] <- as.numeric(ht)
  lambda.store[save.ind,] <- c(h.grid[[sl.h.own]], sig.grid[[sl.sig.own]], h.grid[[sl.h.other]], sig.grid[[sl.sig.other]])
  
  if (Q>0) gamma.store[save.ind,] <- as.numeric(gamma.draw)
}
setTxtProgressBar(pb, irep)
}
  

ret.obj <- list("F" = f.store, "sigma2"= sigma2.store, "ht" = ht.store, "lambda"=lambda.store, "gamma"=gamma.store)
  
print(paste0("Eq ",nr,"/",MM," completed."))
return (ret.obj)
}


###--------------------------------------------------------------------------###
###-------------------------- Auxiliary functions ---------------------------### 
###----------------------------- for GP-VAR ---------------------------------### 
###--------------------------------------------------------------------------###

GaussKernel <- function(Xin,Xstar,h,sigf,sc){
  #print(1/sc)
  nstar <- ncol(Xstar)
  nin <- ncol(Xin)
  p <- nrow(Xstar)
  K <- matrix(0, nstar, nin) #This is K(X*, X)
  for (ii in 1:nstar){
    K[ii, ] <- sigf*exp((-h/(2*sc))*apply(matrix(abs((Xstar[,ii]-Xin[,])^2), p, nin),2, sum))
  }
  return(K)
}

LogGaussKernel <- function(Xin,Xstar,sc){ #,sigf
  print(1/sc)
  nstar <- ncol(Xstar)
  nin <- ncol(Xin)
  p <- nrow(Xstar)
  K <- matrix(0, nstar, nin) #This is K(X*, X)
  for (ii in 1:nstar){
    K[ii, ] <- (-1/(2*sc))*apply(matrix(abs((Xstar[,ii]-Xin[,])^2), p, nin),2, sum) # + log(sigf) 
  }
  return(K)
}

# Get conditional prediction with GP without average variance
get.GPpred <- function(Kstar, KKn.out, Kn.inv, cfe, sigma2, tt.pos){
  sigma2.t <- mean(sigma2)
  KKnKninv <- crossprod(t(KKn.out), Kn.inv)
  fhat.out <- as.numeric(KKnKninv%*%cfe) # fhat.out <- as.numeric(sqrt(sigma2.t)*KKn.out %*% Kn.inv%*%diag(sqrt(1/sigma2))%*%cfe)
  Vhat.out <- as.numeric(Kstar - KKnKninv%*%t(KKn.out))*sigma2.t  # Vhat.out <- as.numeric(Kstar - KKn.out %*%Kn.inv%*% t(KKn.out))*sigma2.t
  out.draw <- rnorm(1, fhat.out, sqrt(Vhat.out))
  return(out.draw)
}

#Function creating the lags of y (alternative to embed)
mlag <- function(X,lag)
{
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)  
}

# Get posteriors for the horseshoe prior (see Makalic & Schmidt, 2015)
get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1   
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
  }
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

###############################################################################
###-------------------------- Auxiliary functions --------------------------### 
###---------------------- for stochastic volatility ------------------------### 
###############################################################################
get.prior <- function(f.draw, det.sig, sigma.inv, N){
  -1/2*N*log(2*pi) -1/2 *  (det.sig) - 1/2 * sum(crossprod(f.draw,sigma.inv)*as.numeric(f.draw))
}

log.lik.func <- function(eps, htt, KIinv){
  omega    <- exp(htt)
  eps.norm <- eps/sqrt(omega)
  eval.lik <- -1/2*sum(htt) - 1/2*t(eps.norm)%*%KIinv%*%eps.norm
  return(eval.lik)  
}

log.pr.func  <- function(htt, HH, h0, h_sig2){
  htt.norm <- htt - h0
  eval.pr <-  -1/(2*h_sig2)*t(htt.norm)%*%HH%*%htt.norm   
  return(eval.pr)  
}

prop.corr.func  <- function(htt, prop_mean, prop_prec){
  htt.norm <- htt - prop_mean
  eval.prop <-  -1/2*t(htt.norm)%*%prop_prec%*%htt.norm   
  return(eval.prop)  
}

get.svpara <- function(ht, svpara, svprior){
  # ----------------------------------------------------------------------------
  # Extract previous draw
  # ----------------------------------------------------------------------------
  ht <- as.numeric(ht)
  T  <- length(ht)
  h_mu   <- svpara["h_mu"]
  h_rho  <- svpara["h_rho"]
  h_sig2 <- svpara["h_sig2"]
  h0  <- svpara["h0"]
  
  # ----------------------------------------------------------------------------
  # Draw parameters of the state equation
  # ----------------------------------------------------------------------------
  steq <- svprior$steq
  
  # 1.) Sample initial state h0
  h0_pr <- svprior$h0
  if(abs(h_rho) < 0.95){
    h0_m <- 0
    h0_V <- h_sig2/(1 - h_rho^2)
  }else{
    h0_m <-  h0_pr["h0_m"]
    h0_V  <- h0_pr["h0_V"]
  }
  
  h0_Vp <- 1/(1/h0_V+h_rho^2/h_sig2)
  h0_mp <- h0_Vp*(h0_m/h0_V + (h_rho*ht[1])/h_sig2)
  
  h0 <- rnorm(1,h0_mp,sqrt(h0_Vp))
  
  # 2.) Sample AR(1) parameter 
  # Prior: State innovation variance
  sig2_pr  <- svprior$sig2
  sig2_m   <- sig2_pr["sig2_m"]
  sig2_V   <- sig2_pr["sig2_V"]
  
  sig2_a0 <- (sig2_m^2/sig2_V)+2
  sig2_b0 <- sig2_m*(sig2_a0-1)
  
  # Design matrices
  ht_lag <- c(h0, ht[1:(T-1)])
  
  if(steq == "RW"){
    h_rho <- 1         
  }else if(steq == "AR1"){
    # Prior: AR(1) parameter
    rho_pr  <- svprior$rho
    rho_m   <- rho_pr["rho_m"]
    rho_V   <- rho_pr["rho_V"]
    
    # Non-conjugate with Beta prior on AR1 coefficient
    # Update AR1 parameter (slice sampler)
    rho_m <- (rho_m + 1)/2
    rho_a <- ((1-rho_m)*rho_m- rho_V)*rho_m/rho_V
    rho_b <- rho_a*(1-rho_m)/rho_m
    
    rho_tr <- (h_rho + 1)/2 # Transform h_rho to be bounded between 0 and 1
    rho_tr <- uni.slice(rho_tr, g = function(x){
      -0.5/h_sig2*sum((ht - (2*x - 1)*ht_lag)^2) +
        dbeta(x, shape1 = rho_a, shape2 = rho_b, log = TRUE)
    }, lower = 0.05, upper = 0.99)[1]
    
    h_rho <- 2*rho_tr - 1 # Map transformed h_rho back
  }
  
  # 3.) Sample state innovation variance 
  sig2_a1 <- sig2_a0 + T/2
  
  ssr <- sum((ht - h_rho*ht_lag)^2)
  sig2_b1 <- sig2_b0 + ssr/2
  
  h_sig2 <- 1/rgamma(1,sig2_a1,sig2_b1)
  
  svpara <- c("h_rho"=h_rho,"h_sig2"=h_sig2, "h0" = h0)
  return("svpara"=svpara)
}

# Chan et al algorithm with independent MH and Newton Raphson: 
get.logvolas <- function(eps=yt, ht = ht, svpara = svpara, svmhnr = svmhnr){
  # eps=yt; ht = ht; svpara = svpara; svmhnr = svmhnr
  # ----------------------------------------------------------------------------
  # Extract previous draw
  # ----------------------------------------------------------------------------
  ht <- htt_acc <- as.numeric(ht)
  T  <- length(ht)
  h_rho  <- svpara["h_rho"]
  h_sig2 <- svpara["h_sig2"]
  h0  <- svpara["h0"]
  h0c.vec <- c(h_rho*h0, rep(0, T-1)) 
  
  # ----------------------------------------------------------------------------
  # Draw log-volatilities
  # ----------------------------------------------------------------------------
  KIinv <- svmhnr$KIinv
  nr_maxitr <- svmhnr$maxitr
  nr_tol    <- svmhnr$tol
  mhnr_fast   <- svmhnr$fast
  
  htt <- ht # Define auxiliary ht for NR and MH step
  # Pre-computation of prior quantities
  R <- diag(1, T)
  diag(R[2:T,1:(T-1)]) <- -h_rho
  R <- as(R, "TsparseMatrix")
  #R <- matrix(0, T, T)
  #R[c(T*(0:(T-1)) + 1:T)] <- 1
  #R[c(T*(0:(T-2)) + 2:T)] <- -h_rho
  delta <- solve(R, h0c.vec) 
  
  HH <- crossprod(R) #Matrix::crossprod(R)
  HH <- as(HH, "TsparseMatrix")
  nr_crit <- 100; nr_count <- 0
  sps.id <- as.vector(HH == 0)
  
  hess.pr <- -1/h_sig2*HH # Negative Hessian of  log prior density
  start <- Sys.time()
  while(nr_crit > nr_tol & nr_count < nr_maxitr){
    omega <- exp(htt)
    eps.norm <- as.numeric(eps/sqrt(omega))
    
    # Prior moments for proposal
    grad.pr <- as.numeric(hess.pr%*%(htt - delta)) # Gradient of log prior density
    KIeps <- KIinv%*%eps.norm
    KIeps <- as.numeric(KIeps*eps.norm)
    
    # Likelihood moments for proposal
    grad.lik <- 1/2*(KIeps -1) # Gradient of log likelihood 
    # Negative Hessian  of log likelihood (product rule)
    hess.lik <- as(-1/4*t((KIinv*eps.norm))*eps.norm, "TsparseMatrix")
    sl.diag <- hess.lik@i==hess.lik@j
    hess.lik@x[sl.diag] <-   hess.lik@x[sl.diag] -1/4*KIeps
    
    if(mhnr_fast & !Matrix::isDiagonal(hess.lik)){
      #hess.lik1 <- Matrix::Diagonal(eps.norm,n = T)%*%KIinv%*%Matrix::Diagonal(eps.norm,n = T)
      #hess.lik1 <- t((KIinv*eps.norm))*eps.norm # KIinv*eps.norm^2#diag(eps.norm)%*%KIinv%*%diag(eps.norm) #CAREFUL HERE: what happens if KIinv is non-diagonal
      #hess.lik2 <- Matrix::Diagonal(KIeps,n = T)
      #hess.lik.backup <-  -1/4*(hess.lik1 + hess.lik2)
      #hess.lik.backup[sps.id] <- 0
      sl.Udiag <- hess.lik@i == (hess.lik@j +1)
      sl.Ldiag <- hess.lik@i==  (hess.lik@j -1)
      sl.LUdiag <- (sl.diag + sl.Udiag + sl.Ldiag)
      hess.lik@x[!sl.LUdiag] <- 0
    }
    grad.post <- grad.pr + grad.lik
    hess.post <- hess.pr + hess.lik
    htt_new   <- htt - Matrix::solve(hess.post, grad.post)
    
    # Check convergences of Newton-Raphson algorithm
    nr_crit <- max(abs(htt_new - htt))
    nr_count <- nr_count + 1  
    
    htt <- htt_new
  }
  end <- Sys.time()
  
  omega <- exp(htt)
  eps.norm <- as.numeric(eps/sqrt(omega))
  
  hess.lik <- as(-1/4*t((KIinv*eps.norm))*eps.norm, "TsparseMatrix")
  sl.diag <- hess.lik@i==hess.lik@j
  hess.lik@x[sl.diag] <-   hess.lik@x[sl.diag] -1/4*KIeps
  
  if(mhnr_fast & !Matrix::isDiagonal(hess.lik)){
    sl.Udiag <- hess.lik@i == (hess.lik@j +1)
    sl.Ldiag <- hess.lik@i==  (hess.lik@j -1)
    sl.LUdiag <- (sl.diag + sl.Udiag + sl.Ldiag)
    hess.lik@x[!sl.LUdiag] <- 0
  }
  
  prop_prec <- (-1)*(hess.pr + hess.lik) # Negative Hessian  of the posterior
  prop_var <- solve(prop_prec)  # Proposal variance
  prop_mean <- htt # Proposal mean
  htt_prop <- prop_mean + t(chol(as.matrix(prop_var)))%*%rnorm(T)
  
  alp_prop <- log.lik.func(eps = eps, htt = htt_prop, KIinv = KIinv) + log.pr.func(htt = htt_prop, HH = HH, h0 = h0, h_sig2 = h_sig2) + prop.corr.func(htt = htt_acc , prop_mean = prop_mean, prop_prec = prop_prec)
  alp_acc  <- log.lik.func(eps = eps, htt = htt_acc , KIinv = KIinv) + log.pr.func(htt = htt_acc , HH = HH, h0 = h0, h_sig2 = h_sig2) + prop.corr.func(htt = htt_prop, prop_mean = prop_mean, prop_prec = prop_prec) 
  
  alp_MH <- as.numeric(alp_prop - alp_acc)
  
  if(alp_MH > log(runif(1))){
    ht <- matrix(as.numeric(htt_prop), T, 1)
    accept <- TRUE
  }else{
    ht <- htt_acc
    accept <- FALSE
  }
  return(list("ht"=ht, "nr_count" = nr_count, "accept" = accept))
  
}

# Univariate Slice Sampler from Neal (2008)
uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL)
{
  # Check the validity of the arguments.
  if (!is.numeric(x0) || length(x0)!=1
      || !is.function(g)
      || !is.numeric(w) || length(w)!=1 || w<=0
      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
      || !is.numeric(lower) || length(lower)!=1 || x0<lower
      || !is.numeric(upper) || length(upper)!=1 || x0>upper
      || upper<=lower
      || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
  {
    stop ("Invalid slice sampling argument")
  }
  
  # Find the log density at the initial point, if not already known.
  
  if (is.null(gx0))
  { #uni.slice.evals <<- uni.slice.evals + 1
    gx0 <- g(x0)
  }
  
  # Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)
  
  # Find the initial interval to sample from.
  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff
  
  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.
  
  if (is.infinite(m))  # no limit on number of steps
  {
    repeat
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
    }
    
    repeat
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
    }
  }
  
  else if (m>1)  # limit on steps, bigger than one
  {
    J <- floor(runif(1,0,m))
    K <- (m-1) - J
    
    while (J>0)
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }
    
    while (K>0)
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }
  
  # Shrink interval to lower and upper bounds.
  
  if (L<lower)
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }
  
  # Sample from the interval, shrinking it on each rejection.
  
  repeat
  {
    x1 <- runif(1,L,R)
    
    #uni.slice.evals <<- uni.slice.evals + 1
    gx1 <- g(x1)
    
    if (gx1>=logy) break
    
    if (x1>x0)
    { R <- x1
    }
    else
    { L <- x1
    }
  }
  
  # Return the point sampled, with its log density attached as an attribute.
  attr(x1,"log.density") <- gx1
  return (x1)
}
