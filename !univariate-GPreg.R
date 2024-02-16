###--------------------------------------------------------------------------###
###-- Gaussian process vector autoregressions and macroeconomic uncertainty -###
###---------------- Hauzenberger, Huber, Marcellino, & Petz -----------------###
###--------------- Journal of Business & Economic Statistics-----------------###
###--------------------------------------------------------------------------###
###------------ Simple univariate example with a GP regression --------------###
###--------------------------------------------------------------------------###
rm(list = ls())
w.dir <- ""
library(ggplot2)
library(stochvol)
library(Rcpp)
library(latex2exp)
library(cowplot)
library(zoo)
library(scales)

figs.dir  <-  paste0(w.dir, "figs/"); dir.create(figs.dir, showWarnings = FALSE)
func.dir  <-  paste0(w.dir, "gpvar_funcs/")

###--------------------------------------------------------------------------###
###-----------------------Auxiliary functions -------------------------------###
###--------------------------------------------------------------------------###
sourceCpp(paste0(func.dir, "sqexp_kernel.cpp"))

###--------------------------------------------------------------------------###
###--------------------------------- Set-up ---------------------------------###
###--------------------------------------------------------------------------###
set.seed(2000)          # Set seed
start.date <-  2005     # Specify start of sample
end.date   <-  2014+3/4 # Specify end of sample
ntot       <-  1000     # No. of draws
# Varying the kernels
kernel.h <- c("linear"     = 1, 
              "tight"      = 0.01, 
              "loose"      = 0.1,
              "very loose" = 4) 
no.kns <- length(kernel.h) # No. of kernels
kernel.sig <- 5            # Kernel scaling

###--------------------------------------------------------------------------###
###-------------------------- Plot preliminaries ----------------------------###
###--------------------------------------------------------------------------###
color.est <- "darkred"
color.obs <- "black"
no.drws <- 3 # No. of draws from the prior shown in the graph
ls.est <- 1
ls.draw <- 0.5
ps.obs <- 1
ylim_prior <- c(-6, 10)
ylim_post  <- c(-6, 10)  

###--------------------------------------------------------------------------###
###------------------------------ Data set-up -------------------------------###
###--------------------------------------------------------------------------###
load("data//MacroUncQ.rda") # Load data
# Define explanatory variable: lag of macro uncertainty
Xraw <- yraw[,"JLN2015_ORG"]
Xraw <- window(Xraw, start = start.date-0.25, end = end.date-0.25) 
mean.x <- mean(Xraw)
X <- as.matrix(Xraw-mean.x)
# Define dependent variable: real GDP
yraw <- yraw[,"GDPC1",drop = F] 
yraw <- window(yraw, start = start.date, end = end.date) 
mean.y <- mean(yraw)
y <- as.matrix(yraw-mean.y)
# Get dimensions
K  <- ncol(X)
MM <- ncol(y)
N  <- nrow(y)

# Variance starting values
sigma2 <- 1e-1
sigma2.draw <- rep(sigma2, N)

###--------------------------------------------------------------------------###
###----------------------------- Prior set-up -------------------------------###
###--------------------------------------------------------------------------###
run <- 1
ggp_prior <- ggp_post <- ggp_dens <- list()
for(run in 1:no.kns){
  
kns.slct <- kernel.h[run]
kns.lbl <- names(kns.slct)
# Define GP prior based on run
if(kns.lbl == "linear"){
  Kn <- cbind(X, 1)%*%t(cbind(X, 1))*kns.slct
}else{
  Kn <- GKcpp(t(X),t(X), sig = kernel.sig, h = kns.slct,sc = 1) 
}
Kn.chol <- solve(t(chol(Kn+diag(N))))

###--------------------------------------------------------------------------###
###------------------------------- Storage ----------------------------------###
###--------------------------------------------------------------------------###
f.store <- array(0,c(N,ntot))
prior.store <- array(0,c(N,ntot))
pred.store <- array(0,c(N,ntot))

###--------------------------------------------------------------------------###
###---------------------------- Start MCMC loop -----------------------------###
###--------------------------------------------------------------------------###
irep <- 1
for(irep in seq_len(ntot)){

###------------------- Step 1.: Draw from GP prior --------------------------###
prior.store[,irep] <- MASS::mvrnorm(1,rep(0,N), Kn) + mean.y

###------------------- Step 2.: Draw from GP posterior -----------------------###
Dinv <- solve(Kn + diag(sigma2.draw))
Vhat <- Kn - Kn %*% Dinv %*% Kn
fhat <- Kn %*% Dinv %*% y
  
f.draw <- try(fhat + t(chol(Vhat + diag(N)*1e-10))%*%rnorm(N), silent=FALSE)
if (is(f.draw, "try-error")) f.draw <- try(t(mvtnorm::rmvnorm(1,fhat, Vhat)), silent = TRUE)
if (is(f.draw, "try-error")) f.draw <- t(mvtnorm::rmvnorm(1,fhat, Matrix::forceSymmetric(Vhat)))

f.store[,irep] <- f.draw + mean.y
  
###------------------------- Step 3.: Draw sigma ----------------------------###
shocks <- y - f.draw 
a.star <- 10 + N/2
b.star <- 0.1 + sum(shocks^2)/2
sigma2 <- 1/rgamma(1, a.star, b.star)
sigma2.draw <- rep(sigma2, N)
if(irep %% 100 == 0) print(irep)
}
  
###--------------------------------------------------------------------------###
###------------------ Obtain posterior/prior quantiles ----------------------###
###--------------------------------------------------------------------------###
qu_prior <- t(apply(prior.store,1,function(x) quantile(x,c(0.05,0.95))))
qu_post <- t(apply(f.store,1,function(x) quantile(x,c(0.01,0.5,0.99))))

###--------------------------------------------------------------------------###
###------------------ Select prior and posterior draws ----------------------###
###--------- Ensure that draws are not outside of credible intervals --------###
###--------------------------------------------------------------------------###
f.draw.slct <- array(0,c(N,no.drws,2))
for(n in 1:no.drws){
  ct = 0 
  while(ct == 0){
      temp <- prior.store[,round(runif(1,1,ntot),0)]
      if(prod(temp <= qu_prior[,2]) & prod(temp >= qu_prior[,1])){
        f.draw.slct[,n,1] <- temp
        ct <- 1
      }
  }
}

###--------------------------------------------------------------------------###
###-------------------------------- Plots -----------------------------------###
###--------------------------------------------------------------------------###

df <- as.data.frame(cbind(y+mean.y,qu_prior,qu_post,X+mean.x))
colnames(df) <- c("y","qu_prior_5","qu_prior_95","qu_post_5","qu_post_50","qu_post_95","X")

df.draw.prior <- as.data.frame(f.draw.slct[,,1])
df.draw.post <- as.data.frame(f.draw.slct[,,2])
colnames(df.draw.prior) <- colnames(df.draw.post) <- paste0("Draw", 1:no.drws)
  
if(run == 1){
  title.temp <- TeX("(a) Linear kernel: \\textbf{$XX'$}")
}else if(run == 2){
 title.temp <- TeX("(b) Gaussian kernel: $\\kappa = 0.01$")
}else if(run == 3){
 title.temp <- TeX("(c) Gaussian kernel: $\\kappa = 0.1$")
}else if(run == 4){
 title.temp <- TeX("(d) Gaussian kernel: $\\kappa = 4$")
}

datevec <- seq(as.Date("2005-01-01"),by="quarter",length.out= T)

# Scatter plot: Prior
ggp_prior[[run]] <- ggplot(df, aes(x=X)) +
    geom_ribbon(aes(ymin = qu_prior_5,ymax = qu_prior_95), alpha = 0.25,fill=color.est) +
    geom_line(data = df.draw.prior,aes(x=df$X,y=Draw1),color = color.est, size = ls.draw, linetype = "dashed") +
    geom_line(data = df.draw.prior,aes(x=df$X,y=Draw2),color = color.est, size = ls.draw, linetype = "dashed") +
    geom_line(data = df.draw.prior,aes(x=df$X,y=Draw3),color = color.est, size = ls.draw, linetype = "dashed") +
    geom_point(aes(y=y), size = ps.obs, color = color.obs) +
    ggtitle(title.temp) +
    ylim(ylim_prior) +
    theme_cowplot(font_size = 14) +
    theme(legend.position="none",plot.margin = unit(c(3,0,0,0),"mm"),plot.title = element_text(family = "sans", size = 18, margin=margin(0,0,1,0)),
        axis.title.x = element_blank(),axis.title.y= element_blank()) + xlab("Macroeconomic uncertainty indicator") #+ {if(run==4)theme(axis.title.x = element_text())} +

# Scatter plot: Posterior
ggp_post[[run]] <- ggplot(df, aes(x=X)) +
    geom_ribbon(aes(ymin = qu_post_5,ymax = qu_post_95), alpha = 0.25,fill= color.est) +
    geom_line(aes(y=qu_post_50),size = ls.est, color = color.est) +
    geom_point(aes(y=y), size = ps.obs, color = color.obs) +
    ggtitle("") +
    ylim(ylim_post) +
    theme_cowplot(font_size = 14) +
    theme(legend.position="none",plot.margin = unit(c(3,0,0,0),"mm"),plot.title = element_text(family = "sans", size = 18, margin=margin(0,0,1,0)),
        axis.title.x = element_blank(),axis.title.y= element_blank()) + xlab("Macroeconomic uncertainty indicator")+ #{if(run==4){theme(axis.title.x = element_text())} } + 
    theme(axis.title.x = element_text()) + theme(plot.margin = unit(c(3,0,3,0),"mm"))

}

ggp <- list()
for(ii in 1:no.kns){
  ggp[[ii]] <- ggp_prior[[ii]]
  ggp[[ii+4]] <- ggp_post[[ii]]
}

# PDF output
pdf(file = paste0(figs.dir, "illustr-",start.date,"-",end.date, "_GDPUNC.pdf"),width = 15,height = 8)
gridExtra::grid.arrange(grobs = lapply(ggp, function(x) x),ncol = 4,widths=c(1.5,1.5,1.5,1.5))
dev.off()
