# Set working directory as the folder downloaded from 
# https://github.com/ZilongXie/LatentRegMissingKnockoff
setwd('') 

source('./R_code/gci_estimate.R')
source('./R_code/knockoff_construct.R')
source('./R_code/knockoff_select.R')
source('./R_code/stem.R')

bootstrap <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind='default')
  
  X <- read.csv('./Real_data/USA_X.csv')
  Y <- read.csv('./Real_data/USA_Y.csv')
  load('./Real_data/Item_parameters_2015.RData')
  load('./Real_data/Real_data_estimate.RData')
  
  N <- dim(X)[1]
  ind <- sample(1:N, size = N, replace = T)
  X.resample <- X[ind, ]
  Y.resample <- Y[ind, ]
  gci.resampled <- GCI.estimate(X=X.resample,
                                bin.ind=bin.ind,
                                ord.ind=ord.ind,
                                con.ind=con.ind,
                                max.iter = 2000,
                                burn.in = 1000,
                                print.iter = T,
                                is.params.trace = T,
                                is.return.under = T,
                                Sigma0 = Sigma.es)
  
  X.under.resampled <- gci.resampled$X.under
  con.params.resampled <- list()
  con.params.resampled$mean <- gci.resampled$params.avg$con.params.mean
  con.params.resampled$sd <- gci.resampled$params.avg$con.params.sd
  bin.params.resampled <- gci.resampled$params.avg$bin.params
  ord.params.resampled <- gci.resampled$params.avg$ord.params
  Sigma.resampled <- gci.resampled$params.avg$Sigma
  
  latreg.resampled <- StEM(X=X.resample, X.under=X.under.resampled, X.knock=NULL, X.knock.under=NULL, 
                           Y=Y.resample, Sigma=Sigma.resampled, S=NULL, 
                           a.vec=a.vec, d1.vec=d1.vec, d2.vec=d2.vec, 
                           bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind, 
                           bin.params=bin.params.resampled, ord.params=ord.params.resampled, con.params=con.params.resampled, 
                           max_iter=1000, burn_in=500, is.beta.trace=T, is.print.iter=T, is.knockoff=F,
                           theta.interp0 = theta.interp.es, theta.var0 = theta.var.es, theta0 = latreg.es$theta,
                           lambda=0)
  
  beta.resampled <- latreg.resampled$beta
  save(gci.resampled, latreg.resampled,
       beta.resampled, con.params.resampled, 
       file = paste0('./Real_data/bootstrap_', seed, '.RData'))
}

# Parallelization
library(parallel)
library(foreach)
library(doParallel)
library(doRNG)
cl <- makeCluster(100)
registerDoParallel(cl)
foreach (seed = 1:100 + 123456, .combine=cbind,  .packages = c('MASS', 'psych','truncnorm',
                                                               'glasso', 'Matrix','arstheta')) %dorng% bootstrap(seed)
foreach (seed = 101:200 + 123456, .combine=cbind,  .packages = c('MASS', 'psych','truncnorm',
                                                     'glasso', 'Matrix','arstheta')) %dorng% bootstrap(seed)
stopCluster(cl)