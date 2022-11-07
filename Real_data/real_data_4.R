bootstrap <- function(seed) {
  set.seed(seed)
  X <- read.csv('./Real_data/USA_X.csv')
  Y <- read.csv('./Real_data/USA_Y.csv')
  load('./Real_data/Item_parameters_2015.RData')
  load('./Real_data/Real_data_estimate.RData')
  
  N <- dim(X)[1]
  ind <- sample(1:N, size = N, replace = T)
  X.resample <- X[ind, ]
  Y.resample <- Y[ind, ]
  gci.resampled <- GCI.estimate(X.resample,
                                bin.ind,
                                ord.ind,
                                con.ind,
                                max.iter = 1000,
                                burn.in = 500,
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
  
  latreg.resampled <- StEM(X.resample, X.under.resampled, NULL, NULL, Y.resample, Sigma.resampled, NULL,
                           a.vec, d1.vec, d2.vec, 
                           bin.ind, ord.ind, con.ind, 
                           bin.params.resampled, ord.params.resampled, con.params.resampled, 
                           max_iter=1000, burn_in=500, is.beta.trace=T, is.print.iter=T,
                           beta0 = beta.es, theta.interp0 = theta.interp.es,
                           theta.var0 = theta.var.es, theta0 = latreg.es$theta)
  
  beta.resampled <- latreg.resampled$beta
  save(beta.resampled, con.params.resampled, file = paste0('./Real_data/bootstrap_', seed, '.RData'))
}

source('./R_code/gci_estimate.R')
source('./R_code/stem.R')
source('./R_code/knockoff_construct.R')
source('./R_code/knockoff_select.R')
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