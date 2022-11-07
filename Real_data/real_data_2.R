source('./R_code/gci_estimate.R')
source('./R_code/knockoff_construct.R')
source('./R_code/knockoff_select.R')
source('./R_code/stem.R')

run <- function(seed) {
  load('./Real_data/Real_data_estimate.RData')
  load("./Real_data/Item_parameters_2015.RData")
  Y <- read.csv('./Real_data/USA_Y.csv')
  
  # Construct knockoff copy
  Gibbs <- X.Gibbs.sample(X, X.under, XX, XX.ind, Y, Omega.es, a.vec, d1.vec, d2.vec,
                          beta.es, theta.interp.es, theta.var.es,
                          bin.ind, ord.ind, con.ind, 
                          bin.params.es, ord.params.es, con.params.es, max_iter=5)
  
  X.under <- Gibbs$X.under.Gibbs
  X.impute <- Impute(X, X.under, bin.ind, ord.ind, con.ind, bin.params.es, ord.params.es, con.params.es)
  res.knocks <- Knockoffs.construct(X.under, Omega.es, S, bin.ind, ord.ind, con.ind,
                                    bin.params.es, ord.params.es, con.params.es)
  X.knock.under <- res.knocks$X.knocks.under
  X.knock <- res.knocks$X.knocks.impute
  X.knock[is.na(X)] <- NA
  
  # MMLE of joint latent regression model
  joint.latreg.es <- StEM(X, X.under, X.knock, X.knock.under, Y, Sigma.es, S, a.vec, d1.vec, d2.vec,
                          bin.ind, ord.ind, con.ind, bin.params.es, ord.params.es, con.params.es, 
                          eps = 1e-2, max_iter = 2000, burn_in = 1000, is.beta.trace=T,
                          is.print.iter=T,  is.knockoff = T)
  save(mmle.res, file = paste0('./Real_data/Knockoff_real_data_final_result_seed_', seed, '.RData'))
}

# Parallelization for 31 knockoff copies
library(parallel)
library(foreach)
library(doParallel)
library(doRNG)
cl <- makeCluster(31)
registerDoParallel(cl)
foreach (seed = 1:31, 
         .combine=cbind, 
         .packages = c('MASS', 'psych','truncnorm','glasso', 'Matrix','arstheta')
         ) %dorng% run(seed)
stopCluster(cl)