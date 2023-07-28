# Set working directory as the folder downloaded from 
# https://github.com/ZilongXie/LatentRegMissingKnockoff
setwd('') 

source('./R_code/gci_estimate.R')
source('./R_code/knockoff_construct.R')
source('./R_code/knockoff_select.R')
source('./R_code/stem.R')


run <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind='default')
  
  load('./Real_data/Real_data_estimate.RData')
  load("./Real_data/Item_parameters_2015.RData")
  Y <- read.csv('./Real_data/USA_Y.csv')

  # Construct knockoff copy
  Gibbs <- X.Gibbs.sample(X=X, X.under=X.under, XX=XX, XX.ind=XX.ind, Y=Y, Omega=Omega.es,
                          a.vec=a.vec, d1.vec=d1.vec, d2.vec=d2.vec,
                          beta=beta.es, theta.interp=theta.interp.es, theta.var=theta.var.es,
                          bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind,
                          bin.params=bin.params.es, ord.params=ord.params.es, con.params=con.params.es,
                          max_iter=5000)
  X.under <- Gibbs$X.under.Gibbs
  save(X.under, .Random.seed,
       file = paste0('./Real_data/Knockoff_real_data_final_Gibbs_seed_', seed, '.RData'))
  
  X.impute <- Impute(X=X, X.under=X.under, 
                     bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind,
                     bin.params=bin.params.es, ord.params=ord.params.es, con.params=con.params.es)
  res.knocks <- Knockoffs.construct(X.under=X.under, Omega=Omega.es, S=S, 
                                    bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind,
                                    bin.params=bin.params.es, ord.params=ord.params.es, con.params=con.params.es)
  X.knock.under <- res.knocks$X.knocks.under
  X.knock <- res.knocks$X.knocks.impute
  X.knock[is.na(X)] <- NA
  
  # MMLE of joint latent regression model
  joint.latreg.es <- StEM(X=X, X.under=X.under, X.knock=X.knock, X.knock.under=X.knock.under, Y=Y, 
                          Sigma=Sigma.es, S=S, a.vec=a.vec, d1.vec=d1.vec, d2.vec=d2.vec, 
                          bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind,
                          bin.params=bin.params.es, ord.params=ord.params.es, con.params=con.params.es, 
                          max_iter=2000, burn_in=1000, is.beta.trace=T, is.print.iter=T, is.knockoff=T,
                          theta0=Gibbs$theta.Gibbs, theta.interp0=theta.interp.es, theta.var0=theta.var.es,
                          lambda=0)  
  
  save(joint.latreg.es, file = paste0('./Real_data/Knockoff_real_data_final_result_seed_', seed, '.RData'))
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