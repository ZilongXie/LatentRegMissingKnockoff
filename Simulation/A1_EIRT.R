# Set working directory as the folder downloaded from 
# https://github.com/ZilongXie/LatentRegMissingKnockoff
setwd('') 

################################################################################
source('./R_code/gci_estimate.R')
source('./R_code/knockoff_construct.R')
source('./R_code/knockoff_select.R')
source('./R_code/stem.R')

# Utility function
expit <- function(x) {1/(1+exp(-x))}

# Data generation
data.generate <- function(N, p, J, Sigma.true, bin.ind, con.ind, bin.params.true,
                          con.params.true, beta.true, a.vec, d.vec, booklets) {
  # Mixed data
  X.under <- mvtnorm::rmvnorm(n = N, mean = rep(0, p), Sigma.true)
  X <- X.under
  X[,bin.ind] <- t(t(X.under[,bin.ind]) > bin.params.true) * 1
  X[,con.ind] <- t(t(X.under[,con.ind]) * con.params.true$sd + con.params.true$mean)
  
  # Latent trait
  theta <- as.vector(X %*% beta.true + rnorm(N))
  
  # Item response
  Y.prob <- expit(outer(theta, a.vec) + outer(rep(1, N), d.vec))
  Y <- matrix(rbinom(n = prod(dim(Y.prob)), size = 1, prob = Y.prob),
              nrow = dim(Y.prob)[1], ncol = dim(Y.prob)[2])
  
  miss.prob <- matrix(0, nrow = N, ncol = p)
  true.ind <- which(beta.true != 0)
  
  # Missing at random
  groups <- sample(1:5, N, replace = T)
  for (i in 1:5) {
    rows <- which(groups == i)
    miss.prob[rows,] <- outer(expit(rowMeans(X[rows,true.ind[(i-1)*2 + 1:2]]) - 1), rep(1, p))
    miss.prob[rows, true.ind[(i-1)*2 + 1:2]] <- 0
  }
  
  miss.mask <- (matrix(rbinom(n = prod(dim(miss.prob)), size = 1, prob = miss.prob), 
                       nrow = dim(miss.prob)[1], ncol = dim(miss.prob)[2]) == 1)
  X[miss.mask] <- NA
  
  # Booklets design
  Y.miss.mask <- booklets[sample(1:3, N, replace = T),] == 1
  Y[Y.miss.mask] <- NA
  
  return(list(X=X, Y=Y))
}

run <- function(seed, N) {
  set.seed(seed)
  load('./Simulation/simulation_params.RData')
  
  data <- data.generate(N, p, J, Sigma.true, bin.ind, con.ind, bin.params.true,
                        con.params.true, beta.true, a.vec, d1.vec, booklets)
  X <- data$X
  Y <- data$Y

  # Estimation of Gaussian copula model
  estimate.res <- GCI.estimate(X,   
                               bin.ind = bin.ind,
                               ord.ind = ord.ind,
                               con.ind = con.ind,
                               max.iter = 1000,
                               burn.in = 500,
                               print.iter = T,
                               is.params.trace = T,
                               is.return.under = T,
                               Sigma0 = Sigma.true,
                               bin.params0 = bin.params.true,
                               con.params0 = con.params.true
  )
  
  con.params.es <- list()
  con.params.es$mean <- estimate.res$params.avg$con.params.mean
  con.params.es$sd <- estimate.res$params.avg$con.params.sd
  bin.params.es <- estimate.res$params.avg$bin.params
  ord.params.es <- estimate.res$params.avg$ord.params
  Sigma.es <- estimate.res$params.avg$Sigma
  Omega.es <- solve(Sigma.es)
  X.under0 <- estimate.res$X.under
  S <- MVR.knockoffs(Sigma.es)
  
  latreg.es <- StEM(X, X.under0, NULL, NULL, Y, Sigma.es, NULL, NULL, NULL, NULL, #==Estimate IRT model
                    bin.ind, ord.ind, con.ind, bin.params.es, ord.params.es, con.params.es, 
                    max_iter = 1000, burn_in = 500, is.beta.trace=T,
                    is.print.iter=T, is.knockoff = F)
  
  beta.es <- latreg.es$beta
  X.under <- latreg.es$X.under
  XX <- latreg.es$XX
  XX.ind <- latreg.es$XX.ind
  
  theta.interp.es <- latreg.es$interp
  theta.var.es <- latreg.es$var
  a.vec.es <- latreg.es$a.vec
  d1.vec.es <- latreg.es$d1.vec
  d2.vec.es <- latreg.es$d2.vec
  
  knockoff.stat.record <- NULL
  num.copy <- 31
  # Derandomized knockoff
  for (copy in 1:num.copy) {
    Gibbs <- X.Gibbs.sample(X, X.under, XX, XX.ind, Y, Omega.es, a.vec.es, d1.vec.es, d2.vec.es,
                            beta.es, theta.interp.es, theta.var.es,
                            bin.ind, ord.ind, con.ind,
                            bin.params.es, ord.params.es, con.params.es, max_iter = 1000)
    
    # Construct knockoff copy
    X.under <- Gibbs$X.under.Gibbs
    X.impute <- Impute(X, X.under, bin.ind, ord.ind, con.ind, bin.params.es, ord.params.es, con.params.es)
    res.knocks <- Knockoffs.construct(X.under, Omega.es, S, bin.ind, ord.ind, con.ind,
                                      bin.params.es, ord.params.es, con.params.es)
    X.knock.under <- res.knocks$X.knocks.under
    X.knock <- res.knocks$X.knocks.impute
    X.knock[is.na(X)] <- NA
    
    # MMLE
    # MMLE of joint latent regression model
    joint.latreg.es <- StEM(X, X.under, X.knock, X.knock.under, Y, Sigma.es, S, a.vec.es, d1.vec.es, d2.vec.es,
                            bin.ind, ord.ind, con.ind, bin.params.es, ord.params.es, con.params.es, 
                            eps = 1e-2, max_iter = 1000, burn_in = 500, is.beta.trace=T,
                            is.print.iter=T,  is.knockoff = T)
    beta.avg <- joint.latreg.es$beta
    
    # Construct knockoff statistics
    X.scale <- params.to.scale(bin.ind, ord.ind, con.ind, bin.params.es, ord.params.es, con.params.es)
    stats <- knockoff.stat(beta.avg, X.scale)
    knockoff.stat.record <- rbind(knockoff.stat.record, stats)
  }  
  
  select.res <- matrix(0, nrow = 5, ncol = p)  
  for (v in 1:5) {
    select.v <- NULL
    for (i in 1:num.copy) {
      stats <- knockoff.stat.record[i,]
      select.v <- rbind(select.v, knockoff.select(stats, v = v))
    }
    select.res[v,] <- (colMeans(select.v) > 0.5) * 1
  }
  save(estimate.res, latreg.es, knockoff.stat.record, select.res, file = paste0('./Simulation/Base_EIRT_result_', N, '_', seed, '.RData'))
}

library(parallel)
library(foreach)
library(doParallel)
library(doRNG)

cl <- makeCluster(100)
registerDoParallel(cl)
foreach (seed = 1:100, .combine=cbind,  .packages = c('MASS', 'psych','truncnorm', 'Matrix','arstheta')) %dorng% run(seed = seed, N = 1000)
stopCluster(cl)

cl <- makeCluster(100)
registerDoParallel(cl)
foreach (seed = 1:100, .combine=cbind,  .packages = c('MASS', 'psych','truncnorm', 'Matrix','arstheta')) %dorng% run(seed = seed, N = 2000)
stopCluster(cl)

cl <- makeCluster(100)
registerDoParallel(cl)
foreach (seed = 1:100, .combine=cbind,  .packages = c('MASS', 'psych','truncnorm', 'Matrix','arstheta')) %dorng% run(seed = seed, N = 4000)
stopCluster(cl)

