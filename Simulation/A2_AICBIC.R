# Set working directory as the folder downloaded from 
# https://github.com/ZilongXie/LatentRegMissingKnockoff
setwd('') 

################################################################################
source('Knockoff_Simulation/gci_estimate.R')
source('Knockoff_Simulation/knockoff_construct.R')
source('Knockoff_Simulation/knockoff_select.R')
source('Knockoff_Simulation/stem.R')

expit <- function(x) {1/(1+exp(-x))}

# Data generation
data.generate <- function(N, p, J, Sigma.true, bin.ind, con.ind, bin.params.true,
                          con.params.true, beta.true, a.vec, d.vec, booklets) {
  # Mixed data
  X.under <- mvtnorm::rmvnorm(n = N, mean = rep(0, p), Sigma.true)
  X <- X.under
  X[,bin.ind] <- t(t(X.under[,bin.ind]) > bin.params.true) * 1
  X[,con.ind] <- t(t(X.under[,con.ind]) * con.params.true$sd + con.params.true$mean)
  X.complete <- X
  
  # Latent trait
  theta <- as.vector(X %*% beta.true + rnorm(N))
  theta.complete <- theta
  
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
  
  return(list(X.complete=X.complete, X.under=X.under, theta.complete=theta.complete, X=X, Y=Y))
}

# Automatic stepwise variable selection based on AIC/BIC
auto.stepwise.selection <- function(X, theta, criterion = 'AIC') {
  N <- dim(X)[1]
  p <- dim(X)[2]
  
  df <- data.frame(cbind(X, theta))
  names(df) <- c(paste0('X', 1:p), 'theta')
  full <- lm('theta~.', data = df)
  
  temp = full
  var.len = length(temp$coefficients)
  var.set = sort(sapply(strsplit(names(temp$coefficients)[-1], split = 'X'), 
                        function(x){as.numeric(x[2])}))
  
  for (i in 1:10000) {
    var.len.old <- var.len
    var.set.old <- var.set
    
    if (criterion == 'AIC') {
      multiply <- 2
    }
    
    if (criterion == 'BIC') {
      multiply <- log(N)
    }
    
    temp = MASS::stepAIC(temp,
                         direction = 'both',
                         scope = list(lower = 'theta~1', upper = 'theta~.'),
                         steps = 1,
                         k=multiply,
                         trace = F)
    
    var.len = length(temp$coefficients)
    var.set = sort(sapply(strsplit(names(temp$coefficients)[-1], split = 'X'), 
                          function(x){as.numeric(x[2])}))
    
    if (var.len == var.len.old) {
      if (identical(var.set, var.set.old)) {
        break()
      }
    }
  }
  
  return(list(var.set = var.set, num.steps = i))
}

run <- function(seed, N) {
  set.seed(seed)
  load('./Simulation/simulation_params.RData')
  
  data <- data.generate(N, p, J, Sigma.true, bin.ind, con.ind, bin.params.true,
                        con.params.true, beta.true, a.vec, d1.vec, booklets)
  X <- data$X
  Y <- data$Y
  X.complete <- data$X.complete
  theta.complete <- data$theta.complete
  X.under <- data$X.under
  
  # Stepwise AIC/BIC selection on complete data
  AIC.complete <- auto.stepwise.selection(X = X.complete, theta = theta.complete, criterion = 'AIC')
  BIC.complete <- auto.stepwise.selection(X = X.complete, theta = theta.complete, criterion = 'BIC')
  
  # Gibbs sampling from true model
  Omega.true <- solve(Sigma.true)
  theta.interp.true <- 0
  theta.var.true <- 1
  XX <- X.complete  
  XX.ind <- list()
  for (j in 1:p) {XX.ind[[j]] <- j}
  
  Gibbs <- X.Gibbs.sample(X, X.under, XX, XX.ind, Y, Omega.true, 
                          a.vec, d1.vec, d2.vec,
                          beta.true, theta.interp.true, theta.var.true,
                          bin.ind, ord.ind, con.ind, 
                          bin.params.true, ord.params.true, con.params.true, 
                          max_iter = 1000)
  
  X.Gibbs <- Gibbs$XX.Gibbs
  theta.Gibbs <- Gibbs$theta.Gibbs
  
  # Stepwise AIC/BIC selection on imputed data
  AIC.impute.true <- auto.stepwise.selection(X = X.Gibbs, theta = theta.Gibbs, criterion = 'AIC')
  BIC.impute.true <- auto.stepwise.selection(X = X.Gibbs, theta = theta.Gibbs, criterion = 'BIC')
  
  
  # Gibbs sampling from estimated model
  load(paste0('Knockoff_Simulation/final_derand_V2_', N, '_', seed, '.RData'))
  
  con.params.es <- list()
  con.params.es$mean <- estimate.res$params.avg$con.params.mean
  con.params.es$sd <- estimate.res$params.avg$con.params.sd
  bin.params.es <- estimate.res$params.avg$bin.params
  ord.params.es <- estimate.res$params.avg$ord.params
  Omega.es <- solve(estimate.res$params.avg$Sigma)    
  beta.es <- latreg.es$beta
  theta.interp.es <- latreg.es$interp
  theta.var.es <- latreg.es$var
  X.under.es <- latreg.es$X.under
  XX.es <- latreg.es$XX
  XX.ind <- latreg.es$XX.ind
  
  Gibbs <- X.Gibbs.sample(X, X.under.es, XX.es, XX.ind, Y, Omega.es,
                          a.vec, d1.vec, d2.vec,
                          beta.es, theta.interp.es, theta.var.es,
                          bin.ind, ord.ind, con.ind, 
                          bin.params.es, ord.params.es, con.params.es, 
                          max_iter = 1000)
  
  X.Gibbs <- Gibbs$XX.Gibbs
  theta.Gibbs <- Gibbs$theta.Gibbs
  
  # Stepwise AIC/BIC selection on imputed data
  AIC.impute.es <- auto.stepwise.selection(X = X.Gibbs, theta = theta.Gibbs, criterion = 'AIC')
  BIC.impute.es <- auto.stepwise.selection(X = X.Gibbs, theta = theta.Gibbs, criterion = 'BIC')
  
  
  save(AIC.complete, BIC.complete, AIC.impute.true, BIC.impute.true, AIC.impute.es, BIC.impute.es,
       file = paste0('./Simulation/Base_AICBIC_result_', N, '_', seed, '.RData'))
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