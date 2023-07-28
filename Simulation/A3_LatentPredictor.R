# Set working directory as the folder downloaded from 
# https://github.com/ZilongXie/LatentRegMissingKnockoff
setwd('')


#############
# Functions #
#############
source('./R_code/gci_estimate.R')
source('./R_code/knockoff_construct.R')
source('./R_code/knockoff_select.R')
source('./R_code/stem.R')

# Utility function
expit <- function(x) {1/(1+exp(-x))}

###########################
## Parameters generation ##  
###########################
J <- 60
p <- 100

J <- 60
p <- 100
set.seed(123456789, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind='default')

# Block-wise underlying correlation matrix
block1 <- 1:20
block2 <- 21:40 
block3 <- 41:60
block4 <- 61:80
block5 <- 81:100
Sigma.temp <- matrix(0, nrow = p, ncol = p)

Sigma.temp[block1, block1] <- 0.6
Sigma.temp[block2, block2] <- 0.6
Sigma.temp[block3, block3] <- 0.6
Sigma.temp[block4, block4] <- 0.3
Sigma.temp[block5, block5] <- 0.3

Sigma.temp[block1, block2] <- 0.15
Sigma.temp[cbind(block1, block2)] <- 0.3
Sigma.temp[block1, block3] <- 0.15
Sigma.temp[block2, block3] <- 0.15

Sigma.temp[block1, block4] <- 0.15
Sigma.temp[cbind(block1, block4)] <- 0.3
Sigma.temp[block2, block4] <- 0.15
Sigma.temp[cbind(block2, block4)] <- 0.15
Sigma.temp[block3, block4] <- 0.15

# Sigma.temp[block1, block1] <- 0.4
# Sigma.temp[block2, block2] <- 0.4
# Sigma.temp[block3, block3] <- 0.4
# Sigma.temp[block4, block4] <- 0.2
# Sigma.temp[block5, block5] <- 0.2

# Sigma.temp[block1, block2] <- 0.1
# Sigma.temp[cbind(block1, block2)] <- 0.2
# Sigma.temp[block1, block3] <- 0.1
# Sigma.temp[block2, block3] <- 0.1

# Sigma.temp[block1, block4] <- 0.1
# Sigma.temp[cbind(block1, block4)] <- 0.2
# Sigma.temp[block2, block4] <- 0.1
# Sigma.temp[cbind(block2, block4)] <- 0.1
# Sigma.temp[block3, block4] <- 0.1

Sigma.temp[c(block1, block2, block3, block4), block5] <- matrix(runif(80 * 20, min = 0.1, max = 0.2), nrow = 80)

Sigma.temp[lower.tri(Sigma.temp)] <- t(Sigma.temp)[lower.tri(Sigma.temp)]
diag(Sigma.temp) <- 1

########################
## Unknown parameters ##
########################
bin.ind.temp <- c(1:10, 21:30, 41:50, 61:70, 81:90)
con.ind.temp <- c(11:20, 31:40, 51:60, 71:80, 91:100)
mea.ind.temp <- c(11, 13, 15, 17, 19,
                  31, 33, 35, 37, 39,
                  51, 53, 55, 57, 59,
                  71, 73, 75, 77, 79,
                  91, 93, 95, 97, 99)
con.obs.ind.temp <- con.ind.temp[which(!(con.ind.temp %in% mea.ind.temp))]
beta.temp <- rep(0, 100)
beta.temp[c(1,11,22,32,43,53,64,74,85,95)] <- c(1, -1, 1, -1, 1, -1, 1, -1, 1, -1)/2

#######################
## Rearrange indices ##
#######################
bin.ind <- 1:length(bin.ind.temp)
con.ind <- length(bin.ind.temp) + 1:length(con.obs.ind.temp)
mea.ind <- length(bin.ind.temp) + length(con.obs.ind.temp) + 1:length(mea.ind.temp)
ord.ind <- NULL
beta.true <- beta.temp[c(bin.ind.temp, con.obs.ind.temp, mea.ind.temp)]
Sigma.true <- Sigma.temp[c(bin.ind.temp, con.obs.ind.temp, mea.ind.temp),][,c(bin.ind.temp, con.obs.ind.temp, mea.ind.temp)]

con.params.true <- list()
con.params.true$mean <- rep(0, 25)
con.params.true$sd <- rep(1, 25)
bin.params.true <- as.vector(matrix(c(-1.2, -0.3, 0, 0.3, 1.2), nrow = 5, ncol = 10))
ord.params.true <- NULL
mea.params.true <- list()
mea.params.true$mean <-list()
mea.params.true$loadings <- list()
mea.params.true$errvar <- list()
for (l in 1:length(mea.ind.temp)) {
  mea.params.true$mean[[l]] <- rep(0, 5)
  mea.params.true$loadings[[l]] <- rep(1, 5)
  mea.params.true$errvar[[l]] <- rep(1, 5)
}

####################
## IRT parameters ##
####################
a.vec <- runif(J, 0.5, 1)
d1.vec <- runif(J, -2, 0)
d2.vec <- NA * d1.vec

##############
## Booklets ##
##############
booklets <- matrix(0, nrow = 3, ncol = J)
booklets[1, 1:20] <- 1
booklets[2, 21:40] <- 1
booklets[3, 41:60] <- 1


#####################
## Data generation ##
#####################
data.generate.latent <- function(N, p, J, Sigma.true, bin.ind, con.ind, mea.ind,
                                 bin.params.true, con.params.true, mea.params.true,
                                 beta.true, a.vec, d.vec, booklets) {
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
  
  # Factor models for predictors
  W.list <- list()
  for (l in 1:length(mea.ind)) {
    W <- outer(X.under[,mea.ind[l]], mea.params.true$loadings[[l]])
    W <- W + outer(rep(1, N), mea.params.true$mean[[l]])
    W <- W + mvtnorm::rmvnorm(n = N, mean = rep(0, ncol(W)), 
                              sigma = diag(mea.params.true$errvar[[l]]))
    W.list[[l]] <- W
  }
  
  # SMAR condition
  miss.prob <- matrix(0, nrow = N, ncol = p)
  true.ind <- which(beta.true != 0)
  true.ind <- true.ind[c(1, 8,
                         2, 6,
                         3, 9,
                         4, 7,
                         5, 10)]
  
  X.temp <- X 
  for (l in 1:length(mea.ind)) {
    X.temp[,mea.ind[l]] <- rowMeans(W.list[[l]])/2 # Substitution with mean of W  
  }
  
  groups <- sample(1:5, N, replace = T)
  for (i in 1:5) {
    rows <- which(groups == i)
    
    miss.prob[rows,] <- outer(expit(rowMeans(X.temp[rows,true.ind[(i-1)*2 + 1:2]]) - 1), rep(1, p))
    miss.prob[rows, true.ind[(i-1)*2 + 1:2]] <- 0
    miss.prob[rows, mea.ind] <- 0
  }
  miss.mask <- (matrix(rbinom(n = prod(dim(miss.prob)), size = 1, prob = miss.prob), 
                       nrow = dim(miss.prob)[1], ncol = dim(miss.prob)[2]) == 1)
  X.temp[miss.mask] <- NA
  X <- X.temp[,c(bin.ind, con.ind)]
  
  # Booklets design
  Y.miss.mask <- booklets[sample(1:3, N, replace = T),] == 1
  Y[Y.miss.mask] <- NA
  
  return(list(X=X, 
              W.list=W.list, 
              Y=Y))
}

run <- function(seed, N) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind='default')
  
  data <- data.generate.latent(N=N, p=p, J=J, Sigma.true=Sigma.true, bin.ind=bin.ind, 
                               con.ind=con.ind, mea.ind=mea.ind, 
                               bin.params.true=bin.params.true,
                               con.params.true=con.params.true,
                               mea.params.true=mea.params.true,
                               beta.true=beta.true, 
                               a.vec=a.vec, 
                               d.vec=d1.vec, 
                               booklets=booklets)
  X <- data$X
  Y <- data$Y
  W.list <- data$W.list
  
  
  # Estimation of Gaussian copula model
  estimate.res <- GCI.estimate(X = X,
                               bin.ind = bin.ind,
                               ord.ind = ord.ind,
                               con.ind = con.ind,
                               max.iter = 2000,
                               burn.in = 1000,
                               print.iter = T,
                               is.params.trace = T,
                               print.max.change = F,
                               is.return.under = T,
                               W.list = W.list,
                               Sigma0 = Sigma.true,
                               bin.params0 = bin.params.true,
                               con.params0 = con.params.true,
                               mea.params = mea.params.true)

  con.params.es <- list()
  con.params.es$mean <- estimate.res$params.avg$con.params.mean
  con.params.es$sd <- estimate.res$params.avg$con.params.sd
  bin.params.es <- estimate.res$params.avg$bin.params
  ord.params.es <- estimate.res$params.avg$ord.params
  Sigma.es <- estimate.res$params.avg$Sigma
  Omega.es <- solve(Sigma.es)
  X.under0 <- estimate.res$X.under
  S <- MVR.knockoffs(Sigma.es)
  
  # Estimation of latent regression IRT model
  latreg.es <- StEM(X=X, X.under=X.under0, X.knock=NULL, X.knock.under=NULL, 
                    Y=Y, Sigma=Sigma.es, S=NULL,  
                    a.vec=a.vec, d1.vec=d1.vec, d2.vec=d2.vec,
                    bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind, 
                    bin.params=bin.params.es, ord.params=ord.params.es, 
                    con.params=con.params.es, max_iter=1000, burn_in=500, 
                    is.beta.trace=F, is.print.iter=T, is.knockoff=F,
                    W.list=W.list, mea.params=mea.params.true, 
                    lambda = sqrt(1/N))
  beta.es <- latreg.es$beta
  X.under <- latreg.es$X.under
  theta.interp.es <- latreg.es$interp
  theta.var.es <- latreg.es$var
  XX <- latreg.es$XX
  XX.ind <- latreg.es$XX.ind
  
  
  knockoff.stat.record <- NULL
  num.copy <- 31
  # Derandomized knockoff
  for (copy in 1:num.copy) {
    Gibbs <- X.Gibbs.sample(X=X, X.under=X.under, XX=XX, XX.ind=XX.ind, Y=Y, Omega=Omega.es,
                            a.vec=a.vec, d1.vec=d1.vec, d2.vec=d2.vec,
                            beta=beta.es, theta.interp=theta.interp.es, theta.var=theta.var.es,
                            bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind, 
                            bin.params=bin.params.es, ord.params=ord.params.es, 
                            con.params=con.params.es, W.list=W.list, 
                            mea.ind=mea.ind, mea.params=mea.params.true,
                            max_iter=1000)
    
    # Construct knockoff copy
    X.under <- Gibbs$X.under.Gibbs
    X.impute <- Impute(X=X, X.under=X.under, bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind, 
                       bin.params=bin.params.es, ord.params=ord.params.es, con.params=con.params.es)
    res.knocks <- Knockoffs.construct(X.under=X.under, Omega=Omega.es, S=S, 
                                      bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind,
                                      bin.params=bin.params.es, ord.params=ord.params.es, 
                                      con.params=con.params.es, mea.ind = mea.ind,
                                      mea.params = mea.params.true)
    X.knock.under <- res.knocks$X.knocks.under
    X.knock <- res.knocks$X.knocks.impute
    W.list.knock <- res.knocks$W.list.knock
    X.knock[is.na(X)] <- NA
    
    # MMLE of joint latent regression model
    joint.latreg.es <- StEM(X=X, X.under=X.under, X.knock=X.knock, X.knock.under=X.knock.under, 
                            Y=Y, Sigma=Sigma.es, S=S, a.vec=a.vec, d1.vec=d1.vec, d2.vec=d2.vec,
                            bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind, 
                            bin.params=bin.params.es, ord.params=ord.params.es, con.params=con.params.es, 
                            max_iter=1000, burn_in=500, is.beta.trace=T, is.print.iter=T, is.knockoff=T,
                            W.list=W.list, mea.params=mea.params.true, W.list.knock=W.list.knock,
                            lambda = sqrt(1/N))
    
    beta.avg <- joint.latreg.es$beta
    
    # Construct knockoff statistics
    X.scale <- params.to.scale(bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind,
                               bin.params=bin.params.es, ord.params=ord.params.es,
                               con.params=con.params.es, mea.ind=mea.ind)
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
  save(estimate.res, latreg.es, knockoff.stat.record, select.res, file = paste0('./Simulation/Base_Latent_result_', N, '_', seed, '.RData'))
}

#######################
# Parallel simulation #
#######################
library(parallel)
library(foreach)
library(doParallel)
library(doRNG)

cl <- makeCluster(100)
registerDoParallel(cl)
foreach (seed = 1:100, .combine=cbind,  .packages = c('MASS', 'psych','truncnorm', 'Matrix','arstheta','glmnet')) %dorng% run(seed = seed, N = 1000)
stopCluster(cl)

cl <- makeCluster(100)
registerDoParallel(cl)
foreach (seed = 1:100, .combine=cbind,  .packages = c('MASS', 'psych','truncnorm', 'Matrix','arstheta','glmnet')) %dorng% run(seed = seed, N = 2000)
stopCluster(cl)

cl <- makeCluster(100)
registerDoParallel(cl)
foreach (seed = 1:100, .combine=cbind,  .packages = c('MASS', 'psych','truncnorm', 'Matrix','arstheta','glmnet')) %dorng% run(seed = seed, N = 4000)
stopCluster(cl)

###########
# Results #
###########
N = 1000 # Modify N here, N = 1000, 2000, 4000
###########################
## Base selection result ##
###########################
select.1 <- NULL
select.2 <- NULL
select.3 <- NULL
select.4 <- NULL
select.5 <- NULL
for (seed in 1:100) {
  load(paste0('./Simulation/Base_Latent_result_', N, '_', seed,'.RData'))
  
  select.1 <- rbind(select.1, knockoff.select(knockoff.stat.record[1,], 1))
  select.2 <- rbind(select.2, knockoff.select(knockoff.stat.record[1,], 2))
  select.3 <- rbind(select.3, knockoff.select(knockoff.stat.record[1,], 3))
  select.4 <- rbind(select.4, knockoff.select(knockoff.stat.record[1,], 4))
  select.5 <- rbind(select.5, knockoff.select(knockoff.stat.record[1,], 5))
}   

# PFER
mean(rowSums(t(t(select.1 == 1) * (beta.true == 0))))
mean(rowSums(t(t(select.2 == 1) * (beta.true == 0))))
mean(rowSums(t(t(select.3 == 1) * (beta.true == 0))))
mean(rowSums(t(t(select.4 == 1) * (beta.true == 0))))
mean(rowSums(t(t(select.5 == 1) * (beta.true == 0))))

# TPR
mean(rowSums(t(t(select.1 == 1) * (beta.true != 0))))
mean(rowSums(t(t(select.2 == 1) * (beta.true != 0))))
mean(rowSums(t(t(select.3 == 1) * (beta.true != 0))))
mean(rowSums(t(t(select.4 == 1) * (beta.true != 0))))
mean(rowSums(t(t(select.5 == 1) * (beta.true != 0))))

###################################
## Derandomized selection result ##
###################################
select.1 <- NULL
select.2 <- NULL
select.3 <- NULL
select.4 <- NULL
select.5 <- NULL
for (seed in 1:100) {
  load(paste0('./Simulation/Base_Latent_result_', N, '_', seed,'.RData'))
  
  select.1 <- rbind(select.1, select.res[1,])
  select.2 <- rbind(select.2, select.res[2,])
  select.3 <- rbind(select.3, select.res[3,])
  select.4 <- rbind(select.4, select.res[4,])
  select.5 <- rbind(select.5, select.res[5,])
}

# PFER
mean(rowSums(t(t(select.1 == 1) * (beta.true == 0))))
mean(rowSums(t(t(select.2 == 1) * (beta.true == 0))))
mean(rowSums(t(t(select.3 == 1) * (beta.true == 0))))
mean(rowSums(t(t(select.4 == 1) * (beta.true == 0))))
mean(rowSums(t(t(select.5 == 1) * (beta.true == 0))))

# Power
mean(rowSums(t(t(select.1 == 1) * (beta.true != 0))))
mean(rowSums(t(t(select.2 == 1) * (beta.true != 0))))
mean(rowSums(t(t(select.3 == 1) * (beta.true != 0))))
mean(rowSums(t(t(select.4 == 1) * (beta.true != 0))))
mean(rowSums(t(t(select.5 == 1) * (beta.true != 0))))
