# Set working directory as the folder downloaded from 
# https://github.com/ZilongXie/LatentRegMissingKnockoff
setwd('') 

source('./R_code/knockoff_construct.R')
source('./R_code/knockoff_select.R')
load('./Real_data/Real_data_estimate.RData')

###################
# Selected result #
###################
X.scale <- params.to.scale(bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind,
                           bin.params=bin.params.es, ord.params=ord.params.es, con.params=con.params.es)

v <- 1
select.res <- NULL
for (seed in 1:31) {
  load(paste0('./Real_data/Knockoff_real_data_final_result_seed_',seed,'.RData'))
  stats <- knockoff.stat(beta=joint.latreg.es$beta, est.scale=X.scale, 
                         bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind)
  select <- knockoff.select(stats, v)
  select.res <- as.data.frame(rbind(select.res, select))
  names(select.res) <- names(X)
}
select.ind <- which(colMeans(select.res) > 0.5)
selected.names <- names(select.ind)
# Selection probabilities
ranking.prob <- sort(colMeans(select.res), decreasing = T)

v <- 2
select.res <- NULL
for (seed in 1:31) {
  load(paste0('./Real_data/Knockoff_real_data_final_result_seed_',seed,'.RData'))
  stats <- knockoff.stat(beta=joint.latreg.es$beta, est.scale=X.scale, 
                         bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind)
  select <- knockoff.select(stats, v)
  select.res <- as.data.frame(rbind(select.res, select))
  names(select.res) <- names(X)
}
select.ind <- which(colMeans(select.res) > 0.5)
selected.names <- names(select.ind)

v <- 3
select.res <- NULL
for (seed in 1:31) {
  load(paste0('./Real_data/Knockoff_real_data_final_result_seed_',seed,'.RData'))
  stats <- knockoff.stat(beta=joint.latreg.es$beta, est.scale=X.scale, 
                         bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind)
  select <- knockoff.select(stats, v)
  select.res <- as.data.frame(rbind(select.res, select))
  names(select.res) <- names(X)
}

# Variable ranking
ranking.prob[which(names(ranking.prob) %in% selected.names)]








