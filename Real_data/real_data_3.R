source('./R_code/knockoff_construct.R')
source('./R_code/knockoff_select.R')

load('./Real_data/Real_data_estimate.RData')
X.scale <- params.to.scale(bin.ind, ord.ind, con.ind, bin.params.es, ord.params.es, con.params.es)

###################
# Selected result #
###################
v <- 1
select.res <- NULL
for (seed in 1:31) {
  load(paste0('./Real_data/Knockoff_real_data_final_result_seed_',seed,'.RData'))
  stats <- knockoff.stat(joint.latreg.es$beta, X.scale, bin.ind = bin.ind, ord.ind = ord.ind, con.ind = con.ind)
  select <- knockoff.select(stats, v)
  select.res <- as.data.frame(rbind(select.res, select))
  names(select.res) <- names(X)
}
select.ind <- which(colMeans(select.res) > 0.5)
selected.names <- names(select.ind)

v <- 2
select.res <- NULL
for (seed in 1:31) {
  load(paste0('./Real_data/Knockoff_real_data_final_result_seed_',seed,'.RData'))
  stats <- knockoff.stat(joint.latreg.es$beta, X.scale, bin.ind = bin.ind, ord.ind = ord.ind, con.ind = con.ind)
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
  stats <- knockoff.stat(joint.latreg.es$beta, X.scale, bin.ind = bin.ind, ord.ind = ord.ind, con.ind = con.ind)
  select <- knockoff.select(stats, v)
  select.res <- as.data.frame(rbind(select.res, select))
  names(select.res) <- names(X)
}

select.ind <- which(colMeans(select.res) > 0.5)
selected.names <- names(select.ind)

##################
# Rank variables #
##################
v <- 1
select.res <- NULL
for (seed in 1:31) {
  load(paste0('./Real_data/Knockoff_real_data_final_result_seed_',seed,'.RData'))
  stats <- knockoff.stat(joint.latreg.es$beta, X.scale, bin.ind = bin.ind, ord.ind = ord.ind, con.ind = con.ind)
  select <- knockoff.select(stats, v)
  select.res <- as.data.frame(rbind(select.res, select))
  names(select.res) <- names(X)
}
sort(colMeans(select.res),decreasing = T) 
