source('./R_code/knockoff_select.R')
###############
# Single copy #
###############
load('./Simulation/simulation_params.RData')

###########################
## Base selection result ##
###########################
# Specify N here
N <- 1000
# N <- 2000
# N <- 4000

select.1 <- NULL
select.2 <- NULL
select.3 <- NULL
select.4 <- NULL
select.5 <- NULL
for (seed in 1:100) {
  load(paste0('./Simulation/Base_result_', N, '_', seed,'.RData'))
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
# Specify N here
N <- 1000
# N <- 2000
# N <- 4000

select.1 <- NULL
select.2 <- NULL
select.3 <- NULL
select.4 <- NULL
select.5 <- NULL
for (seed in 1:100) {
  load(paste0('./Simulation/Base_result_', N, '_', seed,'.RData'))
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

