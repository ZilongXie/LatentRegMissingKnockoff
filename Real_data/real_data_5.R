load('./Real_data/Real_data_estimate.RData')

# Extended indices
beta.names <- names(beta.es)
cc.ind <- c()
len.temp <- 0
ord.max <- apply(X[,ord.ind], 2, max, na.rm = T)
for (j in 1:dim(X)[2]) {
  if (j %in% con.ind) {
    cc.ind <- c(cc.ind, len.temp + 1)
    len.temp <- len.temp + 1
  }
  if (j %in% bin.ind) {
    len.temp <- len.temp + 1
  }
  if (j %in% ord.ind) {
    loc <- which(ord.ind == j)
    len.temp <- len.temp + ord.max[loc]
  }
}

# Standardized effect size
beta.es[cc.ind] <- beta.es[cc.ind] * con.params.es$sd
names(beta.es) <- beta.names
round(beta.es, 4)

beta.bootstrap <- NULL
for (seed in 1:200 + 123456) {
  load(paste0('./Real_data/bootstrap_',seed,'.RData'))
  beta.resampled[cc.ind] <- beta.resampled[cc.ind] * con.params.resampled$sd
  setNames(beta.resampled, beta.names)
  beta.bootstrap <- rbind(beta.bootstrap, beta.resampled)
}

beta.sd <- sqrt(apply(beta.bootstrap, 2, var))
beta.zscore <- abs(beta.es)/beta.sd
standardized.res <- as.data.frame(cbind(beta.es, beta.sd, beta.zscore))
standardized.res <- apply(standardized.res, 2, round, digit=4)
standardized.res

