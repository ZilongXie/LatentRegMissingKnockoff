# Truncation function
trunc <- function(x, a, b) {
  pmax(pmin(x, b), a)
}

# Compute row min
rowmin <- function(x) {
  apply(x, 1, min, na.rm = T)
}

# Convert ordinal variables to dummy variables
ord.to.dum <- function(X.ord) {
  # Convert to matrix if the input is a vector
  if (is.vector(X.ord)) {
    X.ord <- as.matrix(X.ord)
  }
  
  N <- dim(X.ord)[1]
  p <- dim(X.ord)[2]
  
  if (!all(apply(X.ord, 2, min, na.rm = T) == 0)) {
    stop('Ordinal variables should start from 0!', call. = F)
  }
  
  # Maximum of each variable
  ord.max <- apply(X.ord, 2, max, na.rm = T)
  
  # Convert ordinal variables to dummy variables
  X.ord.dum <- NULL
  for (j in 1:p) {
    dum <- matrix(0, nrow = N, ncol = ord.max[j])
    mask <- (outer(rep(1, N), (1:ord.max[j])) <= outer(X.ord[,j], rep(1, ord.max[j])))
    dum[mask] <- 1
    X.ord.dum <- cbind(X.ord.dum, dum)
  }
  return(X.ord.dum)
}

# A fast weighted sampling function with unstandardized probabilies
weighted.sampling <- function(prob.mat, value.mat) {
  N <- dim(prob.mat)[1]
  p <- dim(prob.mat)[2]
  
  # Standardize probabiliy matrix
  log.prob.mat <- log(prob.mat)
  prob.mat <- exp(log.prob.mat - rowmin(log.prob.mat))
  prob.mat <- prob.mat/rowSums(prob.mat)
  
  # Calculate cumsum of each row without apply
  cumprob.mat <- prob.mat %*% upper.tri(diag(p), diag = TRUE)
  
  # Random weighted sampling
  U <- runif(N)
  index <- rowSums(U > cumprob.mat) + 1L
  
  return(value.mat[cbind(1:N, index)])
}

# Sampling the latent trait theta by adaptive rejection sampling
theta.sample <- function(XX, Y, beta, theta.interp, theta.var, a.vec, d1.vec, d2.vec) {
  Y <- as.matrix(Y)
  XX <- as.matrix(XX)
  mu <- XX %*% beta + theta.interp
  theta.sampler <- arstheta::ThetaPosterior$new(Y, a.vec, d1.vec, d2.vec, mu, theta.var)
  theta.ars <- theta.sampler$sample(mu - 10, mu + 10)
  
  return(theta.ars)
}


# Sampling the underlying variables (for both true and knockoffs) and impute the 
# standardized covariate matrix 
XX.sample <- function(XX, X, X.under, X.knock.under=NULL, theta, beta, 
                      Omega, XX.ind, bin.params=NULL, ord.params=NULL, 
                      con.params=NULL, bin.ind=NULL, ord.ind=NULL, con.ind=NULL, 
                      thresh.bin=NULL, thresh.ord=NULL, is.knockoff=F, ...) {
  N <- dim(X.under)[1]
  p <- dim(X.under)[2]
  XX <- as.matrix(XX)
  args <- list(...)
  
  pp <- dim(Omega)[1]
  # Allowed to provide quantities to save computation
  if (is.null(args$meantrans.mat)) {
    matrix(0, nrow = pp-1 , ncol = pp)
    for (i in 1:pp) {
      meantrans.mat[, i] <- -Omega[i, -i]/Omega[i, i]
    }
  } else {
    meantrans.mat <- args$meantrans.mat
  }
  
  
  if (is.null(args$miss.mask)) {
    miss.mask <- is.na(X)
  } else {
    miss.mask <- args$miss.mask
  }
  
  if (!is.null(ord.ind)) {
    if (is.null(args$ord.max)) {
      ord.max <- ord.max <- apply(X[,ord.ind], 2, max, na.rm = T)
    } else {
      ord.max <- args$ord.max
    }
  }
  
  # One-round Gibbs sampling with random scanning
  random <- sample(1:pp)
  for (j in random) {
    cond.mean <- cbind(X.under, X.knock.under)[, -j] %*% meantrans.mat[, j]
    cond.sd <- sqrt(1/Omega[j, j])
    jj <- XX.ind[[j]]
    
    # True covariates
    if (j <= p) {
      miss.mask.j <- miss.mask[,j]
      miss.num.j <- sum(miss.mask.j)
      
      # Continuous variables
      if(j %in% con.ind) {
        if(miss.num.j > 0) {
          loc <- which(con.ind == j)
          mean.miss <- (cond.sd^2 * beta[jj] * (theta - XX[,-jj] %*% beta[-jj]) + cond.mean)[miss.mask.j]/(1 + beta[jj]^2 * cond.sd^2)
          sd.miss <- cond.sd/sqrt(1 + cond.sd^2 * beta[jj]^2)
          
          # Update X.under
          X.under[miss.mask.j,j] <- mean.miss + rnorm(n = miss.num.j, mean = 0, sd = sd.miss)
          
          # Update XX
          XX[miss.mask.j,jj]  <- X.under[miss.mask.j,j] * con.params$sd[loc] + con.params$mean[loc]
        }
      }
      
      # Binary variables
      if(j %in% bin.ind) {
        loc <- which(bin.ind == j)
        X.under[!miss.mask.j, j] <- truncnorm::rtruncnorm(n = 1,
                                                          a = thresh.bin$true$lower[!miss.mask.j,loc],
                                                          b = thresh.bin$true$upper[!miss.mask.j,loc],
                                                          mean = cond.mean[!miss.mask.j],
                                                          sd = cond.sd)
        
        if(miss.num.j > 0) {
          # Weighted sampling
          temp <- (theta - XX[,-jj] %*% beta[-jj])[miss.mask.j]
          prob.mat <- cbind(exp(-temp^2/2), exp(-(beta[jj] - temp)^2/2))
          value.mat <- cbind(rep(0, miss.num.j), rep(1, miss.num.j))
          weighted.bin <- weighted.sampling(prob.mat, value.mat)
          
          # Set new thresholds
          threshjj <- array(0, dim = c(miss.num.j, 2))
          threshjj[, 1] <- -Inf
          threshjj[, 2] <- Inf
          threshjj[weighted.bin == 1, 1] <- bin.params[loc]
          threshjj[weighted.bin == 0, 2] <- bin.params[loc]
          
          # Update X.under
          X.under[miss.mask.j, j] <- truncnorm::rtruncnorm(n = 1,
                                                           a = threshjj[,1],
                                                           b = threshjj[,2],
                                                           mean = cond.mean[miss.mask.j],
                                                           sd = cond.sd)
          
          # Update XX
          XX[miss.mask.j,jj] <- weighted.bin
        }
      }
      
      # Ordinal variables
      if(j %in% ord.ind) {
        loc <- which(ord.ind == j)
        X.under[!miss.mask.j, j] <- truncnorm::rtruncnorm(n = 1,
                                                          a = thresh.ord$true$lower[!miss.mask.j,loc],
                                                          b = thresh.ord$true$upper[!miss.mask.j,loc],
                                                          mean = cond.mean[!miss.mask.j],
                                                          sd = cond.sd)
        
        if(miss.num.j > 0) {
          # Weighted sampling
          temp <- (theta - XX[,-jj] %*% beta[-jj])[miss.mask.j]
          temp.matr <- outer(rep(1, miss.num.j), c(0,cumsum(beta[jj]))) - outer(temp, rep(1, ord.max[loc]+1))
          prob.mat <- exp(-temp.matr^2/2)
          value.mat <- outer(rep(1, miss.num.j), (1:(ord.max[loc] + 1)))
          weighted.ord <- weighted.sampling(prob.mat, value.mat)
          
          # Set new thresholds
          extend.ord.params.j <- c(-Inf, ord.params[loc,], Inf)
          threshjj <- matrix(0, nrow = miss.num.j, ncol = 2)
          threshjj[,1] <- extend.ord.params.j[weighted.ord]
          threshjj[,2] <- extend.ord.params.j[weighted.ord + 1]
          
          # Update X.under
          X.under[miss.mask.j, j] <- truncnorm::rtruncnorm(n = 1,
                                                           a = threshjj[,1],
                                                           b = threshjj[,2],
                                                           mean = cond.mean[miss.mask.j],
                                                           sd = cond.sd)
          
          # Update XX
          dum <- matrix(0, nrow = miss.num.j, ncol = ord.max[loc])
          ord <- outer(rep(1, miss.num.j), (1:(ord.max[loc]))) <= outer((weighted.ord - 1), rep(1, ord.max[loc]))
          dum[ord] <- 1
          XX[miss.mask.j,jj] <- dum
        }
      }
    }
    
    # Knockoff covariates
    if (j > p) {
      j.ad <- j-p # The adujsted index
      miss.mask.j <- miss.mask[,j.ad]
      miss.num.j <- sum(miss.mask.j)
      
      # Continuous variables
      if(j.ad %in% con.ind) {
        loc <- which(con.ind == j.ad)
        if(miss.num.j > 0) {
          mean.miss <- (cond.sd^2 * beta[jj] * (theta - XX[,-jj] %*% beta[-jj])
                        + cond.mean)[miss.mask.j]/(1 + beta[jj]^2 * cond.sd^2)
          sd.miss <- cond.sd/sqrt(1 + cond.sd^2 * beta[jj]^2)
          
          # Update X.knock.under
          X.knock.under[miss.mask.j,j.ad] <- mean.miss + rnorm(n = miss.num.j, mean = 0, sd = sd.miss)
          
          # Update XX
          XX[miss.mask.j,jj] <- X.knock.under[miss.mask.j,j.ad] * con.params$sd[loc] + con.params$mean[loc]
        }
      }
      
      # Binary variables
      if(j.ad %in% bin.ind) {
        loc <- which(bin.ind == j.ad)
        X.knock.under[!miss.mask.j, j.ad] <- truncnorm::rtruncnorm(n = 1,
                                                                   a = thresh.bin$knock$lower[!miss.mask.j,loc],
                                                                   b = thresh.bin$knock$upper[!miss.mask.j,loc],
                                                                   mean = cond.mean[!miss.mask.j],
                                                                   sd = cond.sd)
        
        if(miss.num.j > 0) {
          # Weighted sampling
          temp <- (theta - XX[,-jj] %*% beta[-jj])[miss.mask.j]
          prob.mat <- cbind(exp(-temp^2/2), exp(-(beta[jj] - temp)^2/2))
          value.mat <- cbind(rep(0, miss.num.j), rep(1, miss.num.j))
          weighted.bin.knock <- weighted.sampling(prob.mat, value.mat)
          
          # Set new indices 
          threshjj <- array(0, dim = c(miss.num.j, 2))
          threshjj[, 1] <- -Inf
          threshjj[, 2] <- Inf
          threshjj[weighted.bin.knock == 1, 1] <- bin.params[loc]
          threshjj[weighted.bin.knock == 0, 2] <- bin.params[loc]
          
          # Update X.knock.under
          X.knock.under[miss.mask.j, j.ad] <- truncnorm::rtruncnorm(n = 1,
                                                                    a = threshjj[, 1],
                                                                    b = threshjj[, 2],
                                                                    mean = cond.mean[miss.mask.j],
                                                                    sd = cond.sd)
          
          # Update XX
          XX[miss.mask.j, jj] <- weighted.bin.knock
        }
      }
      
      # Ordinal variables
      if(j.ad %in% ord.ind) {
        loc <- which(ord.ind == j.ad)
        X.knock.under[!miss.mask.j, j.ad] <- truncnorm::rtruncnorm(n = 1,
                                                                   a = thresh.ord$knock$lower[!miss.mask.j,loc],
                                                                   b = thresh.ord$knock$upper[!miss.mask.j,loc],
                                                                   mean = cond.mean[!miss.mask.j],
                                                                   sd = cond.sd)
        
        if(miss.num.j > 0) {
          # Weighted sampling
          temp <- (theta - XX[,-jj] %*% beta[-jj])[miss.mask.j]
          temp.matr <- outer(rep(1, miss.num.j), c(0,cumsum(beta[jj]))) - outer(temp, rep(1, ord.max[loc]+1))
          prob.mat <- exp(-temp.matr^2/2)
          value.mat <- outer(rep(1, miss.num.j), (1:(ord.max[loc] + 1)))
          weighted.knock.ord <- weighted.sampling(prob.mat, value.mat)
          
          # Set new thresholds
          extend.ord.params.j <- c(-Inf, ord.params[loc,], Inf)
          threshjj <- matrix(0, nrow = miss.num.j, ncol = 2)
          threshjj[,1] <- extend.ord.params.j[weighted.knock.ord]
          threshjj[,2] <- extend.ord.params.j[weighted.knock.ord + 1]
          
          # Update X.knock.under
          X.knock.under[miss.mask.j, j.ad] <- truncnorm::rtruncnorm(n = 1,
                                                                    a = threshjj[,1],
                                                                    b = threshjj[,2],
                                                                    mean = cond.mean[miss.mask.j],
                                                                    sd = cond.sd)
          
          # Update XX
          dum.knock <- matrix(0, nrow = miss.num.j, ncol = ord.max[loc])
          ord <- outer(rep(1, miss.num.j), (1:(ord.max[loc]))) <= outer((weighted.knock.ord - 1), rep(1, ord.max[loc]))
          dum.knock[ord] <- 1
          XX[miss.mask.j,jj] <- dum.knock
        }
      }
    }
  }
  
  results <- list(X.under = X.under, XX = XX)
  if (is.knockoff) {
    results$X.knock.under <- X.knock.under
  }
  
  return(results)
}

# Initialize underlying variables
X.under.init <- function(X, bin.ind=NULL, ord.ind=NULL, con.ind=NULL, 
                         thresh.bin=NULL, thresh.ord=NULL, con.params=NULL) {
  X.under <-  matrix(nrow = dim(X)[1], ncol = dim(X)[2])
  
  if (length(bin.ind) + length(ord.ind) > 0) {
    if (!requireNamespace('truncnorm', quietly=T)) {
      stop('truncnorm is not installed', call.=F)
    }  
  }
  
  if (length(bin.ind) > 0) {
    for (j in 1:length(bin.ind)) {
      bj <- bin.ind[j]
      X.under[, bj] <- truncnorm::rtruncnorm(a = thresh.bin$lower[,j],
                                             b = thresh.bin$upper[,j],
                                             n = 1)
    }
  }
  
  if (length(ord.ind) > 0) {
    for (j in 1:length(ord.ind)) {
      oj <- ord.ind[j]
      X.under[, oj] <- truncnorm::rtruncnorm(a = thresh.ord$lower[,j],
                                             b = thresh.ord$upper[,j],
                                             n = 1)
    }
  }
  
  if (length(con.ind > 0)) {
    for (j in 1:length(con.ind)) {
      cj <- con.ind[j]
      miss.mask.j <- is.na(X[,cj])
      X.under[!miss.mask.j, cj] <- (X[!miss.mask.j, cj] - con.params$mean[j])/con.params$sd[j]
      X.under[miss.mask.j, cj] <- rnorm(n = sum(miss.mask.j))
    }
  }
  
  return(X.under)
}

# Stochastic EM algorithm
StEM <- function(X, X.under=NULL, X.knock=NULL, X.knock.under=NULL, Y, Sigma, S=NULL, a.vec, d1.vec, d2.vec,
                 bin.ind=NULL, ord.ind=NULL, con.ind=NULL, bin.params=NULL, 
                 ord.params=NULL, con.params=NULL, eps=1e-2, max_iter=1000, 
                 burn_in=500, is.beta.trace=F, is.print.iter=F, 
                 is.knockoff=F, ...) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  J <- dim(Y)[2]
  miss.mask <- is.na(X)
  args <- list(...)
  
  if (is.knockoff) {
    pp <- 2*p
    # The augmented correlation matrix
    Sigma.aug <- cbind(rbind(Sigma, Sigma - S), rbind(Sigma - S, Sigma))
    # The inverse of augmented correlation matrix
    Omega <- solve(Sigma.aug)
  } else {
    pp <- p
    # Inverse matrix
    Omega <- solve(Sigma)
  }
  
  if (is.null(X.under)) {
    # Randomly initialize underlying values if initialization is not provided.
    X.under <- X.under.init(X, bin.ind, ord.ind, con.ind, 
                            thresh.bin, thresh.ord, con.params)
  }
  
  # Check if there are polytomous items 
  Y.max <- apply(Y, 2, max, na.rm = T)
  dich.ind <- which(Y.max == 1)
  poly.ind <- which(Y.max > 1)
  
  # Initial values
  if (!requireNamespace('mirt', quietly=T)) {
    stop('mirt is not installed', call.=F)
  }
  
  if (!is.null(args$theta0)) {
    theta0 <- args$theta0
  } else {
    # Use mirt to estimate initial values of theta
    if (length(poly.ind) > 0) {
      # TODO: Extend to polytomous items with more than 3 categories
      Y.model<- mirt::mirt(as.data.frame(Y), 1, itemtype= 'gpcm', SE = F, pars = 'values')
      Y.model$value[Y.model$name == 'a1'] <- a.vec
      Y.model$value[Y.model$name == 'd1'] <- d1.vec
      Y.model$value[Y.model$name == 'd2'] <- (d1.vec + d2.vec)[which(!is.na(d2.vec))]
      Y.model$est <- FALSE
      theta0 <- mirt::fscores(mirt::mirt(as.data.frame(Y), 1, itemtype= 'gpcm', pars=Y.model), rotate = F)  
    } else {
      Y.model <- mirt::mirt(as.data.frame(Y), 1, itemtype= '2PL', SE = F, pars = 'values')
      Y.model$value[Y.model$name == 'a1'] <- a.vec
      Y.model$value[Y.model$name == 'd1'] <- d1.vec
      Y.model$est<-FALSE
      theta0 <- mirt::fscores(mirt::mirt(as.data.frame(Y), 1, pars=Y.model), rotate = F)
    }
  }
  
  if(!is.null(bin.ind)) {
    # Thresholds for binary variables
    thresh.bin <- list()
    thresh.bin$true <- list(lower = matrix(-Inf, nrow = N, ncol = length(bin.ind)),
                            upper = matrix(Inf, nrow = N, ncol = length(bin.ind)))
    if (is.knockoff) {
      thresh.bin$knock <- list(lower = matrix(-Inf, nrow = N, ncol = length(bin.ind)),
                               upper = matrix(Inf, nrow = N, ncol = length(bin.ind))) 
    }
    
    for (j in 1:length(bin.ind)) {
      bj <- bin.ind[j]
      thresh.bin$true$lower[X[,bj] == 1, j] <- bin.params[j]
      thresh.bin$true$upper[X[,bj] == 0, j] <- bin.params[j]
      if (is.knockoff) {
        thresh.bin$knock$lower[X.knock[,bj] == 1, j] <- bin.params[j]
        thresh.bin$knock$upper[X.knock[,bj] == 0, j] <- bin.params[j]  
      }
    } 
  }
  
  ord.max <- NULL
  if(!is.null(ord.ind)) {
    ord.max <- apply(X[,ord.ind], 2, max, na.rm = T)
    thresh.ord <- list()
    thresh.ord$true <- list(lower = matrix(-Inf, nrow = N, ncol = length(ord.ind)),
                            upper = matrix(Inf, nrow = N, ncol = length(ord.ind)))
    if (is.knockoff) {
      thresh.ord$knock <- list(lower = matrix(-Inf, nrow = N, ncol = length(ord.ind)),
                               upper = matrix(Inf, nrow = N, ncol = length(ord.ind))) 
    }
    
    for (j in 1:length(ord.ind)) {
      oj <- ord.ind[j]
      for (v in 1:(ord.max[j]+1)) {
        if (v > 1) {
          thresh.ord$true$lower[which(X[,oj] == v), j] <- ord.params[j, v-1]
          if (is.knockoff) {
            thresh.ord$knock$lower[which(X.knock[,oj] == v), j] <- ord.params[j, v-1]
          }
        }
        if (v < (ord.max[j]+1)) {
          thresh.ord$true$upper[which(X[,oj] == v), j] <- ord.params[j, v]
          if (is.knockoff) {
            thresh.ord$knock$upper[which(X.knock[,oj] == v), j] <- ord.params[j, v]
          }
        }
      }
    }
  }
  
  # Impute true and knockoff covariates
  X.impute <- Impute(X, X.under, bin.ind, ord.ind, con.ind, bin.params, ord.params, con.params)
  if (is.knockoff) {
    X.knock.impute <- Impute(X.knock, X.knock.under, bin.ind, ord.ind, con.ind, bin.params, ord.params, con.params)  
  }
  
  # Make names of X
  if (!is.null(names(X))) {
    X.names <- names(X)
  } else {
    X.names <- paste0('X', 1:p)
  }
  
  # Combine all covariates and record the indices
  XX <- NULL
  XX.names <- NULL
  XX.ind <- list()
  len.temp <- 0
  for (j in 1:pp) {
    if (j <= p) {
      # Add true contiunous or binary variables into XX
      if (!(j %in% ord.ind)) {
        len.temp <- len.temp + 1
        XX <- cbind(XX, X.impute[,j])
        XX.ind[[j]] <- len.temp
        XX.names <- c(XX.names, X.names[j])
      }
      # Add true ordinal dummy variables into XX
      if (j %in% ord.ind) {
        loc <- which(ord.ind == j)
        XX <- cbind(XX, ord.to.dum(X.impute[,j]))
        XX.ind[[j]] <- len.temp + 1:ord.max[loc]
        len.temp <- len.temp + ord.max[loc]
        XX.names <- c(XX.names, paste0(paste0(X.names[j], '.dum'), 1:ord.max[loc]))
      }
    }
    
    if (j > p) {
      #Adjusted index
      j.ad <- j-p
      # Add knockoffs contiunous or binary variables into XX
      if (!(j.ad %in% ord.ind)) {
        len.temp <- len.temp + 1
        XX <- cbind(XX, X.knock.impute[,j.ad])
        XX.ind[[j]] <- len.temp
        XX.names <- c(XX.names, paste0(X.names[j.ad], '.knock'))
      }
      # Add knockoffs ordinal dummy variables into XX
      if ((j.ad %in% ord.ind)) {
        loc <- which(ord.ind == j.ad)
        XX <- cbind(XX, ord.to.dum(X.knock.impute[,j.ad]))
        XX.ind[[j]] <- len.temp + 1:ord.max[loc]
        len.temp <- len.temp + ord.max[loc]
        XX.names <- c(XX.names, paste0(paste0(X.names[j.ad], '.knock.dum'), 1:ord.max[loc]))
      }
    }
  }
  
  # Names of XX
  XX <- as.data.frame(XX)
  names(XX) <- XX.names
  
  # Initialize theta and parameters
  theta <- c(theta0)
  lm.model <- lm(theta ~ ., data = as.data.frame(cbind(theta, XX)))
  coeff <- coef(lm.model)
  
  if (!is.null(args$beta0)) {
    beta0 <- args$beta0
  } else {
    beta0 <- coeff[-1]
  }
  
  if (!is.null(args$theta.interp0)) {
    theta.interp0 <- args$theta.interp0
  } else {
    theta.interp0 <- coeff[1]
  }
  
  if (!is.null(args$theta.var0)) {
    theta.var0 <- args$theta.var0
  } else {
    theta.var0 <- mean(lm.model$residuals^2)
  }
  
  beta <- beta0
  theta.interp <- theta.interp0
  theta.var <- theta.var0
  
  # Record of results
  beta.avg <- 0
  if (is.beta.trace) {
    beta.record <- NULL
  }
  theta.interp.avg <- 0
  theta.var.avg <- 0
  
  # There are 2*p columns if there is knockoff copy
  meantrans.mat <- matrix(0, nrow = pp-1 , ncol = pp)
  for (i in 1:pp) {
    meantrans.mat[, i] <- -Omega[i, -i]/Omega[i, i]
  }
  
  # Stochastic EM algorithm
  for(iter in 1:max_iter) {
    if (is.print.iter) {
      print(iter)
    }
    
    #===============#
    # Sampling Step #
    #===============#
    # Sampling theta
    theta <- theta.sample(XX, Y, beta, theta.interp, theta.var, a.vec, d1.vec, d2.vec)
    
    # Sampling X.under, X.knock.under and XX
    res.sample <- XX.sample(XX, X, X.under, X.knock.under, theta, beta, Omega, 
                            XX.ind, bin.params, ord.params, con.params, bin.ind, 
                            ord.ind, con.ind, thresh.bin, thresh.ord, is.knockoff = is.knockoff,
                            miss.mask = miss.mask, meantrans.mat = meantrans.mat,
                            ord.max = ord.max)
    
    X.under <- res.sample$X.under
    XX <- res.sample$XX
    if (is.knockoff) {
      X.knock.under <- res.sample$X.knock.under 
    }
    
    #===================#
    # Maximization Step #
    #===================#
    # Optimization by simple linear regression
    lm.model <- lm(theta ~ ., data = as.data.frame(cbind(theta, XX)))
    coeffs <- coef(lm.model)
    beta <- coeffs[-1]
    theta.interp <- coeffs[1]
    theta.var <- mean(lm.model$residuals^2)
    
    if (is.beta.trace) {
      beta.record <- cbind(beta.record, beta)
    }
    
    if(iter > burn_in) {
      beta.avg <- ((iter-burn_in-1)*beta.avg + beta)/(iter-burn_in)
      theta.interp.avg <- ((iter-burn_in-1)*theta.interp.avg + theta.interp)/(iter-burn_in)
      theta.var.avg <- ((iter-burn_in-1)*theta.var.avg + theta.var)/(iter-burn_in)
    }
    
    # TODO: Add stopping criterion
  }
  
  # Save results
  result <- list()
  result$iter <- iter
  setNames(beta.avg, names(XX))
  result$beta <- beta.avg
  result$theta <- theta
  result$interp <- theta.interp.avg
  result$var <- theta.var.avg
  result$XX <- XX
  result$XX.ind <- XX.ind
  result$X.under <- X.under
  if (is.beta.trace) {
    result$beta.record <- beta.record
  }
  
  return(result)
}

X.Gibbs.sample <- function(X, X.under, XX, XX.ind, Y, Omega, a.vec, d1.vec, d2.vec,
                           beta, theta.interp, theta.var,
                           bin.ind=NULL, ord.ind=NULL, con.ind=NULL, bin.params=NULL, 
                           ord.params=NULL, con.params=NULL, max_iter=1000) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  J <- dim(Y)[2]
  miss.mask <- is.na(X)
  
  if(!is.null(bin.ind)) {
    # Thresholds for binary variables
    thresh.bin <- list()
    thresh.bin$true <- list(lower = matrix(-Inf, nrow = N, ncol = length(bin.ind)),
                            upper = matrix(Inf, nrow = N, ncol = length(bin.ind)))
    
    for (j in 1:length(bin.ind)) {
      bj <- bin.ind[j]
      thresh.bin$true$lower[X[,bj] == 1, j] <- bin.params[j]
      thresh.bin$true$upper[X[,bj] == 0, j] <- bin.params[j]
    }
  }
  
  ord.max <- NULL
  if(!is.null(ord.ind)) {
    ord.max <- apply(X[,ord.ind], 2, max, na.rm = T)
    thresh.ord <- list()
    thresh.ord$true <- list(lower = matrix(-Inf, nrow = N, ncol = length(ord.ind)),
                            upper = matrix(Inf, nrow = N, ncol = length(ord.ind)))
    
    for (j in 1:length(ord.ind)) {
      oj <- ord.ind[j]
      for (v in 1:(ord.max[j]+1)) {
        if (v > 1) {
          thresh.ord$true$lower[which(X[,oj] == v), j] <- ord.params[j, v-1]
        }
        if (v < (ord.max[j]+1)) {
          thresh.ord$true$upper[which(X[,oj] == v), j] <- ord.params[j, v]
        }
      }
    }
  }
  
  meantrans.mat <- matrix(0, nrow = p-1 , ncol = p)
  for (i in 1:p) {
    meantrans.mat[, i] <- -Omega[i, -i]/Omega[i, i]
  }
  
  for (iter in 1:max_iter) {
    # Sampling theta
    theta <- theta.sample(XX, Y, beta, theta.interp, theta.var, a.vec, d1.vec, d2.vec)
    
    # Sampling X.under, X.knock.under and XX
    res.sample <- XX.sample(XX, X, X.under, NULL, theta, beta, Omega,
                            XX.ind, bin.params, ord.params, con.params, bin.ind,
                            ord.ind, con.ind, thresh.bin, thresh.ord, is.knockoff=F,
                            miss.mask = miss.mask, meantrans.mat = meantrans.mat,
                            ord.max = ord.max)
    
    X.under <- res.sample$X.under
    XX <- res.sample$XX
  }
  
  result <- list()
  result$X.under.Gibbs <- X.under
  result$XX.Gibbs <- XX
  
  return(result)
}