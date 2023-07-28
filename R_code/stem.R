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

# Impute mixed data
Impute <- function(X, X.under, bin.ind=NULL, ord.ind=NULL, con.ind=NULL,
                   bin.params=NULL, ord.params=NULL, con.params=NULL) {
  
  X.impute <- X
  N <- dim(X)[1]
  miss.mask <- is.na(X)
  
  #Impute binary variables
  if(!is.null(bin.ind)) {
    X.impute[, bin.ind] <- (X.under[,bin.ind] > outer(rep(1, N), bin.params)) * 1
  }
  
  #Impute ordinal variables
  if(!is.null(ord.ind)) {
    for (j in ord.ind) {
      loc <- which(ord.ind == j)
      miss.mask.j <- miss.mask[,j]
      X.impute[miss.mask.j, j]  <- rowSums(outer(X.under[miss.mask.j,j], rep(1, dim(ord.params)[2]))
                                           > outer(rep(1, sum(miss.mask.j)), ord.params[loc,]))
    }
  }
  
  #Impute continuous variables
  if(!is.null(con.ind)) {
    X.impute[,con.ind] <- t(t(X.under[,con.ind]) * c(con.params$sd) + c(con.params$mean))
  }
  
  return(X.impute)
}

# Sampling the latent trait theta by adaptive rejection sampling
theta.sample <- function(XX, Y, beta, theta.interp, theta.var, a.vec, d1.vec, d2.vec) {
  Y <- as.matrix(Y)
  XX <- as.matrix(XX)
  mu <- XX %*% beta + theta.interp
  theta.sampler <- arstheta::ThetaPosterior$new(Y=Y, a.vec=a.vec, d1.vec=d1.vec, 
                                                d2.vec=d2.vec, mu=mu, theta.var=theta.var)
  
  theta.ars <- theta.sampler$sample(mu);
  
  return(theta.ars)
}

# Given theta, sample the underlying variables (including or not including knockoffs) 
# and impute the standardized covariate matrix 
XX.sample <- function(XX, XX.ind, X, X.under, theta, beta, theta.interp, theta.var, Omega, 
                      bin.params=NULL, ord.params=NULL, con.params=NULL, 
                      bin.ind=NULL, ord.ind=NULL, con.ind=NULL, 
                      thresh.bin=NULL, thresh.ord=NULL, is.knockoff=F, X.knock.under=NULL,
                      W.list=NULL, W.list.knock=NULL, mea.ind=NULL, mea.params=NULL, 
                      is.random.scan=T, ...) {
  N <- nrow(X)
  p <- ncol(X.under)
  pp <- nrow(Omega)
  
  if(!is.null(W.list)) {
    num.mea <- length(W.list)
    p <- ncol(X) + num.mea
    stopifnot(length(mea.ind) == num.mea)
    stopifnot(!is.null(mea.params))
  }
  
  if(is.knockoff) {
    stopifnot(all(!is.null(X.knock.under),
                  dim(X.knock.under) == dim(X.under)))
    if(!is.null(W.list)) {
      stopifnot(all(!is.null(W.list.knock)))
    }
    stopifnot(pp == 2*p)
  }
  
  # The transformed covariates
  XX <- as.matrix(XX) 
  stopifnot(all(length(unlist(XX.ind)) == ncol(XX),
                length(beta) == ncol(XX)))
  
  # Allowed to provide quantities to save computation
  args <- list(...)
  if (is.null(args$meantrans.mat)) {
    meantrans.mat <- matrix(0, nrow = pp-1 , ncol = pp)
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
      ord.max <- apply(X[,ord.ind], 2, max, na.rm = T)
    } else {
      ord.max <- args$ord.max
    }
  }

  # One-round Gibbs sampling with/without random scanning
  if (is.random.scan) {
    scan <- sample(1:pp)
  } else {
    scan <- 1:pp
  }

  for (j in scan) {
    cond.mean <- cbind(X.under, X.knock.under)[, -j] %*% meantrans.mat[, j]
    cond.sd <- sqrt(1/Omega[j, j])
    jj <- XX.ind[[j]]
    
    # True covariates
    if (j <= p) {
      # Unobserved variables
      if(j %in% mea.ind) {
        loc <- which(mea.ind == j)
        W <- W.list[[loc]]
        W.loadings <- mea.params$loadings[[loc]]
        W.errvar <- mea.params$errvar[[loc]]
        W.mean <- mea.params$mean[[loc]]
        
        post.sd <- 1/sqrt(sum(W.loadings**2/W.errvar)
                          + 1/cond.sd**2
                          + beta[jj]**2/theta.var)
        
        post.mean <- (rowSums((W - outer(rep(1, N), W.mean))
                              * outer(rep(1, N), W.loadings/W.errvar))
                      + cond.mean/cond.sd**2
                      + ((theta - XX[,-jj] %*% beta[-jj] - theta.interp) * beta[jj])/theta.var
                      ) * post.sd**2
        
        # Update X.under
        X.under[,j] <- rnorm(n = N,
                             mean = post.mean,
                             sd = post.sd)
        
        # Update XX
        XX[,jj] <- X.under[,j]
        
      } else {
        miss.mask.j <- miss.mask[,j]
        miss.num.j <- sum(miss.mask.j)
        
        # Continuous variables
        if(j %in% con.ind) {
          if(miss.num.j > 0) {
            loc <- which(con.ind == j)
            if (miss.num.j > 0){
              sd.miss <- 1/sqrt(1/cond.sd**2 + (beta[jj] * con.params$sd[loc])**2/theta.var)
              mean.miss <- (beta[jj] * con.params$sd[loc] * 
                              (theta
                               - XX[,-jj] %*% beta[-jj]
                               - theta.interp
                               - beta[jj] * con.params$mean[loc])/theta.var
                            + cond.mean/cond.sd**2)[miss.mask.j] * sd.miss**2
              
              # Update X.under
              X.under[miss.mask.j,j] <- mean.miss + rnorm(n = miss.num.j, mean = 0, sd = sd.miss)
              
              # Update XX
              XX[miss.mask.j,jj]  <- X.under[miss.mask.j,j] * con.params$sd[loc] + con.params$mean[loc]
            }
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
            temp <- (theta - XX[,-jj] %*% beta[-jj] - theta.interp)[miss.mask.j]
            temp.mat <- cbind(exp(-temp^2/(2*theta.var)), exp(-(beta[jj] - temp)^2/(2*theta.var)))
            prob.truc.mat <- matrixStats::rowDiffs(cbind(rep(0, miss.num.j),
                                                         pnorm(q = bin.params[loc],
                                                               mean = cond.mean[miss.mask.j],
                                                               sd = cond.sd),
                                                         rep(1, miss.num.j)))
            prob.mat <- temp.mat * prob.truc.mat
            prob.mat <- pmax(pmin(prob.mat, 1-1e-10), 1e-10)
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
            temp <- (theta - XX[,-jj] %*% beta[-jj] - theta.interp)[miss.mask.j]
            temp.mat <- outer(rep(1, miss.num.j), c(0,cumsum(beta[jj]))) - outer(temp, rep(1, ord.max[loc]+1))
            temp.mat <- exp(-temp.mat^2/(2*theta.var))
            prob.truc.mat <- matrixStats::rowDiffs(cbind(rep(0, miss.num.j),
                                                         pnorm(q = outer(rep(1, miss.num.j), ord.params[loc, 1:ord.max[loc]]),
                                                               mean = outer(cond.mean[miss.mask.j], rep(1, ord.max[loc])),
                                                               sd = cond.sd),
                                                         rep(1, miss.num.j)))
            prob.mat <- exp(log(prob.truc.mat) + log(temp.mat))
            prob.mat <- pmax(pmin(prob.mat, 1-1e-10), 1e-10)
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
    }
    
    # Knockoff covariates
    if (j > p) {
      j.ad <- j-p # The adujsted index
      
      # Unobserved variables
      if(j.ad %in% mea.ind) {
        loc <- which(mea.ind == j.ad)
        W <- W.list[[loc]]
        W.loadings <- mea.params$loadings[[loc]]
        W.errvar <- mea.params$errvar[[loc]]
        W.mean <- mea.params$mean[[loc]]
        
        post.sd <- 1/sqrt(sum(W.loadings**2/W.errvar)
                          + 1/cond.sd**2
                          + beta[jj]**2/theta.var)
        
        post.mean <- (rowSums((W - outer(rep(1, N), W.mean))
                              * outer(rep(1, N), W.loadings/W.errvar))
                      + cond.mean/cond.sd**2
                      + ((theta - XX[,-jj] %*% beta[-jj] - theta.interp) * beta[jj])/theta.var
                      ) * post.sd**2
        
        # Update X.under
        X.knock.under[,j.ad] <- rnorm(n = N, 
                                      mean = post.mean, 
                                      sd = post.sd)
        
        # Update XX
        XX[,jj] <- X.knock.under[,j.ad]
        
      } else {
        miss.mask.j <- miss.mask[,j.ad]
        miss.num.j <- sum(miss.mask.j)
        
        # Continuous variables
        if(j.ad %in% con.ind) {
          loc <- which(con.ind == j.ad)
          if(miss.num.j > 0) {
            sd.miss <- 1/sqrt(1/cond.sd**2 + beta[jj]^2/theta.var)
            mean.miss <- (beta[jj] * con.params$sd[loc] * 
                            (theta
                             - XX[,-jj] %*% beta[-jj]
                             - theta.interp
                             - beta[jj] * con.params$mean[loc])/theta.var
                          + cond.mean/cond.sd**2)[miss.mask.j] * sd.miss**2
            
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
            temp <- (theta - XX[,-jj] %*% beta[-jj] - theta.interp)[miss.mask.j]
            temp.mat <- cbind(exp(-temp^2/(2*theta.var)), exp(-(beta[jj] - temp)^2/(2*theta.var)))
            prob.truc.mat <- matrixStats::rowDiffs(cbind(rep(0, miss.num.j),
                                                         pnorm(q = bin.params[loc],
                                                               mean = cond.mean[miss.mask.j],
                                                               sd = cond.sd),
                                                         rep(1, miss.num.j)))
            prob.mat <- temp.mat * prob.truc.mat
            prob.mat <- pmax(pmin(prob.mat, 1-1e-10), 1e-10)
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
            temp <- (theta - XX[,-jj] %*% beta[-jj] - theta.interp)[miss.mask.j]
            temp.mat <- outer(rep(1, miss.num.j), c(0,cumsum(beta[jj]))) - outer(temp, rep(1, ord.max[loc]+1))
            temp.mat <- exp(-temp.mat^2/(2*theta.var))
            prob.truc.mat <- matrixStats::rowDiffs(cbind(rep(0, miss.num.j),
                                                         pnorm(q = outer(rep(1, miss.num.j), ord.params[loc, 1:ord.max[loc]]),
                                                               mean = outer(cond.mean[miss.mask.j], rep(1, ord.max[loc])),
                                                               sd = cond.sd),
                                                         rep(1, miss.num.j)))
            prob.mat <- exp(log(prob.truc.mat) + log(temp.mat))
            prob.mat <- pmax(pmin(prob.mat, 1-1e-10), 1e-10)
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
  }
  
  results <- list(X.under = X.under, XX = XX)
  if (is.knockoff) {
    results$X.knock.under <- X.knock.under
  }
  
  return(results)
}

# Given Y, sample theta and underlying variables (not including knockoffs) 
# and impute the standardized covariate matrix.
# Used in generating knockoffs.
X.Gibbs.sample <- function(X, X.under, XX, XX.ind, Y, Omega, a.vec, d1.vec, d2.vec,
                           beta, theta.interp, theta.var,
                           bin.ind=NULL, ord.ind=NULL, con.ind=NULL, 
                           bin.params=NULL, ord.params=NULL, con.params=NULL, 
                           W.list=NULL, mea.ind=NULL, mea.params=NULL,
                           max_iter=1000) {
  N <- nrow(X)
  J <- ncol(Y)
  
  if (!is.null(W.list)) {
    num.mea <- length(W.list)
    p <- ncol(X) + num.mea
    stopifnot(length(mea.ind) == num.mea)
    stopifnot(!is.null(mea.params))
  } else {
    p <- ncol(X)
  }
  stopifnot(all(nrow(Omega) == p, ncol(X.under) == p))
  stopifnot(all(length(unlist(XX.ind)) == ncol(XX),
                length(beta) == ncol(XX)))
  
  miss.mask <- is.na(X)
  # Thresholds for binary variables
  thresh.bin <- NULL
  if(!is.null(bin.ind)) {
    thresh.bin <- list()
    thresh.bin$true <- list(lower = matrix(-Inf, nrow = N, ncol = length(bin.ind)),
                            upper = matrix(Inf, nrow = N, ncol = length(bin.ind)))
    
    for (j in 1:length(bin.ind)) {
      bj <- bin.ind[j]
      thresh.bin$true$lower[X[,bj] == 1, j] <- bin.params[j]
      thresh.bin$true$upper[X[,bj] == 0, j] <- bin.params[j]
    }
  }
  
  # Thresholds for ordinal variables
  ord.max <- NULL
  thresh.ord <- NULL
  if(!is.null(ord.ind)) {
    ord.max <- apply(X[,ord.ind], 2, max, na.rm = T)
    thresh.ord <- list()
    thresh.ord$true <- list(lower = matrix(-Inf, nrow = N, ncol = length(ord.ind)),
                            upper = matrix(Inf, nrow = N, ncol = length(ord.ind)))
    
    for (j in 1:length(ord.ind)) {
      oj <- ord.ind[j]
      for (v in (0:ord.max[j])) {
        if (v > 0) {
          thresh.ord$true$lower[which(X[,oj] == v), j] <- ord.params[j, v]
        }
        if (v < ord.max[j]) {
          thresh.ord$true$upper[which(X[,oj] == v), j] <- ord.params[j, v+1]
        }
      }
    }
  }
  
  # Mean transformation matrix
  meantrans.mat <- matrix(0, nrow = p-1 , ncol = p)
  for (i in 1:p) {
    meantrans.mat[, i] <- -Omega[i, -i]/Omega[i, i]
  }
  
  # Gibbs sampler
  for (iter in 1:max_iter) {
    print(iter)
    
    # Sampling theta
    theta <- theta.sample(XX=XX, Y=Y, beta=beta, theta.interp=theta.interp, theta.var=theta.var, 
                          a.vec=a.vec, d1.vec=d1.vec, d2.vec=d2.vec)
    
    # Sampling X.under and XX
    res.sample <- XX.sample(XX=XX, XX.ind=XX.ind, X=X, X.under=X.under, theta=theta, 
                            beta=beta, theta.interp=theta.interp, theta.var=theta.var, Omega=Omega, 
                            bin.params=bin.params, ord.params=ord.params, con.params=con.params, 
                            bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind, 
                            thresh.bin=thresh.bin, thresh.ord=thresh.ord, is.knockoff=F, 
                            X.knock.under=NULL, W.list=W.list, W.list.knock=NULL, mea.ind=mea.ind, 
                            mea.params=mea.params, miss.mask=miss.mask, meantrans.mat=meantrans.mat,
                            ord.max = ord.max)
    
    X.under <- res.sample$X.under
    XX <- res.sample$XX
  }
  
  result <- list()
  result$X.under.Gibbs <- X.under
  result$XX.Gibbs <- XX
  result$theta.Gibbs <- theta
  
  return(result)
}


# Stochastic EM algorithm
StEM <- function(X, X.under, X.knock=NULL, X.knock.under=NULL, Y, Sigma, S=NULL, 
                 a.vec=NULL, d1.vec=NULL, d2.vec=NULL,
                 bin.ind=NULL, ord.ind=NULL, con.ind=NULL, bin.params=NULL, 
                 ord.params=NULL, con.params=NULL, eps=1e-2, max_iter=1000, 
                 burn_in=500, is.beta.trace=F, is.print.iter=F, is.knockoff=F,
                 W.list=NULL, mea.params=NULL, W.list.knock=NULL, lambda=0.1,
                 ...) {
  N <- nrow(X)
  J <- ncol(Y)
  args <- list(...)
  miss.mask <- is.na(X)
  
  # Include unobserved predictors
  if (!is.null(W.list)) {
    num.mea <- length(W.list)
    p <- ncol(X) + num.mea
    mea.ind <- length(bin.ind) + length(ord.ind) + length(con.ind) + 1:num.mea
    stopifnot(!is.null(mea.params))
  } else {
    p <- ncol(X)
    mea.ind <- NULL
  }
  
  # Check dimensions
  stopifnot(ncol(X.under) == p)
  
  # Esimate IRT model or not
  is.estimateIRT <- F
  if (any(c(is.null(a.vec), is.null(d1.vec), is.null(d2.vec)))) {
    is.estimateIRT <- T
  }
  
  # Knockoff or not
  if (is.knockoff) {
    stopifnot(all(!is.null(X.knock),
                  !is.null(X.knock.under),
                  !is.null(S)))
    
    if (!is.null(W.list)) {
      stopifnot(!is.null(W.list.knock))
    }
  }
  
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
    # Use mirt to estimate initial values of theta. 
    
    if (!is.estimateIRT) {
      # Do not estimate IRT parameters.
      if (length(poly.ind) > 0) {
        # With polytomous items
        # TODO: Extend to polytomous items with more than 3 categories
        Y.model<- mirt::mirt(as.data.frame(Y), 1, itemtype= 'gpcm', SE = F, pars = 'values')
        Y.model$value[Y.model$name == 'a1'] <- a.vec
        Y.model$value[Y.model$name == 'd1'] <- d1.vec
        Y.model$value[Y.model$name == 'd2'] <- (d1.vec + d2.vec)[which(!is.na(d2.vec))]
        Y.model$est <- FALSE
        Y.model$est[which(Y.model$name == 'MEAN_1')] <- T
        Y.model$est[which(Y.model$name == 'COV_11')] <- T
        theta0 <- mirt::fscores(mirt::mirt(as.data.frame(Y), 1, itemtype= 'gpcm', pars=Y.model), rotate = F)
        
      } else {
        # Without polytomous items
        Y.model <- mirt::mirt(as.data.frame(Y), 1, itemtype= '2PL', SE = F, pars = 'values')
        Y.model$value[Y.model$name == 'a1'] <- a.vec
        Y.model$value[Y.model$name == 'd'] <- d1.vec
        Y.model$est<-FALSE
        Y.model$est[which(Y.model$name == 'MEAN_1')] <- T
        Y.model$est[which(Y.model$name == 'COV_11')] <- T
        
        theta0 <- mirt::fscores(mirt::mirt(as.data.frame(Y), 1, pars=Y.model), rotate = F)
      }
    } else{
      if (length(poly.ind) > 0) {
        # TODO: Extend to polytomous items with more than 3 categories
        Y.model <- mirt::mirt(as.data.frame(Y), 1, itemtype= 'gpcm', SE = F, pars = 'values')
        Y.model$est[which(Y.model$name == 'MEAN_1')] <- F
        Y.model$est[which(Y.model$name == 'COV_11')] <- F
        
        # Estimate GPCM model
        Y.model <- mirt::mirt(as.data.frame(Y), 1, itemtype= 'gpcm', SE = F, pars = Y.model)
        IRTparams <- as.data.frame(mirt::coef(Y.model, simplify = TRUE)$items)
        a.vec <- IRTparams$a1
        d1.vec <- IRTparams$d1
        d2.vec <- IRTparams$d2
        theta0 <- mirt::fscores(Y.model, rotate = F)
        
      } else {
        Y.model <- mirt::mirt(as.data.frame(Y), 1, itemtype= '2PL', SE = F, pars = 'values')
        Y.model$est[which(Y.model$name == 'MEAN_1')] <- F
        Y.model$est[which(Y.model$name == 'COV_11')] <- F
        
        # Estimate 2PL model
        Y.model <- mirt::mirt(as.data.frame(Y), 1, itemtype= '2PL', SE = F, pars = Y.model)
        IRTparams <- as.data.frame(mirt::coef(Y.model, simplify = TRUE)$items)
        a.vec <- IRTparams$a1
        d1.vec <- IRTparams$d
        d2.vec <- rep(NA, length(d1.vec))
        theta0 <- mirt::fscores(Y.model, rotate = F)
      }
    }
  }
  
  # Thresholds for binary variables
  thresh.bin <- NULL
  if(!is.null(bin.ind)) {
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
  
  # Thresholds for ordinal variables
  ord.max <- NULL
  thresh.ord <- NULL
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
      for (v in (0:ord.max[j])) {
        if (v > 0) {
          thresh.ord$true$lower[which(X[,oj] == v), j] <- ord.params[j, v]
          if (is.knockoff) {
            thresh.ord$knock$lower[which(X.knock[,oj] == v), j] <- ord.params[j, v]
          }
        }
        if (v < ord.max[j]) {
          thresh.ord$true$upper[which(X[,oj] == v), j] <- ord.params[j, v+1]
          if (is.knockoff) {
            thresh.ord$knock$upper[which(X.knock[,oj] == v), j] <- ord.params[j, v+1]
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
    if (!is.null(W.list)) {
      X.names <- c(X.names, paste0('X.mea.', 1:num.mea))
    }
  } else {
    X.names <- paste0('X', 1:p)
  }
  stopifnot(length(X.names) == p)
  
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
        if (j %in% mea.ind) {
          # Add the underlying variable if predictor unobserved
          XX <- cbind(XX, X.under[,j])
        } else {
          XX <- cbind(XX, X.impute[,j])
        }
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
        if (j.ad %in% mea.ind) {
          # Add the underlying variable if predictor unobserved
          XX <- cbind(XX, X.knock.under[,j.ad])
        } else {
          XX <- cbind(XX, X.knock.impute[,j.ad])
        }
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
  
  # Variable indicating whether initial values of intercept and variance 
  # of theta are provided
  thetaparams.init.flag <- F 
  if (is.estimateIRT) {
    # When IRT parameters are unknown, fix the intercept and variance 
    # of latent trait 
    theta.interp0 <- 0
    theta.var0 <- 1
    thetaparams.init.flag <- T
  } else {
    # When IRT parameters are known and initial values of intercept and 
    # variance of theta are inputted 
    if (all(c(!is.null(args$theta.interp0), !is.null(args$theta.var0)))) {
      theta.interp0 <- args$theta.interp0
      theta.var0 <- args$theta.var0
      thetaparams.init.flag <- T
    }
  }
  
  if (is.knockoff & (!thetaparams.init.flag || is.estimateIRT)) {
    stop('Provide the values of all nuisance parameters in precence of knockoff variables.', call. = F)
  }
  
  if (thetaparams.init.flag) {
    if (lambda > 0) {
      # Ridge regression
      ridge.lm <- glmnet::glmnet(x = XX, 
                                 y = theta - theta.interp0,
                                 family = 'gaussian',
                                 alpha = 0, # ridge penalty
                                 lambda = lambda,
                                 intercept = F)
      beta0 <- coef(ridge.lm)[-1]
    }
    
    if (lambda == 0) {
      # Linear regression
      lm.model <- lm((theta - theta.interp0) ~ 0+as.matrix(XX))
      beta0 <- coef(lm.model)
    }
    
  } else {
    if (lambda > 0) {
      ridge.lm <- glmnet::glmnet(x = XX, 
                                 y = theta,
                                 family = 'gaussian',
                                 alpha = 0, # ridge penalty
                                 lambda = lambda)
      
      beta0 <- coef(ridge.lm)[-1]
      theta.interp0 <- coef(ridge.lm)[1]
      theta.var0 <- mean((theta - predict(ridge.lm, as.matrix(XX)))**2)
    }
    
    if (lambda == 0) {
      lm.model <- lm(theta ~ as.matrix(XX))
      coeff <- coef(lm.model)
      beta0 <- coeff[-1]
      theta.interp0 <- coeff[1]
      theta.var0 <- mean(lm.model$residuals ** 2)
    }
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
  a.vec.avg <- 0
  d1.vec.avg <- 0
  d2.vec.avg <- 0
  
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
    
    ## Sampling theta
    theta <- theta.sample(XX=XX, Y=Y, beta=beta, theta.interp=theta.interp, theta.var=theta.var,
                          a.vec=a.vec, d1.vec=d1.vec, d2.vec=d2.vec)
    res.sample <- XX.sample(XX=XX, XX.ind=XX.ind, X=X, X.under=X.under, theta=theta,
                            beta=beta, theta.interp=theta.interp, theta.var=theta.var, Omega=Omega,
                            bin.params=bin.params, ord.params=ord.params, con.params=con.params,
                            bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind,
                            thresh.bin=thresh.bin, thresh.ord=thresh.ord, is.knockoff=is.knockoff,
                            X.knock.under=X.knock.under, W.list=W.list, W.list.knock=W.list.knock,
                            mea.ind=mea.ind, mea.params=mea.params, miss.mask=miss.mask,
                            meantrans.mat=meantrans.mat, ord.max = ord.max)


    X.under <- res.sample$X.under
    XX <- res.sample$XX
    if (is.knockoff) {
      X.knock.under <- res.sample$X.knock.under
    }
    
    
    if((!is.estimateIRT) & (!is.knockoff)) {
      # When IRT parameters are known and there exists knockoff variables, 
      # update intercept and variance of theta.
      
      if (lambda > 0) {
        ridge.lm <- glmnet::glmnet(x = as.matrix(XX),
                                   y = theta,
                                   family = 'gaussian',
                                   alpha = 0, # ridge penalty
                                   lambda = lambda)
        
        beta <- coef(ridge.lm)[-1]
        theta.interp <- coef(ridge.lm)[1]
        theta.var <- mean((theta - predict(ridge.lm, as.matrix(XX)))**2)
      }
      
      if (lambda == 0) {
        lm.model <- lm(theta ~ as.matrix(XX))
        coeffs <- coef(lm.model)
        beta <- coeffs[-1]
        theta.interp <- coeffs[1]
        theta.var <- mean(lm.model$residuals**2)
      }
      
    } else {
      # When IRT parameters are unknown or there exists knockoff variables, 
      # do not update intercept and variance
      
      if (lambda > 0) {
        ridge.lm <- glmnet::glmnet(x = as.matrix(XX),
                                   y = theta - theta.interp,
                                   family = 'gaussian',
                                   alpha = 0, # ridge penalty
                                   lambda = lambda,
                                   intercept = F)
        beta <- coef(ridge.lm)[-1]
      }
      
      if (lambda == 0) {
        lm.model <- lm((theta - theta.interp) ~ 0+as.matrix(XX))
        beta <- coef(lm.model)
      }
    }
    
    # Update item parameters by GLM when IRT parameters are unknown.
    if (is.estimateIRT) {
      for (j in dich.ind) {
        glm.model <- glm(Y[,j] ~ theta, family = 'binomial')
        coeffs <- coef(glm.model)
        d1.vec[j] <- coeffs[1]
        a.vec[j] <- coeffs[2]
      }
      
      #TODO: add update formula for polytomous items
    }
    
    if (is.beta.trace) {
      beta.record <- cbind(beta.record, beta)
    }
    
    if(iter > burn_in) {
      beta.avg <- ((iter-burn_in-1)*beta.avg + beta)/(iter-burn_in)
      theta.interp.avg <- ((iter-burn_in-1)*theta.interp.avg + theta.interp)/(iter-burn_in)
      theta.var.avg <- ((iter-burn_in-1)*theta.var.avg + theta.var)/(iter-burn_in)
      a.vec.avg <- ((iter-burn_in-1)*a.vec.avg + a.vec)/(iter-burn_in)
      d1.vec.avg <- ((iter-burn_in-1)*d1.vec.avg + d1.vec)/(iter-burn_in)
      d2.vec.avg <- ((iter-burn_in-1)*d2.vec.avg + d2.vec)/(iter-burn_in)
    }
    
    # TODO: Add stopping criterion
  }
  
  # Save results
  result <- list()
  result$iter <- iter
  beta.avg <- setNames(beta.avg, XX.names)
  result$beta <- beta.avg
  result$interp <- theta.interp.avg
  result$var <- theta.var.avg
  result$a.vec <- a.vec.avg
  result$d1.vec <- d1.vec.avg
  result$d2.vec <- d2.vec.avg
  result$XX <- XX
  result$XX.ind <- XX.ind
  result$X.under <- X.under
  result$X.knock.under <- X.knock.under
  if (is.beta.trace) {
    result$beta.record <- beta.record
  }
  
  return(result)
}