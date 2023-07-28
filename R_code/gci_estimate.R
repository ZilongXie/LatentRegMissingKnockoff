# Truncation function
trunc <- function(x, a, b) {
  pmax(pmin(x, b), a)
}

# Compute initial values of parameters for binary variables.
bin.params.init <- function(X, bin.ind) {
  return(qnorm(1 - colMeans(X[,bin.ind], na.rm = T)))
}

# Compute initial values of parameters for ordinal variables.
ord.params.init <- function(X, ord.max, ord.ind) {
  miss.mask <- is.na(X)
  ord.params <- list()
  # The threshold parameters
  ord.params$par <- matrix(Inf, nrow = length(ord.ind), ncol = max(ord.max))
  # The reparameterized parameters
  ord.params$repar <- matrix(Inf, nrow = length(ord.ind), ncol = max(ord.max))
  
  for (j in 1:length(ord.ind)) {
    for (k in 1:ord.max[j]) {
      ord.params$par[j, k] <- qnorm(mean((X[,ord.ind[j]] <= k-1), na.rm = T))
    }
    temp <- ord.params$par[j, 1:ord.max[j]] 
    temp <- c(ord.params$par[j, 1], log(diff(temp)))
    ord.params$repar[j, 1:ord.max[j]] <- temp
  }
  return(ord.params)
}

# Compute initial values of parameters for continuous variables.
con.params.init <- function(X, con.ind) {
  con.params <- list()
  # Marginal mean
  con.params$mean <- colMeans(X[,con.ind], na.rm =  T)
  # Marginal standard deviation
  con.params$sd <- apply(X[,con.ind], 2, sd, na.rm = T)
  return(con.params)
}

# Compute initial values of parameters in measurement model for Z.
# We only consider continuous indicators under normality assumption here.
mea.init <- function(W.list, ...) {
  # Number of unobserved variables
  N <- nrow(W.list[[1]])
  num.mea <- length(W.list)
  args <- list(...)
  mea.init <- list()
  
  if (!is.null(args$mea.params)) {
    mea.init$params <- args$mea.params
    
    # Matrix of scores
    mea.scores <- matrix(0, nrow = N, ncol = num.mea)
    for (j in 1:num.mea) {
      W <- W.list[[j]]
      W.loadings <- args$mea.params$loadings[[j]]
      W.errvar <- args$mea.params$errvar[[j]]
      W.mean <- args$mea.params$mean[[j]]
      
      post.sd <- 1/sqrt(sum(W.loadings^2/W.errvar))
      post.mean <- (rowSums((W - outer(rep(1, N), W.mean))
                            * outer(rep(1, N), W.loadings/W.errvar))
                    ) * post.sd^2
      
      mea.scores[,j] <- rnorm(mean = post.mean, 
                              sd = post.sd,
                              n = N)
      
    }
    mea.init$scores <- mea.scores
    
  } else {
    # List of parameters
    mea.params <- list()
    mea.params$mean <- list()
    mea.params$loadings <- list()
    mea.params$errvar <- list()
    
    # Matrix of scores
    mea.scores <- matrix(0, nrow = N, ncol = num.mea)
    
    # Require package `psych`
    if (!requireNamespace('psych', quietly=T)) {
      stop('psych is not installed', call.=F)
    }
    
    for (j in 1:num.mea) {
      W <- W.list[[j]]
      # Perform factor analysis
      fa <- psych::fac(W, nfactors = 1, rotate = 'none', fm = 'ml')
      W.sd <- sqrt(diag(cov(W)))
      mea.params$mean[[j]] <- colMeans(W)
      mea.params$loadings[[j]] <- as.vector(fa$loadings) * W.sd
      mea.params$errvar[[j]]<- as.vector(fa$uniquenesses * W.sd**2)
      mea.scores[,j] <- as.vector(fa$scores)
    }
    mea.init$params <- mea.params
    mea.init$scores <- mea.scores
  }
  
  return(mea.init)
}

# Compute initial values of correlation matrix.
corr.init <- function(X, bin.ind, ord.ind, con.ind, X.mea=NULL) {
  # Require package `psych`
  if (!requireNamespace('psych', quietly=T)) {
    stop('psych is not installed', call.=F)
  }
  
  mea.ind <- NULL
  if (is.null(X.mea)) {
    X.data.frame <- as.data.frame(X)
  } else {
    num.mea <- ncol(X.mea)
    mea.ind <- length(bin.ind) + length(ord.ind) + length(con.ind) + 1:num.mea
    X.data.frame <- as.data.frame(cbind(X, X.mea))
  }
  
  fit <- psych::mixedCor(data = X.data.frame,
                         c = c(con.ind, mea.ind),
                         d = bin.ind,
                         p = ord.ind,
                         smooth = T,
                         correct = T,
                         global = FALSE)
  
  return(fit$rho)
}

# Trace of a square matrix
tr <- function(M) {
  if (is.null(dim(M))) {
    stop('The input is not a matrix!', call. = F)
  }
  if (dim(M)[1] != dim(M)[2]) {
    stop('The input is not a square matrix!', call. = F)
  }
  sum(diag(M), na.rm = T)
}

# Loss of given lower triangular matrix
tri.loss <- function(B, R) {
  S.inv <- solve(B %*% t(B))
  return(tr(R %*% S.inv) - log(sum(diag(B)^2)))
}

# Gradient of given lower triangular matrix
tri.grad <- function(B, R) {
  S.inv <- solve(B %*% t(B))
  grad.B <- 2 * (-S.inv %*% R %*% S.inv %*% B - diag(1/diag(B)))
  grad.B[upper.tri(grad.B)] <- 0
  return(grad.B)
}

# Determine whether the matrix is positive definite, if not, project it onto to 
# the cone of semidefinite matrix by a proximal gradient descent algorithm 
correct.pd <- function(M, tol=1e-4, max.iter=1000) {
  if (!all.equal(M, t(M))) {
    stop('The input matrix is not symmetric!', call. = F)
  }
  
  if (min(eigen(M)$value)> 0) {
    return(M)
  } else {
    M.eigen <- eigen(M, symmetric = T)
    U <- M.eigen$vectors
    D <- diag(M.eigen$values)
    D[D<=0] <- 0
    R <- U %*% D %*% t(U)
    
    B <- diag(dim(R)[1])
    for (iter in 1:max.iter) {
      loss.old <- tri.loss(B, R)
      loss.new <- Inf
      alpha <- 2
      # Search the stepsize
      while(loss.new > loss.old) {
        alpha <- alpha/2
        B.temp <- B - alpha * tri.grad(B, R)
        B.temp <- B.temp/sqrt(rowSums(B.temp^2))
        loss.new <- tri.loss(B.temp, R)
      }
      B <- B.temp
      if (abs(loss.new - loss.old) < tol) {
        break
      }
    }
    
    return(B %*% t(B))
  }
}

# Initialize underlying variables
X.under.init <- function(X, bin.ind=NULL, ord.ind=NULL, con.ind=NULL, 
                         thresh.bin=NULL, thresh.ord=NULL, con.params=NULL,
                         X.mea = NULL) {
  # Initialization
  if (is.null(X.mea)) {
    X.under <- matrix(nrow = dim(X)[1], ncol = dim(X)[2])
  } else {
    num.mea <- ncol(X.mea)
    mea.ind <- length(bin.ind) + length(ord.ind) + length(con.ind) + 1:num.mea
    X.under <- cbind(matrix(nrow = dim(X)[1], ncol = dim(X)[2]), X.mea)
  }
  
  # Require package `truncnorm`
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

# One-round Gibbs sampling for underlying variables
X.under.sample <- function(X, X.under, bin.ind=NULL, ord.ind=NULL, con.ind=NULL, 
                           thresh.bin=NULL, thresh.ord=NULL, con.params=NULL, 
                           mat.params=NULL, W.list=NULL, mea.ind=NULL, mea.params=NULL) {
  
  if (all(!is.null(W.list), !is.null(mea.ind), !is.null(mea.params))) {
    num.mea <- length(mea.ind)
    p <- ncol(X) + num.mea
  } else {
    p <- ncol(X)
  }
  N <- nrow(X)
  miss.mask <- is.na(X)
  
  # Check dimensions
  stopifnot(ncol(X.under) == p)
  stopifnot(all(dim(mat.params$Sigma)[1] == p, dim(mat.params$Omega)[1] == p))
  
  # The transformation matrix that transforms X.under_{-i} to the post-mean
  meantrans.matr <- matrix(nrow = p-1 , ncol = p)
  for (i in 1:p) {
    meantrans.matr[, i] <- -mat.params$Omega[i, -i]/mat.params$Omega[i, i]
  }
  
  # Transform continuous variables
  if (length(con.ind) > 0) {
    X.under.con <- X.under[,con.ind]
    X.under.con.obs <- t((t(X[, con.ind]) - con.params$mean)/con.params$sd)
    X.under.con[!is.na(X[, con.ind])] <- X.under.con.obs[!is.na(X[, con.ind])]
    X.under[,con.ind] <- X.under.con
  }
  
  # Random scan Gibbs sampling
  random.scan <- sample(p, p)
  
  # Gibbs sampler
  for (j in random.scan) {
    # Conditional mean and standard deviation
    cond.sd <- sqrt(1/mat.params$Omega[j, j])
    cond.mean <- X.under[,-j] %*% meantrans.matr[, j]
    
    if (j %in% con.ind) {
      miss.mask.j <- miss.mask[,j]
      X.under[miss.mask.j, j] <- rnorm(mean = cond.mean[miss.mask.j], 
                                       sd = cond.sd,
                                       n = sum(miss.mask.j))
    }
    
    if (j %in% bin.ind) {
      loc <- which(bin.ind == j)
      X.under[, j] <- truncnorm::rtruncnorm(a = thresh.bin$lower[,loc],
                                            b = thresh.bin$upper[,loc],
                                            mean = cond.mean,
                                            sd = cond.sd,
                                            n = 1)
    }
    
    if (j %in% ord.ind) {
      loc <- which(ord.ind == j)
      X.under[, j] <- truncnorm::rtruncnorm(a = thresh.ord$lower[,loc],
                                            b = thresh.ord$upper[,loc],
                                            mean = cond.mean,
                                            sd = cond.sd,
                                            n = 1)
    }
    
    if (j %in% mea.ind) {
      loc <- which(mea.ind == j)
      W <- W.list[[loc]]
      W.loadings <- mea.params$loadings[[loc]]
      W.errvar <- mea.params$errvar[[loc]]
      W.mean <- mea.params$mean[[loc]]
      
      post.sd <- 1/sqrt(sum(W.loadings**2/W.errvar) + 1/cond.sd**2)
      
      post.mean <- (rowSums((W - outer(rep(1, N), W.mean))
                            * outer(rep(1, N), W.loadings/W.errvar))
                    + cond.mean/cond.sd**2) * post.sd**2

      X.under[, j] <- rnorm(mean = post.mean, 
                            sd = post.sd,
                            n = N)
    
    }
  }
  
  return(X.under)
}

# Compute first and second order derivatives of continuous marginal parameters
con.params.derivs <- function(X, X.under, con.ind, con.params, mat.params, is.second = T) {
  stopifnot(length(con.ind) > 0)
  
  # Transform continuous variables
  X.under.con <- X.under[,con.ind]
  X.under.con.obs <- t((t(X[, con.ind]) - con.params$mean)/con.params$sd)
  X.under.con[!is.na(X[, con.ind])] <- X.under.con.obs[!is.na(X[, con.ind])]
  X.under[,con.ind] <- X.under.con
  
  # Quantities to be used
  lognorm.grad.con.obs <- (X.under %*% mat.params$Omega)[,con.ind]
  lognorm.grad.con.obs[is.na(X[,con.ind])] <- NA
  Omega.con <- mat.params$Omega[cbind(con.ind, con.ind)]
  X.under.con.obs <- X.under[,con.ind]
  X.under.con.obs[is.na(X[,con.ind])] <- NA
  
  # First-order derivatives
  derivs <- list()
  derivs$first$mean <- colSums(lognorm.grad.con.obs, na.rm = T)/con.params$sd
  derivs$first$sd <- colSums(lognorm.grad.con.obs * X.under.con.obs - 1, na.rm = T)/con.params$sd
  
  # Second-order derivatives
  if (is.second) {
    derivs$second$mean <- -colSums(!is.na(X[,con.ind])) * Omega.con/con.params$sd^2
    derivs$second$sd <- colSums(-t(t(X.under.con.obs^2) * Omega.con) 
                                - 2 * lognorm.grad.con.obs * X.under.con.obs 
                                + 1, na.rm = T)/con.params$sd^2
  }
  return(derivs)
}

# Compute first and second order derivatives of binary marginal parameters
bin.params.derivs <- function(X, X.under, bin.ind, bin.params, thresh.bin, mat.params, is.second = T) {
  stopifnot(length(bin.ind) > 0)
  
  mask <- is.na(X)
  derivs <- list()
  derivs$first <- rep(0, length(bin.ind))
  if (is.second) {
    derivs$second <- rep(0, length(bin.ind))
  }
  for (j in 1:length(bin.ind)) {
    dj <- bin.ind[j]
    miss.mask.j <- mask[,dj]
    
    # The conditional mean and standard deviation
    cond.mean <- (-X.under[,-dj] %*% mat.params$Omega[dj,-dj]/mat.params$Omega[dj,dj])
    cond.sd <- sqrt(1/mat.params$Omega[dj,dj])
    
    # The truncated normal density at boundaries
    trunc.density <- truncnorm::dtruncnorm(x = bin.params[j],
                                           a = thresh.bin$lower[!miss.mask.j, j],
                                           b = thresh.bin$upper[!miss.mask.j, j],
                                           mean = cond.mean[!miss.mask.j],
                                           sd = cond.sd)
    
    # The first derivatives
    x.obs.j <- X[!miss.mask.j, dj]
    derivs$first[j] <- sum((-1)^x.obs.j * trunc.density)
    
    # The second derivatives
    if (is.second) {
      cond.norm <- -(bin.params[j] - cond.mean[!miss.mask.j])/cond.sd^2
      derivs$second[j] <- sum((-1)^x.obs.j * cond.norm * trunc.density - trunc.density^2)
    }
  }
  return(derivs)
}

# Compute first and second order derivatives of ordinal marginal parameters, for the reparametrized form
ord.params.derivs <- function(X, X.under, ord.ind, ord.max, ord.levels.which, ord.params, thresh.ord, mat.params, is.second = T) {
  stopifnot(length(ord.ind) > 0)
  
  miss.mask <- is.na(X)
  derivs <- list()
  derivs$first <- matrix(0, nrow = length(ord.ind), ncol = max(ord.max))
  
  if (is.second) {
    # Here unused entries are set to 1
    derivs$second <- matrix(1, nrow = length(ord.ind), ncol = max(ord.max)) 
  }
  
  for (j in 1:length(ord.ind)) {
    oj <- ord.ind[j]
    cond.mean <- (-X.under[,-oj] %*% mat.params$Omega[oj,-oj]/mat.params$Omega[oj,oj])
    cond.sd <- sqrt(1/mat.params$Omega[oj,oj])
    lower.trunc.density <- truncnorm::dtruncnorm(x = thresh.ord$lower[,j],
                                                 a = thresh.ord$lower[,j],
                                                 b = thresh.ord$upper[,j],
                                                 mean = cond.mean,
                                                 sd = cond.sd)
    upper.trunc.density <- truncnorm::dtruncnorm(x = thresh.ord$upper[,j],
                                                 a = thresh.ord$lower[,j],
                                                 b = thresh.ord$upper[,j],
                                                 mean = cond.mean,
                                                 sd = cond.sd)
    par.upper.first <- c()
    par.lower.first <- c()
    for (v in 1:(ord.max[j]+1)) {
      if (v < ord.max[j] + 1) {
        par.upper.first[v] <- sum(upper.trunc.density[ord.levels.which[[j]][[v]]])
      }
      if (v > 1) {
        par.lower.first[v-1] <- -sum(lower.trunc.density[ord.levels.which[[j]][[v]]])
      }
    }
    # The first derivatives of original parameters
    par.first <- par.upper.first + par.lower.first
    
    # Chain rule for reparametrized parameters
    par.to.repar <- ord.params$repar[j,1:ord.max[j]]
    par.to.repar[1] <- 1
    par.to.repar[-1] <- exp(par.to.repar[-1])
    par.to.repar.mat1 <- outer(par.to.repar, rep(1, ord.max[j]))
    par.to.repar.mat1[lower.tri(par.to.repar.mat1)] <- 0
    derivs$first[j,1:ord.max[j]] <- par.to.repar.mat1 %*% par.first
    
    # Second derivatives
    if (is.second) {
      par.upper.hess <- matrix(0, nrow = ord.max[j], ncol = ord.max[j]) 
      par.lower.hess <- matrix(0, nrow = ord.max[j], ncol = ord.max[j])
      
      for (v in 1:(ord.max[j]+1)) {
        if (v < ord.max[j] + 1) {
          cond.norm <- -(ord.params$par[j,v] - cond.mean[ord.levels.which[[j]][[v]]])/cond.sd^2
          par.upper.hess[v, v] <- sum((cond.norm - upper.trunc.density[ord.levels.which[[j]][[v]]]) 
                                      * upper.trunc.density[ord.levels.which[[j]][[v]]])
          if (v > 1) {
            par.upper.hess[v, v-1] <- sum(upper.trunc.density[ord.levels.which[[j]][[v]]]
                                          * lower.trunc.density[ord.levels.which[[j]][[v]]])
          }
        }
        if (v > 1) {
          cond.norm <- -(ord.params$par[j,v-1] - cond.mean[ord.levels.which[[j]][[v]]])/cond.sd^2
          par.lower.hess[v-1, v-1] <- sum(-(cond.norm + lower.trunc.density[ord.levels.which[[j]][[v]]]) 
                                          * lower.trunc.density[ord.levels.which[[j]][[v]]])
          if (v < ord.max[j] + 1) {
            par.lower.hess[v-1, v] <- sum(upper.trunc.density[ord.levels.which[[j]][[v]]]
                                          * lower.trunc.density[ord.levels.which[[j]][[v]]])
          }
        }
      }
      # Hessian matrix of the original parameters
      par.hess <- par.upper.hess + par.lower.hess
      
      # Chain rule for reparametrized parameters
      par.to.repar.mat2 <- par.to.repar.mat1
      par.to.repar.mat2[1,] <- 0
      derivs$second[j,1:ord.max[j]] <- (diag(par.to.repar.mat1 %*% par.hess %*% t(par.to.repar.mat1)) + par.to.repar.mat2 %*% par.first)
    }
  }
  return(derivs)
}

# Compute first and second order derivatives of the Cholesky decomposition matrix
mat.params.derivs <- function(X.under, mat.params, is.second = T) {
  N <- dim(X.under)[1]
  derivs <- list()
  inv.chol <- solve(mat.params$Sigma.chol)
  temp <- mat.params$Omega %*% t(X.under) %*% X.under %*% t(inv.chol)
  derivs$first <- (temp - N * diag(diag(inv.chol)))
  derivs$first[upper.tri(derivs$first)] <- 0
  
  if (is.second) {
    derivs$second <- (-2 * t(inv.chol) * temp
                      - outer(diag(mat.params$Omega), rep(1, N)) %*% (X.under %*% t(inv.chol))^2
                      + N * diag(diag(inv.chol)^2))
  }
  
  return(derivs)
}

# Estimation of Gaussian copula model
GCI.estimate <- function(X, bin.ind=NULL, ord.ind=NULL, con.ind=NULL,
                         eps=1e-2, max.iter=500, burn.in=300, print.iter=F,
                         print.max.change=F, is.params.trace=F,
                         is.return.under=F, W.list=NULL, ...) {

  if (length(bin.ind) + length(ord.ind) + length(con.ind) != ncol(X)) {
    stop('Dimensions of indices are wrong! Check your input indices.')
  }

  args <- list(...)
  N <- nrow(X)
  p <- ncol(X)

  mea.ind <- NULL
  if (!is.null(W.list)) {
    num.mea <- length(W.list)
    mea.ind <- length(bin.ind) + length(ord.ind) + length(con.ind) + 1:num.mea
    p <- ncol(X) + num.mea
  }

  miss.mask <- is.na(X)
  col.num_iss <- colSums(miss.mask)
  params.avg <- list()

  # Record the trace of parameters if required
  if (is.params.trace) {
    params.trace <- list()
  }

  # Compute initial values of marginal parameters and thresholds
  # Parameters in measurement models
  mea.params <- NULL
  X.mea <- NULL
  if (length(mea.ind) > 0) {
    if (!is.null(args$mea.params)) {
      mea.params <- args$mea.params
      
      stopifnot(all(length(mea.params$mean) == length(mea.ind),
                    length(mea.params$loadings) == length(mea.ind),
                    length(mea.params$errvar) == length(mea.ind)))
    }
    
    mea.init.res <- mea.init(W.list=W.list, mea.params=mea.params)
    mea.params <- mea.init.res$params
    X.mea <- mea.init.res$scores
    params.avg$mea.params <- mea.params
  }

  # Continuous parameters
  con.params <- NULL
  if (length(con.ind) > 0) {
    if (is.null(args$con.params0)) {
      # Initialize parameters if initial values are not provided
      con.params <- con.params.init(X=X, con.ind=con.ind)
    } else {
      con.params <- args$con.params0
    }

    # Record trace of parameters
    if (is.params.trace) {
      params.trace$con.params.mean <- cbind(params.trace$con.params.mean, con.params$mean)
      params.trace$con.params.sd <- cbind(params.trace$con.params.sd, con.params$sd)
    }

    # Initialize averaged continuous parameters
    params.avg$con.params.mean <- 0
    params.avg$con.params.sd <- 0
  }

  # Binary parameters
  thresh.bin <- NULL
  bin.params <- NULL
  if (length(bin.ind) > 0) {
    if (!all(apply(X[,bin.ind], 2, min, na.rm = T) == 0) ||
        !all(apply(X[,bin.ind], 2, max, na.rm = T) == 1)) {
      stop('Binary variables should be 0-1 coded!', call. = F)
    }

    if (is.null(args$bin.params0)) {
      # Initialize parameters if initial values are not provided
      bin.params <- bin.params.init(X=X, bin.ind=bin.ind)
    } else {
      bin.params <- args$bin.params0
    }

    # Record which samples fall in the same levels
    bin.levels.which <- list()
    for (j in 1:length(bin.ind)) {
      bj <- bin.ind[j]
      bin.levels.which.j <- list()
      for (v in 0:1) {
        bin.levels.which.j[[v+1]] <- which(X[,bj] == v)
      }
      bin.levels.which[[j]] <- bin.levels.which.j
    }

    # Thresholds of binary variables
    thresh.bin <- list(lower = matrix(-Inf, nrow = N, ncol = length(bin.ind)),
                       upper = matrix(Inf, nrow = N, ncol = length(bin.ind)))
    for (j in 1:length(bin.ind)) {
      thresh.bin$lower[bin.levels.which[[j]][[2]], j] <- bin.params[j]
      thresh.bin$upper[bin.levels.which[[j]][[1]], j] <- bin.params[j]
    }

    # Record trace of parameters
    if (is.params.trace) {
      params.trace$bin.params <- cbind(params.trace$bin.params, bin.params)
    }

    # Initialize averaged binary parameters
    params.avg$bin.params <- 0
  }

  # Ordinal parameters
  thresh.ord <- NULL
  ord.params <- NULL
  if (length(ord.ind) > 0) {
    if (!all(apply(X[,ord.ind], 2, min, na.rm = T) == 0)) {
      stop('Ordinal variables should start from 0!', call. = F)
    }

    ord.max <- apply(X[,ord.ind], 2, max, na.rm = T)

    if (is.null(args$ord.params0)) {
      # Initialize parameters if initial values are not provided
      ord.params <- ord.params.init(X=X, ord.max=ord.max, ord.ind=ord.ind)
    } else {
      if(is.null(args$ord.params0$par) || is.null(args$ord.params0$repar)) {
        stop('Please provide initial values fro both the raw parameres and reparametrized parameters.', call. = F)
        ord.params <- ord.params.init(X=X, ord.max=ord.max, ord.ind=ord.ind)
      }
      else {
        ord.params <- args$ord.params0
      }
    }

    # Record which samples fall in the same level
    ord.levels.which <- list()
    for (j in 1:length(ord.ind)) {
      oj <- ord.ind[j]
      ord.levels.which.j <- list()
      for (v in 0:ord.max[j]) {
        ord.levels.which.j[[v+1]] <- which(X[,oj] == v)
      }
      ord.levels.which[[j]] <- ord.levels.which.j
    }

    # Thresholds of ordinal variables
    thresh.ord <- list(lower = matrix(-Inf, nrow = N, ncol = length(ord.ind)),
                       upper = matrix(Inf, nrow = N, ncol = length(ord.ind)))
    for (j in 1:length(ord.ind)) {
      for (v in 1:(ord.max[j]+1)) {
        if (v > 1) {
          thresh.ord$lower[ord.levels.which[[j]][[v]], j] <- ord.params$par[j, v-1]
        }
        if (v < (ord.max[j]+1)) {
          thresh.ord$upper[ord.levels.which[[j]][[v]], j] <- ord.params$par[j, v]
        }
      }
    }

    # Record trace of parameters
    if (is.params.trace) {
      params.trace$ord.params.vec <- cbind(params.trace$ord.params.vec, as.vector(ord.params$par))
    }

    # Initialize averaged ordinal parameters
    params.avg$ord.params <- 0
  }

  # Correlation matrix
  if (is.null(args$Sigma0)) {
    # Initialize correlation matrix if initial values are not provided
    Sigma <- correct.pd(corr.init(X=X, bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind, X.mea=X.mea))
    Sigma.chol <- t(chol(Sigma))
    Omega <- solve(Sigma)
    mat.params <- list(Sigma = Sigma, Omega = Omega, Sigma.chol = Sigma.chol)

  } else {
    Sigma <- correct.pd(args$Sigma0)
    Sigma.chol <- t(chol(Sigma))
    Omega <- solve(Sigma)
    mat.params <- list(Sigma = Sigma, Omega = Omega, Sigma.chol = Sigma.chol)
  }

  # Record trace of parameters
  if (is.params.trace) {
    params.trace$Sigma.lower <- cbind(params.trace$Sigma.lower, as.vector(mat.params$Sigma[lower.tri(mat.params$Sigma)]))
  }

  # Initialize averaged correlation matrix
  params.avg$Sigma <- 0

  # Initialize underlying variables
  X.under <- X.under.init(X=X, bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind,
                          thresh.bin=thresh.bin, thresh.ord=thresh.ord,
                          con.params=con.params, X.mea=X.mea)

  # If required, perform additional Gibbs sampling at the initialization stage
  # for better mixing
  if (!is.null(args$init.Gibbs.round)) {
    for (num in 1:args$init.Gibbs.round) {
      X.under <- X.under.sample(X=X, X.under=X.under, bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind,
                                thresh.bin=thresh.bin, thresh.ord=thresh.ord, con.params=con.params,
                                mat.params=mat.params, W.list=W.list, mea.ind=mea.ind, mea.params=mea.params)
    }
  }

  # Record maximum change
  if (print.max.change) {
    max.change <- c()
    params.old <- c(unlist(con.params), unlist(bin.params),
                    unlist(ord.params), unlist(mat.params))
  }

  # Store second order derivatives
  second.derivs <- list()
  
  # Proximal stochastic gradient algorithm
  for (t in 1:max.iter) {
    # Print iterations
    if (print.iter) {
      print(t)
    }

    # Sampling step
    X.under <- X.under.sample(X=X, X.under=X.under, bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind,
                              thresh.bin=thresh.bin, thresh.ord=thresh.ord, con.params=con.params,
                              mat.params=mat.params, W.list=W.list, mea.ind=mea.ind, mea.params=mea.params)

    # Update approximate Hessian diagonal matrix during burn-in period
    if (t <= burn.in) {
      is.second <- TRUE
    } else {
      is.second <- FALSE
    }

    # Stepsize
    step.size <- (t+100)^(-0.51)

    # Update continuous marginal parameters
    if (length(con.ind) > 0) {
      con.derivs <- con.params.derivs(X=X, X.under=X.under, con.ind=con.ind,
                                      con.params=con.params, mat.params=mat.params,
                                      is.second=is.second)

      if (is.second) {
        second.derivs$con.params.mean <- con.derivs$second$mean
        second.derivs$con.params.sd <- con.derivs$second$sd
        # TODO: truncation
      }

      # Update parameters
      con.params$mean <- con.params$mean - step.size * con.derivs$first$mean/second.derivs$con.params.mean
      con.params$sd <- con.params$sd - step.size * con.derivs$first$sd/second.derivs$con.params.sd

      # Record trace of parameters
      if (is.params.trace) {
        params.trace$con.params.mean <- cbind(params.trace$con.params.mean, con.params$mean)
        params.trace$con.params.sd <- cbind(params.trace$con.params.sd, con.params$sd)
      }
    }

    # Update binary marginal parameters
    if (length(bin.ind) > 0) {
      bin.derivs <- bin.params.derivs(X=X, X.under=X.under, bin.ind=bin.ind, bin.params=bin.params,
                                      thresh.bin=thresh.bin, mat.params=mat.params,
                                      is.second=is.second)

      if (is.second) {
        second.derivs$bin.params <- bin.derivs$second
        # TODO: truncation
      }

      # Update parameters
      bin.params <- bin.params - step.size * bin.derivs$first/second.derivs$bin.params

      # Record trace of parameters
      if (is.params.trace) {
        params.trace$bin.params <- cbind(params.trace$bin.params, bin.params)
      }

      # Update thresholds
      for (j in 1:length(bin.ind)) {
        thresh.bin$lower[bin.levels.which[[j]][[2]], j] <- bin.params[j]
        thresh.bin$upper[bin.levels.which[[j]][[1]], j] <- bin.params[j]
      }
    }

    # Update ordinal marginal parameters
    if (length(ord.ind) > 0) {
      ord.derivs <- ord.params.derivs(X=X, X.under=X.under, ord.ind=ord.ind,
                                      ord.max=ord.max, ord.levels.which=ord.levels.which,
                                      ord.params=ord.params, thresh.ord=thresh.ord,
                                      mat.params=mat.params, is.second=is.second)

      if (is.second) {
        second.derivs$ord.params <- ord.derivs$second
        # TODO: truncation
      }

      # Update parameters
      ord.params$repar <- ord.params$repar - step.size * ord.derivs$first/second.derivs$ord.params
      ord.params$par <- ord.params$repar
      ord.params$par[,-1] <- exp(ord.params$repar[,-1])
      ord.params$par <- t(apply(ord.params$par, 1, cumsum))

      # Record trace of parameters
      if (is.params.trace) {
        params.trace$ord.params.vec <- cbind(params.trace$ord.params.vec, as.vector(ord.params$par))
      }

      # Update thresholds
      for (j in 1:length(ord.ind)) {
        for (v in 1:(ord.max[j]+1)) {
          if (v > 1) {
            thresh.ord$lower[ord.levels.which[[j]][[v]], j] <- ord.params$par[j, v-1]
          }
          if (v < (ord.max[j]+1)) {
            thresh.ord$upper[ord.levels.which[[j]][[v]], j] <- ord.params$par[j, v]
          }
        }
      }
    }

    # Update correlation marginal parameters
    mat.derivs <- mat.params.derivs(X.under=X.under, mat.params=mat.params, is.second=is.second)
    if (is.second) {
      second.derivs$mat.params <- mat.derivs$second
      # TODO: truncation
    }

    mat.params$Sigma.chol <- mat.params$Sigma.chol - step.size * mat.derivs$first/second.derivs$mat.params
    row.norm <- sqrt(rowSums(mat.params$Sigma.chol^2))
    mat.params$Sigma.chol <- mat.params$Sigma.chol/row.norm
    mat.params$Sigma <- mat.params$Sigma.chol %*% t(mat.params$Sigma.chol)
    mat.params$Omega <- solve(mat.params$Sigma)

    # Record trace of parameters
    if (is.params.trace) {
      params.trace$Sigma.lower <- cbind(params.trace$Sigma.lower, as.vector(mat.params$Sigma[lower.tri(mat.params$Sigma)]))
    }

    # Monitoring max-change
    if (print.max.change) {
      temp <- c(unlist(con.params), unlist(bin.params),
                unlist(ord.params), unlist(mat.params))
      max.change <- c(max.change, max(abs(temp-params.old)))
      params.old <- temp
    }

    # Polyak-Ruppert averaging
    if (t >= burn.in) {
      t.avg <- t-burn.in
      if (length(con.ind) > 0) {
        params.avg$con.params.mean <- (t.avg * params.avg$con.params.mean + con.params$mean)/(t.avg + 1)
        params.avg$con.params.sd <- (t.avg * params.avg$con.params.sd + con.params$sd)/(t.avg + 1)
      }
      if (length(bin.ind) > 0) {
        params.avg$bin.params <- (t.avg * params.avg$bin.params + bin.params)/(t.avg + 1)
      }
      if (length(ord.ind) > 0) {
        params.avg$ord.params <- (t.avg * params.avg$ord.params + ord.params$par)/(t.avg + 1)
      }
      params.avg$Sigma <- (t.avg * params.avg$Sigma + mat.params$Sigma)/(t.avg + 1)
    }

    #TODO:Add stopping criterion
  }

  # Results
  results <- list()
  results$params.avg <- params.avg

  # Trace of parameters
  if (is.params.trace) {
    results$params.trace <- params.trace
  }

  # Return underlying variables
  if (is.return.under) {
    results$X.under <- X.under
  }

  return(results)
}