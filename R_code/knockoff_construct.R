# Check if positive definite
is.posdef <- function(X){min(eigen(x = X, symmetric = T, only.values = T)$value) > 0}

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

# MVR knockoffs
MVR.knockoffs <- function(Sigma, max_iter = 1000){
  p <- dim(Sigma)[1]
  I <- diag(p)
  # Inital value
  S <- diag((min(eigen(Sigma, symmetric = T)$value)), p)
  D <- 2 * Sigma - S
  L <- t(chol(D))
  
  for (iter in 1:max_iter){
    S.old <- S
    for (j in 1:p) {
      v_d <- forwardsolve(l = L, x = I[,j])
      v_n <- backsolve(r = t(L), x = v_d)
      c_d <- sum(v_d^2)
      c_n <- -sum(v_n^2)
      
      delta <- optimise(function(x){1/(S[j,j]+x)-x*c_n/(1-x*c_d)},
                        lower = -S[j,j],
                        upper = 1/c_d,
                        maximum = F)$minimum
      
      S[j,j] <- S[j,j] + delta
      L <- t(chol(2*Sigma - S))
    }
    if(max(abs(S - S.old)) < 1e-13){
      break()
    }
  }
  return(S)
}

# Construct knockoff copy
Knockoffs.construct <- function(X.under, Omega, S, bin.ind=NULL, ord.ind=NULL, con.ind=NULL,
                                bin.params=NULL, ord.params=NULL, con.params=NULL,
                                mea.ind=NULL, mea.params=NULL) {
  
  N <- dim(X.under)[1]
  p <- dim(X.under)[2]
  
  # The conditional mean and variance 
  # S is a diagonal matrix
  cond.cov <- 2 * S - S %*% Omega %*% S
  cond.mu <- X.under - X.under %*% Omega %*% S
  
  if (!is.posdef(cond.cov)) {
    stop('The conditional covariance matrix is not positive definite, check the knockoff diagonal!')
  }
  
  # Generate knockoff copy of underlying variables
  X.knocks.under <- cond.mu + MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = cond.cov)
  
  # Generate knockoff copy of true covariates and indicators
  W.list.knock <- NULL
  if ((length(mea.ind) > 0) & !is.null(mea.params)) {
    X.mea.knock <- X.knocks.under[,mea.ind]
    W.list.knock <- list()
    for (l in 1:length(mea.ind)) {
      W.knock <- (outer(X.mea.knock[,l], mea.params$loadings[[l]]) + 
                    outer(rep(1, N), mea.params$mean[[l]]) +
                    matrix(rnorm(N * length(mea.params$loadings[[l]])), nrow = N) *
                    outer(rep(1, N), sqrt(mea.params$errvar[[l]]))
                  )
      W.list.knock[[l]] <- W.knock
    }
    
    X.knocks.impute <- Impute(X=X.knocks.under[,-mea.ind] * NA, X.under=X.knocks.under[,-mea.ind],
                              bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind, bin.params=bin.params, 
                              ord.params=ord.params, con.params=con.params)
  } else {
    X.knocks.impute <- Impute(X=X.knocks.under * NA, X.under=X.knocks.under,
                              bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind, bin.params=bin.params, 
                              ord.params=ord.params, con.params=con.params)
  }
  
  return(list(X.knocks.under = X.knocks.under,
              X.knocks.impute = X.knocks.impute,
              W.list.knock = W.list.knock))
}
