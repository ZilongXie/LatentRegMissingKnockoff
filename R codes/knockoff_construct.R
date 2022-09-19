# Check if positive definite
is.posdef <- function(X){min(eigen(x = X, symmetric = T, only.values = T)$value) > 0}

# # Full Gibbs sampling for underlying variables
# Gibbs.sample <- function(X, Omega, bin.ind=NULL, ord.ind=NULL, con.ind=NULL,
#                          bin.params=NULL, ord.params=NULL, con.params=NULL, 
#                          max_iter=10000, is.print.iter=F, ...) {
#   N <- dim(X)[1]
#   p <- dim(X)[2]
#   miss.mask <- is.na(X)
#   args <- list(...)
#   
#   if (length(bin.ind) + length(ord.ind) + length(con.ind) != dim(X)[2]) {
#     stop('Dimensions of indices are wrong! Check your input indices.')
#   }
# 
#   if(!is.null(bin.ind)) {
#     # Thresholds of binary variables
#     thresh.bin <- list(lower = matrix(-Inf, nrow = N, ncol = length(bin.ind)),
#                        upper = matrix(Inf, nrow = N, ncol = length(bin.ind)))
#     for (j in 1:length(bin.ind)) {
#       bj <- bin.ind[j]
#       thresh.bin$lower[X[,bj] == 1, j] <- bin.params[j]
#       thresh.bin$upper[X[,bj] == 0, j] <- bin.params[j]
#     }
#   }
# 
#   if(!is.null(ord.ind)) {
#     ord.max <- apply(X[,ord.ind], 2, max, na.rm = T)
# 
#     # Thresholds of ordinal variables
#     thresh.ord <- list(lower = matrix(-Inf, nrow = N, ncol = length(ord.ind)),
#                        upper = matrix(Inf, nrow = N, ncol = length(ord.ind)))
#     for (j in 1:length(ord.ind)) {
#       oj <- ord.ind[j]
#       for (v in 1:(ord.max[j]+1)) {
#         if (v > 1){
#           thresh.ord$lower[which(X[,oj] == v), j] <- ord.params[j, v-1]
#         }
#         if (v < (ord.max[j]+1)) {
#           thresh.ord$upper[which(X[,oj] == v), j] <- ord.params[j, v]
#         }
#       }
#     }
#   }
# 
#   # Allowed to provide initial underlying variables
#   if (is.null(args$X.under0)) {
#     X.under <- X.under.init(X, bin.ind, ord.ind, con.ind, thresh.bin, thresh.ord, con.params)
#   } else {
#     X.under <- args$X.under0
#   }
#   
#   # The transformation matrix that transforms X.under_{-i} to the post-mean
#   meantrans.matr <- matrix(nrow = p-1 , ncol = p)
#   for (i in 1:p) {
#     meantrans.matr[, i] <- -Omega[i, -i]/Omega[i, i]
#   }
# 
#   # Random scan Gibbs sampling
#   for (iter in 1:max_iter) {
#     # Random scan Gibbs sampling
#     random.scan <- sample(p, p)
#     for (j in random.scan) {
#       # Conditional mean and standard deviation
#       cond.sd <- sqrt(1/Omega[j, j])
#       cond.mean <- X.under[,-j] %*% meantrans.matr[, j]
# 
#       if (j %in% con.ind) {
#         miss.mask.j <- miss.mask[,j]
#         X.under[miss.mask.j, j] <- rnorm(mean = cond.mean[miss.mask.j],
#                                          sd = cond.sd,
#                                          n = sum(miss.mask.j))
#       }
#       if (j %in% bin.ind) {
#         loc <- which(bin.ind == j)
#         X.under[, j] <- truncnorm::rtruncnorm(a = thresh.bin$lower[,loc],
#                                               b = thresh.bin$upper[,loc],
#                                               mean = cond.mean,
#                                               sd = cond.sd,
#                                               n = 1)
#       }
#       if (j %in% ord.ind) {
#         loc <- which(ord.ind == j)
#         X.under[, j] <- truncnorm::rtruncnorm(a = thresh.ord$lower[,loc],
#                                               b = thresh.ord$upper[,loc],
#                                               mean = cond.mean,
#                                               sd = cond.sd,
#                                               n = 1)
#       }
#     }
#     
#     if (is.print.iter) {
#       print(iter)
#     }
#     
#     #TODO: stopping criterion
#   }
# 
#   return(X.under)
# }

# Imputed mixed data
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
                                bin.params=NULL, ord.params=NULL, con.params=NULL) {
  
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
  
  # Generate knockoff copy of true covariates
  X.knocks.impute <- Impute(X.knocks.under * NA, X.knocks.under, bin.ind, ord.ind, con.ind, 
                            bin.params, ord.params, con.params)
  
  return(list(X.knocks.under = X.knocks.under,
              X.knocks.impute = X.knocks.impute))
}
