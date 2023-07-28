# Compute row max
rowmax <- function(x) {
  apply(x, 1, max, na.rm = T)
}

# Use the estimated parameters to obtain scales of the standarized covariates
params.to.scale <- function(bin.ind=NULL, ord.ind=NULL, con.ind=NULL,
                            bin.params=NULL, ord.params=NULL, con.params=NULL,
                            mea.ind=NULL) {
  
  p <- length(bin.ind) + length(ord.ind) + length(con.ind)
  if ((p + length(mea.ind)) == 0) {
    stop('The input indices are null!', call. = F)
  }
  
  est.scale <- list()
  for (j in 1:p) {
    # The scale of binary variables
    if (j %in% bin.ind) {
      loc <- which(bin.ind == j)
      prob <- pnorm(bin.params[loc])
      est.scale[[j]]<- sqrt(prob * (1-prob))
    }
    
    # The scale of ordinal variables
    if (j %in% ord.ind) {
      loc <- which(ord.ind == j)
      prob <- pnorm(ord.params[loc,])
      prob <- prob[which(prob < 1)]
      
      temp1 <- outer(rep(1, length(prob)), prob)
      temp1[upper.tri(temp1)] <- 0
      temp2 <- outer(1-prob, rep(1, length(prob)))
      temp2[upper.tri(temp2)] <- 0
      
      cov.mat <- temp1 * temp2
      cov.mat <- t(cov.mat)
      cov.mat[lower.tri(cov.mat)] <- (temp1 * temp2)[lower.tri(cov.mat)]
      
      cov.mat.eigen <- eigen(cov.mat)
      D <- diag(cov.mat.eigen$values)
      V <- cov.mat.eigen$vectors
      
      est.scale[[j]] <- V %*% sqrt(D) %*% t(V)
    }
    
    # The scale of continuous variables
    if (j %in% con.ind) {
      loc <- which(con.ind == j)
      est.scale[[j]] <- con.params$sd[loc]
    }
  }
  
  est.scale <- c(est.scale, rep(1, length(mea.ind)))
  return(est.scale)
}

knockoff.stat <- function(beta, est.scale, bin.ind=NULL, ord.ind=NULL, con.ind=NULL, 
                          mea.ind=NULL) {
  pp <- length(beta)
  beta.true <- beta[1:(pp/2)]
  beta.knock <- beta[-(1:(pp/2))]
  
  if (length(ord.ind) == 0) {
    est.scale <- unlist(est.scale)
    w.vec <- sign(abs(beta.true) - abs(beta.knock)) * rowmax(cbind(abs(beta.true), abs(beta.knock)))
    w.vec <- w.vec * est.scale
  }
  
  if (length(ord.ind) > 0) {
    if ((length(bin.ind) + length(con.ind) + length(mea.ind)) == 0) {
      stop('There exist ordinal variables. Indices are needed.', call. = F)
    }
    
    p <- length(bin.ind) + length(con.ind) + length(ord.ind) + length(mea.ind)
    w.vec <- c()
    len.temp <- 0
    for (j in 1:p) {
      if (!(j %in% ord.ind)) {
        ind.temp <- len.temp + 1
        stat.true.j <- abs(est.scale[[j]] * beta.true[ind.temp])
        stat.knock.j <- abs(est.scale[[j]] * beta.knock[ind.temp])
        w.vec[j] <- sign(stat.true.j - stat.knock.j) * max(stat.true.j, stat.knock.j)
        len.temp <- len.temp + 1
      }
      
      if (j %in% ord.ind) {
        ind.temp <- len.temp + 1:dim(est.scale[[j]])[1]
        stat.true.j <- norm(est.scale[[j]] %*% beta.true[ind.temp], type="2")/sqrt(length(ind.temp))
        stat.knock.j <- norm(est.scale[[j]] %*% beta.knock[ind.temp], type="2")/sqrt(length(ind.temp))
        w.vec[j] <- sign(stat.true.j - stat.knock.j) * max(stat.true.j, stat.knock.j)
        len.temp <- len.temp + length(ind.temp)
      }
    }    
  }
  
  return(w.vec)
}

# Knockoff selection
# knockoff.select <- function(w.vec, v=1) {
#   w.sign <- sign(w.vec)
#   # Controlling PFER by selecting proper threshold
#   tau <- 0
#   for(t in sort(abs(w.vec[w.vec != 0]), decreasing = T)){
#     if(sum(w.vec <= -t) == v){
#       tau <- t
#       break()
#     }
#   }
#   
#   # The selection indicator vector
#   select.ind <- (w.vec >= tau) * 1
#   return(select.ind)
# }

# Knockoff selection for controlling PFER
knockoff.select <- function(w.vec, v=1) {
  w.sign <- sign(w.vec)
  # Controlling PFER by selecting proper threshold
  tau <- -Inf
  for(t in sort(abs(w.vec[w.vec != 0]), decreasing = F)){
    if(1 + sum(w.vec < -t) == v){
      tau <- t
      break()
    }
  }
  # The selection indicator vector
  select.ind <- (w.vec > tau) * 1
  return(select.ind)
}

knockoff.select.FDR <- function(w.vec, q=0.1) {
  w.sign <- sign(w.vec)
  # Controlling PFER by selecting proper threshold
  tau <- 0
  for(t in sort(abs(w.vec[w.vec != 0]), decreasing = F)){
    if(((1 + sum(w.vec <= -t))/sum(w.vec >= t)) <= q){
      tau <- t
      break()
    }
  }
  
  # The selection indicator vector
  select.ind <- (w.vec >= tau) * 1
  return(select.ind)
}