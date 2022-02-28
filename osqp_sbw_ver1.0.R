library(Matrix)
library(MASS)
library(matrixStats)
library(osqp)

# The following function implements the OSQP-based SBW (7) 
# using up to the K-th moment and (K-1)-th order interactions 
# as our basis functions.

osqp_sbw <- function(X,A,Y,
                     delta.v=0.005, K=2, interactions=FALSE) {
  
  #################################################
  # X: covariates (n x d matrix)
  # A: treatment (n x 1 vector)
  # Y: outcome (n x 1 vector)
  # delta.v: tolerance level (a single value of a numeric vector);
  #          the numeric vector is used for multiple iterations, each with different tolerance level
  # K: use up to K-th moment
  # interactions: add the interaction terms if TRUE
  #################################################
  
  X_ <- NULL
  res.list <- list()
  for (k in 1:K){
    # add the k-th moment
    X.k <- X^k
    colnames(X.k) <- paste(colnames(X),".",k,sep = "")
    X_ <- cbind(X_, X.k)
  }
  
  if (interactions==TRUE & K >=2) {
    for (k in 2:K){
      # add the kth order interactions
      indx <- combn(colnames(X),k)
      int <- as.data.frame(do.call(cbind,
                                   lapply(split(indx, col(indx)), 
                                          function(x) rowProds(as.matrix(X[,x])))
      ))
      colnames(int) <- apply(indx, 2, function(x) paste(x,collapse="."))
      X_ <- cbind(X_, int)
    }
  }
  
  nX <- ncol(X_)
  n1 <- sum(A); n0 <- sum(1-A)
  Xt <- X_[A==1,]; Xc <- X_[A==0,]
  Yt <- Y[A==1]; Yc <- Y[A==0]
  P.mat <- as(.symDiagonal(n=n0, x=1.), "dgCMatrix")
  q.vec <- rep(-1./n0,n0)
  A.mat <- Matrix(rbind(rep(1.,n0), 
                        P.mat,
                        t(Xc)), sparse = TRUE)
  
  settings <- osqpSettings(alpha = 1.0, verbose = FALSE, warm_start = TRUE)
  for (j in 1:length(delta.v)) {
    l.vec <- c(1., rep(0.,n0), 
               colMeans(Xt) - delta.v[j] * rep(1,nX))
    u.vec <- c(1., rep(1.,n0), 
               colMeans(Xt) + delta.v[j] * rep(1,nX))
    if (j==1) {
      model <- osqp(P.mat, q.vec, A.mat, l.vec, u.vec, settings)
    } else {
      model$Update(l = l.vec, u = u.vec)
    }
    res <- model$Solve()
    if (res$info$status != "solved") {
      warning(res$info$status)
    }
    res.list[[j]] <- list(w=res$x, y_hat=sum(res$x * Yc), t=res$info$solve_time)
  }
  return(res.list)
}    

