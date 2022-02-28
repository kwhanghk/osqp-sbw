# --------------------------------------------------------------------#
#                         example 1 (simple)                          #
# --------------------------------------------------------------------#
n=1e+07
X.1 <- runif(n,-3,3)
X.2 <- rchisq(n,1)
X.3 <- rbinom(n,1,0.5)
X <- cbind(X.1, X.2, X.3)
A <- rbinom(n,1,0.5)
Y <- (X.1 + X.2 + X.3)^2 + rnorm(n,0,1)

# a single delta value
res.0 <- osqp_sbw(X,A,Y,0.005)
res.0[[1]]$t 

# five delta values with factorization caching
delta.v = seq(.001, .05, by=.0005)
res.1 <- osqp_sbw(X,A,Y,delta.v)
cum.t <- 0
for (i in 1:length(res.1)) {
  cum.t <- cum.t + res.1[[i]]$t  
}
print(cum.t)


# --------------------------------------------------------------------#
#                  example 2 (Hainmueller, 2012)                      #
# --------------------------------------------------------------------#
n=1e+07
n_cov = 6
sig.123 <- diag(c(2,1,1))
sig.123[1,2] <- 1; sig.123[1,3] <- -1; sig.123[2,3] <- -0.5;
sig.123 <- forceSymmetric(sig.123)
beta_coef <- c(1,2,-2,-1,-0.5,1)

X.123 <- as.matrix(mvrnorm(n, mu = rep(0,3), Sigma = sig.123))
X.4 <- runif(n,-3,3)
X.5 <- rchisq(n,1)
X.6 <- rbinom(n,1,0.5)
X <- cbind(X.123, X.4, X.5, X.6)
A <- ifelse(X %*% matrix(beta_coef, ncol = 1) + rnorm(n,0,30) > 0,1,0)
Y <- (X.123[,1] + X.123[,2] + X.5)^2 + rnorm(n,0,1)

res.0 <- osqp_sbw(X,A,Y,0.005)
res.0[[1]]$t 
