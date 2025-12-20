
library(MASS)

n <- 300
p <- 200
d <- 2
K <- 3

X <- matrix(1, n, d)
Y <- matrix(1, n, K)
BETA <- matrix(0, p, K)
alpha <- .05
sigma <- 1

p1 <- 5
if(cs == 1){
  rho <- .5
  M <- matrix(rho, p1, p1)
  diag(M) <- 1
}else if(cs == 2){
  rho <- .6
  M <- rho^(abs(outer(1:p1, 1:p1, FUN = "-")))
}else if(cs == 3){
  M <- diag(1, p1)
}

bsign <- rep(c(1, -1), p)


theta <- 0.04
theta2 <- c(0.05, 0.03, 0.01)

if(n == 200){
  theta <- 0.08
  theta2 <- c(0.10, 0.06, 0.02)
}else if(n == 400){
  theta <- 0.04
  theta2 <- c(0.05, 0.03, 0.01)
}


G1 <- matrix(0, n, p)
for(i in 1:(p/p1)){
  G0 <- mvrnorm(n, rep(0, p1), M)
  G1[, (1 + (i - 1) * p1):(p1 * i)] <- G0
}
G <- G1

X[, 2] <- rnorm(n, 0.1 * G[,1], 1.0)

TSigma <- matrix(0.5, K, K)
diag(TSigma) <- 1
Eps <- mvrnorm(n, rep(0, K), TSigma)

Y <- X %*% matrix(1, d, K) + G %*% BETA + Eps

pval <- mvtest(Y, X, G)

pval$pval_T1
pval$pval_T2
