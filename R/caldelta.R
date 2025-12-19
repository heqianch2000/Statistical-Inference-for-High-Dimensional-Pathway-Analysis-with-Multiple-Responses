#' estimate the threshold for power enhancement component.
#'
#' @param X A matrix of adjusting covariates.
#' @param G A matrix of genetic variants.
#' @param B The number of the replicates.
#' @param seed The random seed.
#' @param simple logical. If true, use the sampling method to estimate the threshold.
#' @returns returns the threshold for the power enhancement.
#' @export
caldelta <- function(X, G, B = 2000, seed = 1, simple = FALSE){

  set.seed(seed)
  n <- nrow(X)
  p <- ncol(G)

  if(!simple){
    n1 <- ceiling(n/2)

    bsf <- function(){
      idx <- sample(1:n, size = n1)
      X1 <- X[idx, ]
      G1 <- G[idx, ]
      e <- rnorm(n1)
      Px <- X1%*% solve(t(X1) %*% X1) %*% t(X1)
      IPx <- diag(1, n1) - Px
      halfA <- t(G1) %*% IPx
      vecD <- 1 / diag(halfA %*% t(halfA))
      halfA2 <- sweep(halfA, 1, sqrt(vecD), FUN = "*")
      betam <- halfA2 %*% e
      max(abs(betam))
    }
    tmp <- replicate(B, bsf())
  }else{
    tmp <- replicate(B, max(abs(rnorm(p))))
  }

  return(max(tmp))

}
