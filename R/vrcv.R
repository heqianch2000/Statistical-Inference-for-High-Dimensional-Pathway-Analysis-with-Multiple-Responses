#' estimate the error variance in linear regression using refitted cross validation.
#'
#' @param Y The response vector.
#' @param X A matrix of adjusting covariates.
#' @param Num The number of the replicates.
#' @param sis_prop The proportion in SIS to select relevant features.
#' @param seed The random seed.
#' @returns returns the estimated error variance.
#' @export
vrcv <- function(Y, X, Num = 10, sis_prop = 0.5, seed = NULL){

  n <- length(Y)
  p <- ncol(X)

  if(!is.null(seed) & is(seed, "numeric")){
    set.seed(seed)
  }

  n1 <- floor(n / 2)

  res <- replicate(Num, {
    cv1 <- sample(n, n1)
    X1 <- X[cv1, ]
    Y1 <- Y[cv1]
    X2 <- X[-cv1, ]
    Y2 <- Y[-cv1]


    # SIS
    cor1 <- cor(X1, Y1)
    cor2 <- cor(X2, Y2)

    n2 <- n1 * sis_prop

    ind1 <- order(abs(cor1), decreasing = TRUE)[1:n2]
    ind2 <- order(abs(cor2), decreasing = TRUE)[1:n2]

    X12 <- cbind(1, X1[, ind2])
    X22 <- cbind(1, X2[, ind1])

    # fit1 <- lm(Y1 ~ X12 - 1)
    # fit2 <- lm(Y2 ~ X22 - 1)
    #
    # s1 <- (summary(fit1)$sigma)^2
    # s2 <- (summary(fit2)$sigma)^2

    IP1 <- diag(1, n1) - X12 %*% solve(t(X12) %*% X12) %*% t(X12)
    IP2 <- diag(1, n - n1) - X22 %*% solve(t(X22) %*% X22) %*% t(X22)

    s1 <- (t(Y1) %*% IP1 %*% Y1) / (n1 - n2 - 1)
    s2 <- (t(Y2) %*% IP2 %*% Y2) / (n - n1 - n2 - 1)

    (s1 + s2) / 2

  })

  mean(res)

}
