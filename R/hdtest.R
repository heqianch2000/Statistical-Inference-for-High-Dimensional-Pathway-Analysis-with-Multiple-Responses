#' High-dimensional association test with one response
#' @description
#' Test the association between the outcome Y and genetic variants G under high dimensions.
#'
#' @param Y The response.
#' @param X A matrix of adjusting covariates.
#' @param G A matrix of genetic variants.
#' @param sigma2 The estimated error variance. If null, it will be
#'  estimated by refitted cross validation.
#' @param method the method to test the association. The default method 'high' accounts
#'  the accumulated estimation error under high dimensions.
#' @param ape a logical value indicating whether the power enhancement component is used.
#' @param delta the threshold for the power enhancement.
#' @param boncor a logical value indicating whether the p-value from the Bonferroni
#' correction should be reported.
#' @param edgecor a logical value indicating whether the Edgeworth expansion is used
#' for extreme significance levels.
#' @returns A list containing the test statistics and the associated p-values.
#' \item{TH}{the observed TH test statistic}
#' \item{ph}{the p-value of the TH test}
#' \item{sd1}{the standard deviation of the test statistic}
#' \item{Tpe}{the observed test statistic with the power enhancement}
#' \item{phe}{the p-value of the test with the power enhancement}
#' \item{pb}{the p-value of the marginal models after Bonferroni correction}
#' \item{ph_edge}{the p-value of the TH test using Edgeworth expansion}
#' \item{phe_edge}{the p-value of the test with the power enhancement component using Edgeworth expansion}
#' \item{diff}{difference of p^-1||Sigma||_F^2 - p/n to verify a condition in the paper}
#' @export
hdtest <- function(Y, X, G, sigma2 = NULL, method = "high",
                  ape = TRUE, delta = NULL,
                  boncor = TRUE, edgecor = FALSE, ...){

  n <- nrow(X)
  d <- ncol(X)
  p <- ncol(G)

  if(is.null(sigma2)){
    fit <- lm(Y ~ X - 1)
    sigma2 <- (summary(fit)$sigma)^2
  }

  Px <- X %*% solve(t(X) %*% X) %*% t(X)
  IPx <- diag(1, n) - Px
  halfA <- t(G) %*% IPx
  vecD <- 1 / diag(halfA %*% t(halfA))
  halfA2 <- sweep(halfA, 1, vecD, FUN = "*")
  A <- t(halfA) %*% halfA2

  diff <- NA
  if(method == "low"){
    A1 <- A
    Q <- c(t(Y) %*% A1 %*% Y)
    sd1 <- sqrt(2 * sum(A1^2))
    TH <- (Q - p * sigma2)/(sd1 * sigma2)
  }else if(method == "high"){
    A1 <- A - p / (n - d) * IPx
    fnorm <- sum(A1^2)
    sd1 <- sqrt(2 * fnorm)
    Q1 <- c(t(Y) %*% A1 %*% Y)
    TH <- Q1 / (sd1 * sigma2)
    diff <- fnorm / p
  }
  ph <- 2 * pnorm(abs(TH), lower.tail = FALSE) # pvalue for TH


  Tpe <- NA
  phe <- NA
  pb <- NA
  if(ape){

    if(is.null(delta)){
      delta <- log(log(n)) * sqrt(log(p))
    }

    # marginal coefficients and sds
    betam <- halfA2 %*% Y
    sdm <- sqrt(sigma2 * vecD)

    hatS <- which(abs(betam) > sdm * delta)

    if(length(hatS) >= 1){
      T0 <- tsign(TH) * sqrt(p) * sum( betam[hatS]^2 / (sdm[hatS])^2 )
      Tpe <- TH + T0
    }else{
      Tpe <- TH
    }

    if(boncor){
      Tb <- max(abs(betam) / sdm)
      pb <- 2 * pnorm(Tb, lower.tail = FALSE) * p
    }

    phe <- 2 * pnorm(abs(Tpe), lower.tail = FALSE)

  }

  # edgeworth correction
  ph_edge <- rep(NA, 2)
  phe_edge <- rep(NA, 2)
  if(edgecor){

    THa <- abs(TH)
    M1 <- sum(diag(A1))
    M2 <- 2 * sum(A1^2)
    M3 <- 8 * sum(diag(A1 %*% A1 %*% A1))
    M4 <- 48 * sum(diag(A1 %*% A1 %*% A1 %*% A1))

    E1 <- pnorm(THa) + M3 * (1 - THa^2) * dnorm(THa) / 6 / (sqrt(M2))^3
    ph_edge[1] <- 2 * (1 - E1)

    E2 <- E1 - M4 * (THa^3 - 3 * THa) * dnorm(THa) / 24 / M2^2 -
      M3^2 * (THa^5 - 10 * THa^3 + 15 * THa) * dnorm(THa) / 72 / M2^3

    ph_edge[2] <- 2 * (1 - E2)

    if(ape){
      if(length(hatS) >= 1){
        Tpea <- abs(Tpe)
        Ee1 <- pnorm(Tpea) + M3 * (1 - Tpea^2) * dnorm(Tpea) / 6 / (sqrt(M2))^3
        phe_edge[1] <- 2 * (1 - Ee1)

        Ee2 <- Ee1 - M4 * (Tpea^3 - 3 * Tpea) * dnorm(Tpea) / 24 / M2^2 -
          M3^2 * (Tpea^5 - 10 * Tpea^3 + 15 * Tpea) * dnorm(Tpea) / 72 / M2^3

        phe_edge[2] <- 2 * (1 - Ee2)


      }else{
        phe_edge <- ph_edge
      }
    }


  }


  res <- list(TH = TH, ph = ph, sd1 = sd1, Tpe = Tpe, phe = phe,
              pb = pb, ph_edge = ph_edge, phe_edge = phe_edge,
              diff = diff)

  res
}
