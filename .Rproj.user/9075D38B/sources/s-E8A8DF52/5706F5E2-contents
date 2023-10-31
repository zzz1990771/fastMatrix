#' Generate Random Multivariate Normal Samples with Exchangeable Covariance Structure
#'
#' This function generates random samples from a multivariate normal distribution with
#' a given mean vector and an exchangeable covariance structure.
#'
#' @param n Integer. The number of samples to generate.
#' @param mu Numeric vector. The mean vector for the multivariate normal distribution.
#' @param p Integer. The dimension of the multivariate normal distribution.
#' @param rho Numeric. The correlation parameter for the exchangeable covariance structure.
#'
#' @return A numeric matrix with 'n' rows and 'p' columns, where each row represents a
#'   random sample from the specified multivariate normal distribution.
#'
#' @examples
#' \dontrun{
#'   n <- 100
#'   mu <- c(0, 0)
#'   p <- 2
#'   rho <- 0.5
#'   samples <- rmvnorm_tianhai_ex(n, mu, p, rho)
#' }
#'
#' @seealso \code{\link{fast_ex_chol}}, \code{\link{mat_mult_RcppArma}}
#'
#' @export
rmvnorm_tianhai_ex <- function(n, mu, p, rho) {
  # Cholesky decomposition
  L <- fast_ex_chol(p,rho)
  
  # Generate standard normal random variables
  Z <- matrix(rnorm(n * p), n, p)
  
  # Transform samples
  X <- fastMatrix::mat_mult_RcppArma(Z,t(L)) + matrix(mu, n, p, byrow = TRUE)
  
  return(X)
}

#' Generate Cholesky Decomposition for an Exchangeable Covariance Structure
#'
#' This function computes the Cholesky decomposition of a covariance matrix with
#' an exchangeable structure. The matrix is of dimension 'p x p'.
#'
#' @param p Integer. The dimension of the covariance matrix.
#' @param rho Numeric. The correlation parameter for the exchangeable covariance structure.
#'
#' @return A numeric matrix of dimension 'p x p' representing the Cholesky decomposition
#' of the exchangeable covariance matrix.
#'
#' @examples
#' \dontrun{
#'   p <- 5
#'   rho <- 0.5
#'   L <- fast_ex_chol(p, rho)
#' }
#'
#' @export
fast_ex_chol <- function(p,rho){
  diag <- numeric(p)
  tril <- numeric(p)
  
  diag[1] <- 1
  tril[1] <- rho
  
  for (j in 2:p) {
    diag[j] <- sqrt(diag[j - 1]^2 - tril[j - 1]^2)
    tril[j] <- (rho - 1) / diag[j] + diag[j]
  }
  tril[p] <- 0
  L = matrix(tril,p,p,byrow=TRUE)
  
  
  diag(L) = diag
  
  L[upper.tri(L,diag = FALSE)] = 0
  return(L)
}
