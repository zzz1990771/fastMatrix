#' Generate Random Multivariate Normal Samples with R functions
#'
#' This function generates random samples from a multivariate normal distribution with
#' a given mean vector and covariance matrix. It is more efficient compared to 
#' the MASS package.
#'
#' @param n Integer. The number of samples to generate.
#' @param mu Numeric vector. The mean vector for the multivariate normal distribution.
#' @param Sigma Numeric matrix. The covariance matrix for the multivariate normal distribution.
#'
#' @return A numeric matrix with 'n' rows and 'p' columns, where each row represents a
#'   random sample from the specified multivariate normal distribution.
#'
#' @examples
#' \dontrun{
#'   n <- 100
#'   mu <- c(0, 0)
#'   Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#'   samples <- rmvnorm_R(n, mu, Sigma)
#' }
#'
#' @seealso \code{\link{chol}}
#'
#' @export
rmvnorm_R <- function(n, mu, Sigma) {
  p = length(mu)
  
  # Cholesky decomposition
  L <- chol(Sigma)
  
  # Generate standard normal random variables
  Z <- matrix(rnorm(n * p), n, p)
  
  # Transform samples
  X <- Z %*% (L) + matrix(mu, n, p, byrow = TRUE)
  
  return(X)
}


#' Generate Random Multivariate Normal Samples with CPP Armadillo functions
#'
#' This function generates random samples from a multivariate normal distribution with
#' a given mean vector and covariance matrix.
#'
#' @param n Integer. The number of samples to generate.
#' @param mu Numeric vector. The mean vector for the multivariate normal distribution.
#' @param Sigma Numeric matrix. The covariance matrix for the multivariate normal distribution.
#'
#' @return A numeric matrix with 'n' rows and 'p' columns, where each row represents a
#'   random sample from the specified multivariate normal distribution.
#'
#' @examples
#' \dontrun{
#'   n <- 100
#'   mu <- c(0, 0)
#'   Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#'   samples <- rmvnorm_tianhai_cppArma(n, mu, Sigma)
#' }
#'
#' @seealso \code{\link{chol}}
#'
#' @export
rmvnorm_tianhai_cppArma <- function(n, mu, Sigma) {
  p = length(mu)
  
  # Cholesky decomposition
  L <- chol_decomp_RcppArma(Sigma)
  
  # Generate standard normal random variables
  Z <- matrix(rnorm(n * p), n, p)
  
  # Transform samples
  X <- mat_mult_RcppArma(Z,t(L)) + matrix(mu, n, p, byrow = TRUE)
  
  return(X)
}


#' Generate Random Multivariate Normal Samples with CPP functions
#'
#' This function generates random samples from a multivariate normal distribution with
#' a given mean vector and covariance matrix.
#'
#' @param n Integer. The number of samples to generate.
#' @param mu Numeric vector. The mean vector for the multivariate normal distribution.
#' @param Sigma Numeric matrix. The covariance matrix for the multivariate normal distribution.
#'
#' @return A numeric matrix with 'n' rows and 'p' columns, where each row represents a
#'   random sample from the specified multivariate normal distribution.
#'
#' @examples
#' \dontrun{
#'   n <- 100
#'   mu <- c(0, 0)
#'   Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#'   samples <- rmvnorm_tianhai_cpp(n, mu, Sigma)
#' }
#'
#' @seealso \code{\link{chol}}
#'
#' @export
rmvnorm_tianhai_cpp <- function(n, mu, Sigma) {
  p = length(mu)
  
  # Cholesky decomposition
  L <- chol_decomp_Rcpp(Sigma)
  
  # Generate standard normal random variables
  Z <- matrix(rnorm(n * p), n, p)
  
  # Transform samples
  X <- mat_mult_Rcpp(Z,t(L)) + matrix(mu, n, p, byrow = TRUE)
  
  return(X)
}