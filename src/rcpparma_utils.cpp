// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List rcpparma_bothproducts(const arma::colvec & x) {
    arma::mat op = x * x.t();
    double    ip = arma::as_scalar(x.t() * x);
    return Rcpp::List::create(Rcpp::Named("outer")=op,
                              Rcpp::Named("inner")=ip);
}

// [[Rcpp::export]]
arma::mat chol_decomp_RcppArma(arma::mat A) {
  arma::mat L = arma::chol(A, "lower");
  return L;
}

// [[Rcpp::export]]
arma::mat mat_mult_RcppArma(arma::mat A, arma::mat B) {
  arma::mat C = A * B;
  return C;
}

