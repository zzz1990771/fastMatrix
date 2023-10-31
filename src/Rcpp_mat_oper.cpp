#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix chol_decomp_Rcpp(NumericMatrix A) {
  int n = A.nrow();
  NumericMatrix L(n, n);
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = 0;
      for (int k = 0; k < j; ++k) {
        sum += L(i, k) * L(j, k);
      }
      if (i == j) L(i, j) = sqrt(A(i, i) - sum);
      else L(i, j) = (A(i, j) - sum) / L(j, j);
    }
  }
  return L;
}

// [[Rcpp::export]]
NumericMatrix mat_mult_Rcpp(NumericMatrix A, NumericMatrix B) {
  int n = A.nrow(), k = A.ncol(), m = B.ncol();
  NumericMatrix C(n, m);
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      double sum = 0;
      for (int x = 0; x < k; ++x) {
        sum += A(i, x) * B(x, j);
      }
      C(i, j) = sum;
    }
  }
  return C;
}