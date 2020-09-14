#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat fastProd(const arma::vec & X, const arma::vec & Y){
  arma::mat out = X % Y;  
  return out;
}
