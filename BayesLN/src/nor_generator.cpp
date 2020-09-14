#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>



arma::mat mvrnorm(arma::vec mu, arma::mat sigma) {
  int p=sigma.n_cols;
  arma::vec eps = Rcpp::rnorm(p);
  arma::mat output;
    output = arma::chol(sigma) * eps + mu;

  return output;

}





arma::mat mvrnorm_Ind(arma::vec mu, double sigma2) {
  int p=mu.n_elem;
  double sigma = sqrt(sigma2);
  arma::vec eps = Rcpp::rnorm(p);
  arma::mat output;
  output = sigma * arma::eye(p,p) * eps + mu;

  return output;

}
