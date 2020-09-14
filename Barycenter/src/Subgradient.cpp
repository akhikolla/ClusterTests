// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec Subgradient(arma::vec a, arma::vec b, arma::mat M, double lambda, int maxIter = 1000, double tolerance = 0.0001) {

  // Transforming Input, i.e. calculating the kernel

  arma::mat K = exp(-lambda*M);
  arma::mat ainvK = diagmat(1/a) * K;


  //Initialize u, v and the vector next for the stopping criterion
  arma::vec u(a.n_rows);
  u = u.ones()/a.n_rows;
  arma::vec next;

  for(int i=0;i<maxIter;i++){

    u = 1 / (ainvK * (b / (K.t()*u)));  //Sinkhorn`s update

    //The stopping criterion is checked every 20th step, i.e. if u and the next update of u (called next) differ in 2-norm only by a tolerance
    if(i % 20 == 0){
      next = 1 / (ainvK * (b / (K.t()*u)));
      double Criterion = norm(abs(next-u));
      if(Criterion < tolerance){
        u = next;
        break;
      }
      else{
        u = next;
      }
    }

  }

  //Computation of the dual solution (c.f. Proposition 2, Fast Computation of Wasserstein Barycenters)

  arma::vec alpha = (1 / lambda) * log(u) - accu(log(u))/(lambda*(u.n_rows));

  return(alpha);
}
