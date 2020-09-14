

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using  Eigen::Map;
using  Eigen::VectorXd;
typedef  Map<VectorXd>  MapVecd;

// [[Rcpp::export]]
NumericVector fastSum(NumericVector xx) {
  const MapVecd x(as<MapVecd>(xx));
  NumericVector LB(1);
  LB[0] = x.sum();
  return LB;
}
