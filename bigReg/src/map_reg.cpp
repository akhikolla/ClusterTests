/*
  Date: 2016-07-23
  Author: Chibisi Chima-Okereke
*/
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <fstream>

using namespace Rcpp;
using namespace arma;
using namespace std;
//extern "C" void openblas_set_num_threads(int);
//openblas_set_num_threads(4);

/*********************************************************************************************************/

// Function to convert an SEXP matrix to an armadillo mat
arma::mat armaMat(SEXP Ms)
{
  // Convert R matrix to Rcpp Matrix
  Rcpp::NumericMatrix Mr(Ms);
  // The number of rwos and columns
  int nra = Mr.nrow(), nca = Mr.ncol();
  // Convert to armadillo matrix
  arma::mat M(Mr.begin(), nra, nca, false);
  return M;
}

// Function to convert an SEXP vector to an armadillo vec
arma::vec armaVec(SEXP Vs)
{
  Rcpp::NumericVector Vr(Vs);
  int len = Vr.size();
  arma::vec V(Vr.begin(), len, false);
  return V;
}

//' @title Calculation of iterative regression components
//' @description Calculation of iterative regression components
//' @param X design matrix
//' @param y dependent variable
//' @param W weights
//' @export
//[[Rcpp::export]]
Rcpp::List XWXMatrixW(arma::mat X, arma::mat y, arma::vec W)
{
  //Ensure positive non-zero weights
  W = abs(W) + DBL_EPSILON;
  
  mat XW = X;
  XW.each_col() %= W;
  inplace_strans(XW);
  mat XWX = XW * X;
  mat XWy = XW * y;
  mat yW = y;
  yW.each_col() %= W;
  mat yWy = yW.st() * y;
  
  List out;
  out["XWX"] = XWX;
  out["XWy"] = XWy;
  out["yWy"] = yWy;
  out["nobs"] = X.n_rows;
  out["wt"] = sum(log(W));
  out["ySum"] = sum(yW);
  return out;
}

// Singular value decomposition
//' @title Singular value decomposition of the aggregated list from XWXMatrix(W) functions
//' @description Singular value decomposition of the aggregated list from XWXMatrix(W) functions
//' @param out list containing requisite computed values
//' @param epsilon either machine epsilon or user depermined epsilon
//' @export
// [[Rcpp::export]]
Rcpp::List SVD(Rcpp::List out, double epsilon)
{
  arma::mat U;
  arma::vec s;
  arma::mat V;

  arma::mat X = armaMat(out["XWX"]);
  arma::mat Y = armaMat(out["yWy"]);
  int nobs = out["nobs"];
  arma::svd(U, s, V, X);
  int tol = epsilon*s.n_elem*s.max();
  arma::umat aliased = (s <= tol);
  arma::mat covUnscaled = V * diagmat(1/s) * U.st();
  arma::mat XWy = out["XWy"];
  arma::mat beta = covUnscaled * XWy;
  out["XWy"] = XWy;
  out["cov.unscaled"] = covUnscaled;
  out["coefficients"] = beta;
  out["aliased"] = aliased;
  arma::umat rankM = sum(aliased == 0);
  int rank = rankM(0, 0);
  out["rank"] = rank;
  arma::mat mu = X*beta;
  arma::mat numerat = Y - 2*XWy.st()*beta + mu.st()*beta;
  double var = (numerat(0,0)/(nobs - rank));
  out["var"] = var;
  out["cov"] = var * covUnscaled;
  out["deviance"] = numerat(0,0);
  return out;
}

//' @title Calculation of iterative regression components
//' @description Calculation of iterative regression components
//' @param X design matrix
//' @param y dependent variable
//' @export
// [[Rcpp::export]]
Rcpp::List XWXMatrix(arma::mat X, arma::mat y)
{
  arma::mat XW = X;
  inplace_strans(XW);
  arma::mat XWX = XW * X;
  arma::mat XWy = XW * y;
  arma::mat yW = y;
  arma::mat yWy = yW.st() * y;
  
  Rcpp::List out;
  out["XWX"] = XWX;
  out["XWy"] = XWy;
  out["yWy"] = yWy;
  out["nobs"] = X.n_rows;
  out["wt"] = 0;
  out["ySum"] = sum(yW);
  return out;
}

/*********************************************************************************************************/


