#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//some support functions for the main function Greenkhorn
namespace supporters{

arma::mat roundPolytope(arma::mat X, arma::colvec r, arma::rowvec c){
  arma::mat A = X;
  int n = size(A,0);
  double scaling = 0;
  arma::colvec r_A = sum(A,1);
  for(int i = 0; i < n; i++) {
    scaling = std::min(1.0,r[i]/r_A[i]);
    A.row(i) = scaling*A.row(i);
  }
  arma::rowvec c_A = sum(A,0); // additionally added
  n = size(A,1);
  for(int j = 0; j < n; j++) {
    scaling = std::min(1.0,c[j]/c_A[j]);
    A.col(j) = scaling*A.col(j);
  }
  r_A = sum(A,1);
  c_A = sum(A,0);
  arma::colvec err_r = r_A - r;
  double error = sum(abs(err_r));
  if(error != 0){
    arma::rowvec err_c = c_A - c;
    A = A + (err_r*err_c)/sum(abs(err_r));
  }
  return A;
}

int vecmaxInd(arma::vec x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::max_element(x.begin(), x.end());
  // we want the value so dereference
  return it - x.begin();
}

double frobinnerproduct(arma::mat P, arma::mat C){
  return sum(sum(P % C));
}

arma::mat artificalRows(arma::mat A, arma::uvec ids){
  arma::rowvec add = arma::zeros(1,A.n_cols);
  for(int i = 0; i < ids.n_rows; i++){
    A.insert_rows(ids[i],1);
  }
  return A;
}
arma::mat artificalCols(arma::mat A, arma::uvec ids){
  arma::rowvec add = arma::zeros(1,A.n_rows);
  for(int i = 0; i < ids.n_rows; i++){
    A.insert_cols(ids[i],1);
  }
  return A;
}

}

//[[Rcpp::export]]
Rcpp::List Greenkhorn(
    arma::colvec r, arma::rowvec c, arma::mat costm, double lambda = 1,  int maxIter = 10000, double tolerance = 10^(-8)
){

  //initialize the dual solution for r
  arma::vec dual_r;
  dual_r = dual_r.ones(r.n_rows);

  //check for full support of the measures. If not we modify our transportation problem to a smaller instance.
  arma::uvec ids_rows;
  arma::uvec ids_cols;
  arma::uvec ids_rows_temp;
  arma::uvec ids_cols_temp;
  bool ZeroValues_r = false;
  bool ZeroValues_c = false;
  if(all(r)==false){
    ZeroValues_r = true;
    ids_rows = find(r == 0); // Find indices without support
    ids_cols_temp = find(r>0); // Find support for r
    costm = costm.rows(ids_cols_temp);
    r = r.rows(ids_cols_temp);
  }
  if(all(c)==false){
    ZeroValues_c = true;
    ids_cols = find(c == 0); // Find indices without support
    ids_rows_temp = find(c>0); // Find support of c
    costm = costm.cols(ids_rows_temp);
    c = c.cols(ids_rows_temp);
  }

  //kernel matrix to project
  arma::mat P = exp(-(1/lambda)*costm);

  //compute full row and column marginals once
  arma::colvec r_P = sum(P,1);
  arma::rowvec c_P = sum(P,0);

  //compute gains for each row and column
  arma::colvec r_gain = r_P - r + r % log(r/r_P);
  arma::rowvec c_gain = c_P - c + c % log(c/c_P);


  //marginal error
  double error = norm(r_P-r,1) + norm(c_P-c,1);

  //main iteration

  for(int t = 0; t < maxIter; t++){
    //find row or column with maximum gain
    double r_gain_max = max(r_gain); int i = supporters::vecmaxInd(r_gain);
    double c_gain_max = max(c_gain); int j = supporters::vecmaxInd(c_gain.t());
    double scaling ;
    if(r_gain_max > c_gain_max){
      //update row i
      scaling = r[i]/r_P[i];
      dual_r[i] = dual_r[i]*scaling;
      arma::rowvec old_row = P.row(i);
      arma::rowvec new_row = scaling * old_row;
      P.row(i) = new_row;
      //renormalize
      P = P/sum(sum(P));
      //compute full row and column marginals
      r_P = sum(P,1);
      c_P = sum(P,0);
      //compute gais for each row and column
      r_gain = r_P - r + r % log(r/r_P);
      c_gain = c_P - c + c % log(c/c_P);
      //checking marginals in order to termiante
      error = norm(r_P-r,1) + norm(c_P-c,1);
      if(error < tolerance){
        t = maxIter;
      }
    }
    else{
      //update column j
      scaling = c[j]/c_P[j];
      arma::colvec old_col = P.col(j);
      arma::colvec new_col = scaling * old_col;
      P.col(j) = new_col;
      //renormalize
      P = P/sum(sum(P));
      //compute full row and column marginals
      r_P = sum(P,1);
      c_P = sum(P,0);
      //compute gais for each row and column
      r_gain = r_P - r + r % log(r/r_P);
      c_gain = c_P - c + c % log(c/c_P);
      //checking marginals in order to termiante
      error = norm(r_P-r,1) + norm(c_P-c,1);
      if(error < tolerance){
        t = maxIter;
      }
    }
  }

  //Rounding to a feasible point onto the transport polytope
  supporters::roundPolytope(P,r,c);
  //Calculating the Sinkhorn distance
  double sinkhorndistance = supporters::frobinnerproduct(P,costm);

  //We need to add artificial zeros in case r or c originally do not come with full support
  if(ZeroValues_c==true){
    P = supporters::artificalCols(P,ids_cols);
  }
  if(ZeroValues_r==true){
    P = supporters::artificalRows(P,ids_rows);
  }

  return Rcpp::List::create(
    Rcpp::Named("Transportplan") = P,
    Rcpp::Named("Distance") = sinkhorndistance,
    Rcpp::Named("Dual") = lambda * (log(dual_r) - accu(log(dual_r))/(dual_r.n_rows))
  );
}


