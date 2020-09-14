#include <vector>
#include <random>
#include <ctime>
#include <Rcpp.h>
#include <iostream>
#include "model.h"

using namespace std;
using namespace Rcpp;

mt19937 generator(time(0));

//[[Rcpp::export]]
NumericMatrix use_CANL(const NumericMatrix mat, const int N, const double w0, const double k, const double eps) {

  vector<vector<double>> cpp_params;
  int nrows = mat.rows();
  for(int i = 0; i < nrows; i++) {
    NumericVector tmpVec = mat(i, _);
    vector<double> tmp = as<vector<double>>(tmpVec);
    cpp_params.push_back(tmp);
  }

  int n = nrows;

  vector<vector<double>> t(n, vector<double>(N));
  vector<double> tmp(N);

  for(int i = 0; i < n; ++i) {
    tmp = CANL(cpp_params[i], N, w0, k, 0.5);
    t[i] = tmp;
  }
  //t is n x N

  NumericMatrix outMat(n, N);
  for(int i = 0; i < n; i++) {
    // Copy by rows of t
    NumericVector tmpOut(t[i].begin(), t[i].end());
    outMat(i, _) = tmpOut;
  }
  return outMat;
}
