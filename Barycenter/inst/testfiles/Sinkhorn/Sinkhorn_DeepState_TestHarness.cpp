#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List Sinkhorn(arma::mat a, arma::mat b, arma::mat costm, double lambda, int maxIter, double tolerance);

TEST(Barycenter_deepstate_test,Sinkhorn_test){
  std::ofstream a_stream;
  std::ofstream b_stream;
  std::ofstream costm_stream;
  std::ofstream lambda_stream;
  std::ofstream maxIter_stream;
  std::ofstream tolerance_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat a  = RcppDeepState_mat();
  a_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Barycenter/inst/testfiles/Sinkhorn/inputs/a");
  a_stream << a;
  std::cout << "a values: "<< a << std::endl;
  a_stream.close();
  arma::mat b  = RcppDeepState_mat();
  b_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Barycenter/inst/testfiles/Sinkhorn/inputs/b");
  b_stream << b;
  std::cout << "b values: "<< b << std::endl;
  b_stream.close();
  arma::mat costm  = RcppDeepState_mat();
  costm_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Barycenter/inst/testfiles/Sinkhorn/inputs/costm");
  costm_stream << costm;
  std::cout << "costm values: "<< costm << std::endl;
  costm_stream.close();
  double lambda  = RcppDeepState_double();
  lambda_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Barycenter/inst/testfiles/Sinkhorn/inputs/lambda");
  lambda_stream << lambda;
  std::cout << "lambda values: "<< lambda << std::endl;
  lambda_stream.close();
  int maxIter  = RcppDeepState_int();
  maxIter_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Barycenter/inst/testfiles/Sinkhorn/inputs/maxIter");
  maxIter_stream << maxIter;
  std::cout << "maxIter values: "<< maxIter << std::endl;
  maxIter_stream.close();
  double tolerance  = RcppDeepState_double();
  tolerance_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Barycenter/inst/testfiles/Sinkhorn/inputs/tolerance");
  tolerance_stream << tolerance;
  std::cout << "tolerance values: "<< tolerance << std::endl;
  tolerance_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    Sinkhorn(a,b,costm,lambda,maxIter,tolerance);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
