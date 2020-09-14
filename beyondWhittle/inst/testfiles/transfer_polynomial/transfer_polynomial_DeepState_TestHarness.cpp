#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::cx_cube transfer_polynomial(NumericVector lambda, arma::mat coef);

TEST(beyondWhittle_deepstate_test,transfer_polynomial_test){
  std::ofstream lambda_stream;
  std::ofstream coef_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector lambda  = RcppDeepState_NumericVector();
  lambda_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/transfer_polynomial/inputs/lambda");
  lambda_stream << lambda;
  std::cout << "lambda values: "<< lambda << std::endl;
  lambda_stream.close();
  arma::mat coef  = RcppDeepState_mat();
  coef_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/transfer_polynomial/inputs/coef");
  coef_stream << coef;
  std::cout << "coef values: "<< coef << std::endl;
  coef_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    transfer_polynomial(lambda,coef);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
