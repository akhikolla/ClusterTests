#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double randExponential(const double lambda);

TEST(BayesSUR_deepstate_test,randExponential_test){
  std::ofstream lambda_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double lambda  = RcppDeepState_double();
  lambda_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesSUR/inst/testfiles/randExponential/inputs/lambda");
  lambda_stream << lambda;
  std::cout << "lambda values: "<< lambda << std::endl;
  lambda_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    randExponential(lambda);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
