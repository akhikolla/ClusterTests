#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector fastSum(NumericVector xx);

TEST(BayesESS_deepstate_test,fastSum_test){
  std::ofstream xx_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector xx  = RcppDeepState_NumericVector();
  xx_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesESS/inst/testfiles/fastSum/inputs/xx");
  xx_stream << xx;
  std::cout << "xx values: "<< xx << std::endl;
  xx_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    fastSum(xx);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
