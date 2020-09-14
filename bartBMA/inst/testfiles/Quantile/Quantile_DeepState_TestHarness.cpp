#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::NumericVector Quantile(Rcpp::NumericVector x, Rcpp::NumericVector probs);

TEST(bartBMA_deepstate_test,Quantile_test){
  std::ofstream x_stream;
  std::ofstream probs_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector x  = RcppDeepState_NumericVector();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/Quantile/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  NumericVector probs  = RcppDeepState_NumericVector();
  probs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/Quantile/inputs/probs");
  probs_stream << probs;
  std::cout << "probs values: "<< probs << std::endl;
  probs_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    Quantile(x,probs);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
