#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::cube cube_from_NumericVector(NumericVector x);

TEST(beyondWhittle_deepstate_test,cube_from_NumericVector_test){
  std::ofstream x_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector x  = RcppDeepState_NumericVector();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/cube_from_NumericVector/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    cube_from_NumericVector(x);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
