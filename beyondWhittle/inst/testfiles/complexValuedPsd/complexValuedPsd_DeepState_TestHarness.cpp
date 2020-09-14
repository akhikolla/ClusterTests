#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::cx_cube complexValuedPsd(NumericVector f_);

TEST(beyondWhittle_deepstate_test,complexValuedPsd_test){
  std::ofstream f__stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector f_  = RcppDeepState_NumericVector();
  f__stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/complexValuedPsd/inputs/f_");
  f__stream << f_;
  std::cout << "f_ values: "<< f_ << std::endl;
  f__stream.close();
  std::cout << "input ends" << std::endl;
  try{
    complexValuedPsd(f_);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
