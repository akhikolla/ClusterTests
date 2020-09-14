#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector pFromV(NumericVector v);

TEST(beyondWhittle_deepstate_test,pFromV_test){
  std::ofstream v_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector v  = RcppDeepState_NumericVector();
  v_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/pFromV/inputs/v");
  v_stream << v;
  std::cout << "v values: "<< v << std::endl;
  v_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    pFromV(v);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
