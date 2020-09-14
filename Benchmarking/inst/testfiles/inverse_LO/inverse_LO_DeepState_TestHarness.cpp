#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix inverse_LO(const NumericMatrix L);

TEST(Benchmarking_deepstate_test,inverse_LO_test){
  std::ofstream L_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix L  = RcppDeepState_NumericMatrix();
  L_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/inverse_LO/inputs/L");
  L_stream << L;
  std::cout << "L values: "<< L << std::endl;
  L_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    inverse_LO(L);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
