#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double acceptanceRate(NumericVector trace);

TEST(beyondWhittle_deepstate_test,acceptanceRate_test){
  std::ofstream trace_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector trace  = RcppDeepState_NumericVector();
  trace_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/acceptanceRate/inputs/trace");
  trace_stream << trace;
  std::cout << "trace values: "<< trace << std::endl;
  trace_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    acceptanceRate(trace);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
