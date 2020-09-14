#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector unit_trace_x_from_phi(NumericVector phi);

TEST(beyondWhittle_deepstate_test,unit_trace_x_from_phi_test){
  std::ofstream phi_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector phi  = RcppDeepState_NumericVector();
  phi_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/unit_trace_x_from_phi/inputs/phi");
  phi_stream << phi;
  std::cout << "phi values: "<< phi << std::endl;
  phi_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    unit_trace_x_from_phi(phi);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
