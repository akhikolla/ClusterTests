#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector solve_LO(const NumericMatrix L, const NumericVector d);

TEST(Benchmarking_deepstate_test,solve_LO_test){
  std::ofstream L_stream;
  std::ofstream d_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix L  = RcppDeepState_NumericMatrix();
  L_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/solve_LO/inputs/L");
  L_stream << L;
  std::cout << "L values: "<< L << std::endl;
  L_stream.close();
  NumericVector d  = RcppDeepState_NumericVector();
  d_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/solve_LO/inputs/d");
  d_stream << d;
  std::cout << "d values: "<< d << std::endl;
  d_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    solve_LO(L,d);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
