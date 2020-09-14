#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix chol_LO(const NumericMatrix A);

TEST(Benchmarking_deepstate_test,chol_LO_test){
  std::ofstream A_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix A  = RcppDeepState_NumericMatrix();
  A_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/chol_LO/inputs/A");
  A_stream << A;
  std::cout << "A values: "<< A << std::endl;
  A_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    chol_LO(A);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
