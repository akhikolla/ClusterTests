#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix matProdT_LO(const NumericMatrix X);

TEST(Benchmarking_deepstate_test,matProdT_LO_test){
  std::ofstream X_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix X  = RcppDeepState_NumericMatrix();
  X_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/matProdT_LO/inputs/X");
  X_stream << X;
  std::cout << "X values: "<< X << std::endl;
  X_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    matProdT_LO(X);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
