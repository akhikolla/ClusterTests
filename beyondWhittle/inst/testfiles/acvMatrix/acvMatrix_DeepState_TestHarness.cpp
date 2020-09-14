#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix acvMatrix(NumericVector acv);

TEST(beyondWhittle_deepstate_test,acvMatrix_test){
  std::ofstream acv_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector acv  = RcppDeepState_NumericVector();
  acv_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/acvMatrix/inputs/acv");
  acv_stream << acv;
  std::cout << "acv values: "<< acv << std::endl;
  acv_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    acvMatrix(acv);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
