#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix pacf2AR(NumericVector pacf);

TEST(beyondWhittle_deepstate_test,pacf2AR_test){
  std::ofstream pacf_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector pacf  = RcppDeepState_NumericVector();
  pacf_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/pacf2AR/inputs/pacf");
  pacf_stream << pacf;
  std::cout << "pacf values: "<< pacf << std::endl;
  pacf_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    pacf2AR(pacf);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
