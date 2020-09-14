#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double vusC(NumericVector tt, NumericMatrix dd);

TEST(bcROCsurface_deepstate_test,vusC_test){
  std::ofstream tt_stream;
  std::ofstream dd_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector tt  = RcppDeepState_NumericVector();
  tt_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/vusC/inputs/tt");
  tt_stream << tt;
  std::cout << "tt values: "<< tt << std::endl;
  tt_stream.close();
  NumericMatrix dd  = RcppDeepState_NumericMatrix();
  dd_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/vusC/inputs/dd");
  dd_stream << dd;
  std::cout << "dd values: "<< dd << std::endl;
  dd_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    vusC(tt,dd);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
