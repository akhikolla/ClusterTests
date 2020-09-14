#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector calcw_cpp(NumericVector xb, NumericVector y);

TEST(binnednp_deepstate_test,calcw_cpp_test){
  std::ofstream xb_stream;
  std::ofstream y_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector xb  = RcppDeepState_NumericVector();
  xb_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/calcw_cpp/inputs/xb");
  xb_stream << xb;
  std::cout << "xb values: "<< xb << std::endl;
  xb_stream.close();
  NumericVector y  = RcppDeepState_NumericVector();
  y_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/calcw_cpp/inputs/y");
  y_stream << y;
  std::cout << "y values: "<< y << std::endl;
  y_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    calcw_cpp(xb,y);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
