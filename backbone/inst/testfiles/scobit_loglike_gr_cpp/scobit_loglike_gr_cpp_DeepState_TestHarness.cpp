#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector scobit_loglike_gr_cpp(NumericVector x1, NumericVector x2, NumericVector y, NumericVector params);

TEST(backbone_deepstate_test,scobit_loglike_gr_cpp_test){
  std::ofstream x1_stream;
  std::ofstream x2_stream;
  std::ofstream y_stream;
  std::ofstream params_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector x1  = RcppDeepState_NumericVector();
  x1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/backbone/inst/testfiles/scobit_loglike_gr_cpp/inputs/x1");
  x1_stream << x1;
  std::cout << "x1 values: "<< x1 << std::endl;
  x1_stream.close();
  NumericVector x2  = RcppDeepState_NumericVector();
  x2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/backbone/inst/testfiles/scobit_loglike_gr_cpp/inputs/x2");
  x2_stream << x2;
  std::cout << "x2 values: "<< x2 << std::endl;
  x2_stream.close();
  NumericVector y  = RcppDeepState_NumericVector();
  y_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/backbone/inst/testfiles/scobit_loglike_gr_cpp/inputs/y");
  y_stream << y;
  std::cout << "y values: "<< y << std::endl;
  y_stream.close();
  NumericVector params  = RcppDeepState_NumericVector();
  params_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/backbone/inst/testfiles/scobit_loglike_gr_cpp/inputs/params");
  params_stream << params;
  std::cout << "params values: "<< params << std::endl;
  params_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    scobit_loglike_gr_cpp(x1,x2,y,params);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
