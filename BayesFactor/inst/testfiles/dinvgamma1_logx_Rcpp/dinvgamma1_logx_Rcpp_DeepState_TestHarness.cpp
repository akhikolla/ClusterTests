#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double dinvgamma1_logx_Rcpp(const double x, const double a, const double b);

TEST(BayesFactor_deepstate_test,dinvgamma1_logx_Rcpp_test){
  std::ofstream x_stream;
  std::ofstream a_stream;
  std::ofstream b_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double x  = RcppDeepState_double();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesFactor/inst/testfiles/dinvgamma1_logx_Rcpp/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  double a  = RcppDeepState_double();
  a_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesFactor/inst/testfiles/dinvgamma1_logx_Rcpp/inputs/a");
  a_stream << a;
  std::cout << "a values: "<< a << std::endl;
  a_stream.close();
  double b  = RcppDeepState_double();
  b_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesFactor/inst/testfiles/dinvgamma1_logx_Rcpp/inputs/b");
  b_stream << b;
  std::cout << "b values: "<< b << std::endl;
  b_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    dinvgamma1_logx_Rcpp(x,a,b);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
