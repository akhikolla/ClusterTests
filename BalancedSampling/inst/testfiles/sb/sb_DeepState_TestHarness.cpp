#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double sb(NumericVector p, NumericMatrix x, NumericVector s);

TEST(BalancedSampling_deepstate_test,sb_test){
  std::ofstream p_stream;
  std::ofstream x_stream;
  std::ofstream s_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector p  = RcppDeepState_NumericVector();
  p_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/sb/inputs/p");
  p_stream << p;
  std::cout << "p values: "<< p << std::endl;
  p_stream.close();
  NumericMatrix x  = RcppDeepState_NumericMatrix();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/sb/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  NumericVector s  = RcppDeepState_NumericVector();
  s_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/sb/inputs/s");
  s_stream << s;
  std::cout << "s values: "<< s << std::endl;
  s_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    sb(p,x,s);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
