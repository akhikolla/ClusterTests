#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double varFh(double x, int n, NumericVector t, NumericVector p, double h);

TEST(binnednp_deepstate_test,varFh_test){
  std::ofstream x_stream;
  std::ofstream n_stream;
  std::ofstream t_stream;
  std::ofstream p_stream;
  std::ofstream h_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double x  = RcppDeepState_double();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/varFh/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/varFh/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  NumericVector t  = RcppDeepState_NumericVector();
  t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/varFh/inputs/t");
  t_stream << t;
  std::cout << "t values: "<< t << std::endl;
  t_stream.close();
  NumericVector p  = RcppDeepState_NumericVector();
  p_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/varFh/inputs/p");
  p_stream << p;
  std::cout << "p values: "<< p << std::endl;
  p_stream.close();
  double h  = RcppDeepState_double();
  h_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/varFh/inputs/h");
  h_stream << h;
  std::cout << "h values: "<< h << std::endl;
  h_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    varFh(x,n,t,p,h);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
