#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double gaussian_dichotomy(int hn, NumericVector t);

TEST(binnednp_deepstate_test,gaussian_dichotomy_test){
  std::ofstream hn_stream;
  std::ofstream t_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int hn  = RcppDeepState_int();
  hn_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/gaussian_dichotomy/inputs/hn");
  hn_stream << hn;
  std::cout << "hn values: "<< hn << std::endl;
  hn_stream.close();
  NumericVector t  = RcppDeepState_NumericVector();
  t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/gaussian_dichotomy/inputs/t");
  t_stream << t;
  std::cout << "t values: "<< t << std::endl;
  t_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    gaussian_dichotomy(hn,t);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
