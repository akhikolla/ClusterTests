#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double mll_meanvar2(double x, double x2, int n);

TEST(bartBMA_deepstate_test,mll_meanvar2_test){
  std::ofstream x_stream;
  std::ofstream x2_stream;
  std::ofstream n_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double x  = RcppDeepState_double();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/mll_meanvar2/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  double x2  = RcppDeepState_double();
  x2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/mll_meanvar2/inputs/x2");
  x2_stream << x2;
  std::cout << "x2 values: "<< x2 << std::endl;
  x2_stream.close();
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/mll_meanvar2/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    mll_meanvar2(x,x2,n);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
