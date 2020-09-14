#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double update_sigma(double a1, double b, NumericVector resids, int n);

TEST(bartBMA_deepstate_test,update_sigma_test){
  std::ofstream a1_stream;
  std::ofstream b_stream;
  std::ofstream resids_stream;
  std::ofstream n_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double a1  = RcppDeepState_double();
  a1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/update_sigma/inputs/a1");
  a1_stream << a1;
  std::cout << "a1 values: "<< a1 << std::endl;
  a1_stream.close();
  double b  = RcppDeepState_double();
  b_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/update_sigma/inputs/b");
  b_stream << b;
  std::cout << "b values: "<< b << std::endl;
  b_stream.close();
  NumericVector resids  = RcppDeepState_NumericVector();
  resids_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/update_sigma/inputs/resids");
  resids_stream << resids;
  std::cout << "resids values: "<< resids << std::endl;
  resids_stream.close();
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/update_sigma/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    update_sigma(a1,b,resids,n);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
