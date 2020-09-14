#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double secondKindStirlingNumber(int n, int k);

TEST(bartBMA_deepstate_test,secondKindStirlingNumber_test){
  std::ofstream n_stream;
  std::ofstream k_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/secondKindStirlingNumber/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  int k  = RcppDeepState_int();
  k_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/secondKindStirlingNumber/inputs/k");
  k_stream << k;
  std::cout << "k values: "<< k << std::endl;
  k_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    secondKindStirlingNumber(n,k);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
