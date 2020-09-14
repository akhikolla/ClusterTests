#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

long int factorial(int n);

TEST(BayesMallows_deepstate_test,factorial_test){
  std::ofstream n_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesMallows/inst/testfiles/factorial/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    factorial(n);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
