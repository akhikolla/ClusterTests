#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

unsigned int randBernoulli(double pi);

TEST(BayesSUR_deepstate_test,randBernoulli_test){
  std::ofstream pi_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double pi  = RcppDeepState_double();
  pi_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesSUR/inst/testfiles/randBernoulli/inputs/pi");
  pi_stream << pi;
  std::cout << "pi values: "<< pi << std::endl;
  pi_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    randBernoulli(pi);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
