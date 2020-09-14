#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

int randIntUniform(const int a, const int b);

TEST(BayesSUR_deepstate_test,randIntUniform_test){
  std::ofstream a_stream;
  std::ofstream b_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int a  = RcppDeepState_int();
  a_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesSUR/inst/testfiles/randIntUniform/inputs/a");
  a_stream << a;
  std::cout << "a values: "<< a << std::endl;
  a_stream.close();
  int b  = RcppDeepState_int();
  b_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesSUR/inst/testfiles/randIntUniform/inputs/b");
  b_stream << b;
  std::cout << "b values: "<< b << std::endl;
  b_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    randIntUniform(a,b);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
