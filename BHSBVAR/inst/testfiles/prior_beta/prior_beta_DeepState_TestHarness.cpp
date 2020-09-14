#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double prior_beta(const double a1, const double sh1, const double sh2);

TEST(BHSBVAR_deepstate_test,prior_beta_test){
  std::ofstream a1_stream;
  std::ofstream sh1_stream;
  std::ofstream sh2_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double a1  = RcppDeepState_double();
  a1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BHSBVAR/inst/testfiles/prior_beta/inputs/a1");
  a1_stream << a1;
  std::cout << "a1 values: "<< a1 << std::endl;
  a1_stream.close();
  double sh1  = RcppDeepState_double();
  sh1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BHSBVAR/inst/testfiles/prior_beta/inputs/sh1");
  sh1_stream << sh1;
  std::cout << "sh1 values: "<< sh1 << std::endl;
  sh1_stream.close();
  double sh2  = RcppDeepState_double();
  sh2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BHSBVAR/inst/testfiles/prior_beta/inputs/sh2");
  sh2_stream << sh2;
  std::cout << "sh2 values: "<< sh2 << std::endl;
  sh2_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    prior_beta(a1,sh1,sh2);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
