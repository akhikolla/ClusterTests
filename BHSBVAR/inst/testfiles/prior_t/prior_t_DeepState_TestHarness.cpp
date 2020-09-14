#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double prior_t(const double a1, const double p1, const double sigma1, const double nu);

TEST(BHSBVAR_deepstate_test,prior_t_test){
  std::ofstream a1_stream;
  std::ofstream p1_stream;
  std::ofstream sigma1_stream;
  std::ofstream nu_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double a1  = RcppDeepState_double();
  a1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BHSBVAR/inst/testfiles/prior_t/inputs/a1");
  a1_stream << a1;
  std::cout << "a1 values: "<< a1 << std::endl;
  a1_stream.close();
  double p1  = RcppDeepState_double();
  p1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BHSBVAR/inst/testfiles/prior_t/inputs/p1");
  p1_stream << p1;
  std::cout << "p1 values: "<< p1 << std::endl;
  p1_stream.close();
  double sigma1  = RcppDeepState_double();
  sigma1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BHSBVAR/inst/testfiles/prior_t/inputs/sigma1");
  sigma1_stream << sigma1;
  std::cout << "sigma1 values: "<< sigma1 << std::endl;
  sigma1_stream.close();
  double nu  = RcppDeepState_double();
  nu_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BHSBVAR/inst/testfiles/prior_t/inputs/nu");
  nu_stream << nu;
  std::cout << "nu values: "<< nu << std::endl;
  nu_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    prior_t(a1,p1,sigma1,nu);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
