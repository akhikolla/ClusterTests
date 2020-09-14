#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double l_const_uniwnorm(double k);

TEST(BAMBI_deepstate_test,l_const_uniwnorm_test){
  std::ofstream k_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double k  = RcppDeepState_double();
  k_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/l_const_uniwnorm/inputs/k");
  k_stream << k;
  std::cout << "k values: "<< k << std::endl;
  k_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    l_const_uniwnorm(k);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
