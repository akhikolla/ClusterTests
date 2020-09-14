#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

vec cgetC(double e, int k);

TEST(bayesm_deepstate_test,cgetC_test){
  std::ofstream e_stream;
  std::ofstream k_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double e  = RcppDeepState_double();
  e_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bayesm/inst/testfiles/cgetC/inputs/e");
  e_stream << e;
  std::cout << "e values: "<< e << std::endl;
  e_stream.close();
  int k  = RcppDeepState_int();
  k_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bayesm/inst/testfiles/cgetC/inputs/k");
  k_stream << k;
  std::cout << "k values: "<< k << std::endl;
  k_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    cgetC(e,k);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
