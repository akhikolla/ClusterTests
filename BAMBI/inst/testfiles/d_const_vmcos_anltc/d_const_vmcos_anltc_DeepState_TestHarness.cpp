#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::vec d_const_vmcos_anltc(double k1, double k2, double k3);

TEST(BAMBI_deepstate_test,d_const_vmcos_anltc_test){
  std::ofstream k1_stream;
  std::ofstream k2_stream;
  std::ofstream k3_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double k1  = RcppDeepState_double();
  k1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/d_const_vmcos_anltc/inputs/k1");
  k1_stream << k1;
  std::cout << "k1 values: "<< k1 << std::endl;
  k1_stream.close();
  double k2  = RcppDeepState_double();
  k2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/d_const_vmcos_anltc/inputs/k2");
  k2_stream << k2;
  std::cout << "k2 values: "<< k2 << std::endl;
  k2_stream.close();
  double k3  = RcppDeepState_double();
  k3_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/d_const_vmcos_anltc/inputs/k3");
  k3_stream << k3;
  std::cout << "k3 values: "<< k3 << std::endl;
  k3_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    d_const_vmcos_anltc(k1,k2,k3);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
