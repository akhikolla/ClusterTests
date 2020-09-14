#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double const_vmsin(double k1, double k2, double lambda);

TEST(BAMBI_deepstate_test,const_vmsin_test){
  std::ofstream k1_stream;
  std::ofstream k2_stream;
  std::ofstream lambda_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double k1  = RcppDeepState_double();
  k1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/const_vmsin/inputs/k1");
  k1_stream << k1;
  std::cout << "k1 values: "<< k1 << std::endl;
  k1_stream.close();
  double k2  = RcppDeepState_double();
  k2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/const_vmsin/inputs/k2");
  k2_stream << k2;
  std::cout << "k2 values: "<< k2 << std::endl;
  k2_stream.close();
  double lambda  = RcppDeepState_double();
  lambda_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/const_vmsin/inputs/lambda");
  lambda_stream << lambda;
  std::cout << "lambda values: "<< lambda << std::endl;
  lambda_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    const_vmsin(k1,k2,lambda);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
