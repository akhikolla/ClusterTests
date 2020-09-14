#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double do_rgig1(double lambda, double chi, double psi);

TEST(BGVAR_deepstate_test,do_rgig1_test){
  std::ofstream lambda_stream;
  std::ofstream chi_stream;
  std::ofstream psi_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double lambda  = RcppDeepState_double();
  lambda_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BGVAR/inst/testfiles/do_rgig1/inputs/lambda");
  lambda_stream << lambda;
  std::cout << "lambda values: "<< lambda << std::endl;
  lambda_stream.close();
  double chi  = RcppDeepState_double();
  chi_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BGVAR/inst/testfiles/do_rgig1/inputs/chi");
  chi_stream << chi;
  std::cout << "chi values: "<< chi << std::endl;
  chi_stream.close();
  double psi  = RcppDeepState_double();
  psi_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BGVAR/inst/testfiles/do_rgig1/inputs/psi");
  psi_stream << psi;
  std::cout << "psi values: "<< psi << std::endl;
  psi_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    do_rgig1(lambda,chi,psi);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
