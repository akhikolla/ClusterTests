#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::mat updateC(arma::mat Z, double theta, double tau);

TEST(BAREB_deepstate_test,updateC_test){
  std::ofstream Z_stream;
  std::ofstream theta_stream;
  std::ofstream tau_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat Z  = RcppDeepState_mat();
  Z_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updateC/inputs/Z");
  Z_stream << Z;
  std::cout << "Z values: "<< Z << std::endl;
  Z_stream.close();
  double theta  = RcppDeepState_double();
  theta_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updateC/inputs/theta");
  theta_stream << theta;
  std::cout << "theta values: "<< theta << std::endl;
  theta_stream.close();
  double tau  = RcppDeepState_double();
  tau_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updateC/inputs/tau");
  tau_stream << tau;
  std::cout << "tau values: "<< tau << std::endl;
  tau_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    updateC(Z,theta,tau);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
