#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::vec dwnorm2_manyx_onepar(arma::mat x, double k1, double k2, double k3, double mu1, double mu2, arma::mat omega_2pi);

TEST(BAMBI_deepstate_test,dwnorm2_manyx_onepar_test){
  std::ofstream x_stream;
  std::ofstream k1_stream;
  std::ofstream k2_stream;
  std::ofstream k3_stream;
  std::ofstream mu1_stream;
  std::ofstream mu2_stream;
  std::ofstream omega_2pi_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat x  = RcppDeepState_mat();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/dwnorm2_manyx_onepar/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  double k1  = RcppDeepState_double();
  k1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/dwnorm2_manyx_onepar/inputs/k1");
  k1_stream << k1;
  std::cout << "k1 values: "<< k1 << std::endl;
  k1_stream.close();
  double k2  = RcppDeepState_double();
  k2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/dwnorm2_manyx_onepar/inputs/k2");
  k2_stream << k2;
  std::cout << "k2 values: "<< k2 << std::endl;
  k2_stream.close();
  double k3  = RcppDeepState_double();
  k3_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/dwnorm2_manyx_onepar/inputs/k3");
  k3_stream << k3;
  std::cout << "k3 values: "<< k3 << std::endl;
  k3_stream.close();
  double mu1  = RcppDeepState_double();
  mu1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/dwnorm2_manyx_onepar/inputs/mu1");
  mu1_stream << mu1;
  std::cout << "mu1 values: "<< mu1 << std::endl;
  mu1_stream.close();
  double mu2  = RcppDeepState_double();
  mu2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/dwnorm2_manyx_onepar/inputs/mu2");
  mu2_stream << mu2;
  std::cout << "mu2 values: "<< mu2 << std::endl;
  mu2_stream.close();
  arma::mat omega_2pi  = RcppDeepState_mat();
  omega_2pi_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/dwnorm2_manyx_onepar/inputs/omega_2pi");
  omega_2pi_stream << omega_2pi;
  std::cout << "omega_2pi values: "<< omega_2pi << std::endl;
  omega_2pi_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    dwnorm2_manyx_onepar(x,k1,k2,k3,mu1,mu2,omega_2pi);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
