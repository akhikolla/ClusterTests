#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double sldmvnorm(arma::mat z_t, arma::mat Sigma);

TEST(beyondWhittle_deepstate_test,sldmvnorm_test){
  std::ofstream z_t_stream;
  std::ofstream Sigma_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat z_t  = RcppDeepState_mat();
  z_t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/sldmvnorm/inputs/z_t");
  z_t_stream << z_t;
  std::cout << "z_t values: "<< z_t << std::endl;
  z_t_stream.close();
  arma::mat Sigma  = RcppDeepState_mat();
  Sigma_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/sldmvnorm/inputs/Sigma");
  Sigma_stream << Sigma;
  std::cout << "Sigma values: "<< Sigma << std::endl;
  Sigma_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    sldmvnorm(z_t,Sigma);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
