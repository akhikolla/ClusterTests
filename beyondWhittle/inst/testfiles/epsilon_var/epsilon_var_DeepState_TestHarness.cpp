#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::mat epsilon_var(arma::mat zt, arma::mat ar);

TEST(beyondWhittle_deepstate_test,epsilon_var_test){
  std::ofstream zt_stream;
  std::ofstream ar_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat zt  = RcppDeepState_mat();
  zt_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/epsilon_var/inputs/zt");
  zt_stream << zt;
  std::cout << "zt values: "<< zt << std::endl;
  zt_stream.close();
  arma::mat ar  = RcppDeepState_mat();
  ar_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/epsilon_var/inputs/ar");
  ar_stream << ar;
  std::cout << "ar values: "<< ar << std::endl;
  ar_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    epsilon_var(zt,ar);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
