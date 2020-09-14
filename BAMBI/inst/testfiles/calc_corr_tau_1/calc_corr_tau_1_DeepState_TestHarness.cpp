#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double calc_corr_tau_1(arma::mat samp_mat);

TEST(BAMBI_deepstate_test,calc_corr_tau_1_test){
  std::ofstream samp_mat_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat samp_mat  = RcppDeepState_mat();
  samp_mat_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/calc_corr_tau_1/inputs/samp_mat");
  samp_mat_stream << samp_mat;
  std::cout << "samp_mat values: "<< samp_mat << std::endl;
  samp_mat_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    calc_corr_tau_1(samp_mat);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
