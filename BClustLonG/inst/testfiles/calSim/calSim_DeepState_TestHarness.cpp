#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::mat calSim(arma::mat mat);

TEST(BClustLonG_deepstate_test,calSim_test){
  std::ofstream mat_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat mat  = RcppDeepState_mat();
  mat_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BClustLonG/inst/testfiles/calSim/inputs/mat");
  mat_stream << mat;
  std::cout << "mat values: "<< mat << std::endl;
  mat_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    calSim(mat);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
