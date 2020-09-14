#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::vec log_const_vmsin_all(arma::mat par_mat);

TEST(BAMBI_deepstate_test,log_const_vmsin_all_test){
  std::ofstream par_mat_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat par_mat  = RcppDeepState_mat();
  par_mat_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/log_const_vmsin_all/inputs/par_mat");
  par_mat_stream << par_mat;
  std::cout << "par_mat values: "<< par_mat << std::endl;
  par_mat_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    log_const_vmsin_all(par_mat);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
