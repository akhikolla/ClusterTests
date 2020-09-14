#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

List vmsin_var_corr_mc(double k1, double k2, double k3, arma::mat uni_rand, int ncores);

TEST(BAMBI_deepstate_test,vmsin_var_corr_mc_test){
  std::ofstream k1_stream;
  std::ofstream k2_stream;
  std::ofstream k3_stream;
  std::ofstream uni_rand_stream;
  std::ofstream ncores_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double k1  = RcppDeepState_double();
  k1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/vmsin_var_corr_mc/inputs/k1");
  k1_stream << k1;
  std::cout << "k1 values: "<< k1 << std::endl;
  k1_stream.close();
  double k2  = RcppDeepState_double();
  k2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/vmsin_var_corr_mc/inputs/k2");
  k2_stream << k2;
  std::cout << "k2 values: "<< k2 << std::endl;
  k2_stream.close();
  double k3  = RcppDeepState_double();
  k3_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/vmsin_var_corr_mc/inputs/k3");
  k3_stream << k3;
  std::cout << "k3 values: "<< k3 << std::endl;
  k3_stream.close();
  arma::mat uni_rand  = RcppDeepState_mat();
  uni_rand_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/vmsin_var_corr_mc/inputs/uni_rand");
  uni_rand_stream << uni_rand;
  std::cout << "uni_rand values: "<< uni_rand << std::endl;
  uni_rand_stream.close();
  int ncores  = RcppDeepState_int();
  ncores_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/vmsin_var_corr_mc/inputs/ncores");
  ncores_stream << ncores;
  std::cout << "ncores values: "<< ncores << std::endl;
  ncores_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    vmsin_var_corr_mc(k1,k2,k3,uni_rand,ncores);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
