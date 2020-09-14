#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix start_tree(double start_mean, double start_sd);

TEST(bartBMA_deepstate_test,start_tree_test){
  std::ofstream start_mean_stream;
  std::ofstream start_sd_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double start_mean  = RcppDeepState_double();
  start_mean_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/start_tree/inputs/start_mean");
  start_mean_stream << start_mean;
  std::cout << "start_mean values: "<< start_mean << std::endl;
  start_mean_stream.close();
  double start_sd  = RcppDeepState_double();
  start_sd_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/start_tree/inputs/start_sd");
  start_sd_stream << start_sd;
  std::cout << "start_sd values: "<< start_sd << std::endl;
  start_sd_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    start_tree(start_mean,start_sd);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
