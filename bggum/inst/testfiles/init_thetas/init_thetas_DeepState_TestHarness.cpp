#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector init_thetas(int n, double mean, double sd);

TEST(bggum_deepstate_test,init_thetas_test){
  std::ofstream n_stream;
  std::ofstream mean_stream;
  std::ofstream sd_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bggum/inst/testfiles/init_thetas/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  double mean  = RcppDeepState_double();
  mean_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bggum/inst/testfiles/init_thetas/inputs/mean");
  mean_stream << mean;
  std::cout << "mean values: "<< mean << std::endl;
  mean_stream.close();
  double sd  = RcppDeepState_double();
  sd_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bggum/inst/testfiles/init_thetas/inputs/sd");
  sd_stream << sd;
  std::cout << "sd values: "<< sd << std::endl;
  sd_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    init_thetas(n,mean,sd);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
