#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

List rcpp_wt_bases_dog(const NumericVector k, const double scale, const int param);

TEST(biwavelet_deepstate_test,rcpp_wt_bases_dog_test){
  std::ofstream k_stream;
  std::ofstream scale_stream;
  std::ofstream param_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector k  = RcppDeepState_NumericVector();
  k_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/biwavelet/inst/testfiles/rcpp_wt_bases_dog/inputs/k");
  k_stream << k;
  std::cout << "k values: "<< k << std::endl;
  k_stream.close();
  double scale  = RcppDeepState_double();
  scale_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/biwavelet/inst/testfiles/rcpp_wt_bases_dog/inputs/scale");
  scale_stream << scale;
  std::cout << "scale values: "<< scale << std::endl;
  scale_stream.close();
  int param  = RcppDeepState_int();
  param_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/biwavelet/inst/testfiles/rcpp_wt_bases_dog/inputs/param");
  param_stream << param;
  std::cout << "param values: "<< param << std::endl;
  param_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    rcpp_wt_bases_dog(k,scale,param);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
