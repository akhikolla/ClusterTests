#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector rcpp_row_quantile(NumericMatrix data, const double q);

TEST(biwavelet_deepstate_test,rcpp_row_quantile_test){
  std::ofstream data_stream;
  std::ofstream q_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix data  = RcppDeepState_NumericMatrix();
  data_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/biwavelet/inst/testfiles/rcpp_row_quantile/inputs/data");
  data_stream << data;
  std::cout << "data values: "<< data << std::endl;
  data_stream.close();
  double q  = RcppDeepState_double();
  q_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/biwavelet/inst/testfiles/rcpp_row_quantile/inputs/q");
  q_stream << q;
  std::cout << "q values: "<< q << std::endl;
  q_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    rcpp_row_quantile(data,q);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
