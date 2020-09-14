#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double det_chol_downdate(const NumericMatrix L, const NumericVector v);

TEST(Benchmarking_deepstate_test,det_chol_downdate_test){
  std::ofstream L_stream;
  std::ofstream v_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix L  = RcppDeepState_NumericMatrix();
  L_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/det_chol_downdate/inputs/L");
  L_stream << L;
  std::cout << "L values: "<< L << std::endl;
  L_stream.close();
  NumericVector v  = RcppDeepState_NumericVector();
  v_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/det_chol_downdate/inputs/v");
  v_stream << v;
  std::cout << "v values: "<< v << std::endl;
  v_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    det_chol_downdate(L,v);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
