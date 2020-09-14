#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector collectC(IntegerVector xs, NumericVector ys, int n);

TEST(Bestie_deepstate_test,collectC_test){
  std::ofstream xs_stream;
  std::ofstream ys_stream;
  std::ofstream n_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  IntegerVector xs  = RcppDeepState_IntegerVector();
  xs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Bestie/inst/testfiles/collectC/inputs/xs");
  xs_stream << xs;
  std::cout << "xs values: "<< xs << std::endl;
  xs_stream.close();
  NumericVector ys  = RcppDeepState_NumericVector();
  ys_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Bestie/inst/testfiles/collectC/inputs/ys");
  ys_stream << ys;
  std::cout << "ys values: "<< ys << std::endl;
  ys_stream.close();
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Bestie/inst/testfiles/collectC/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    collectC(xs,ys,n);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
