#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector takelast(NumericVector xs, int pos, int n);

TEST(BiDAG_deepstate_test,takelast_test){
  std::ofstream xs_stream;
  std::ofstream pos_stream;
  std::ofstream n_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector xs  = RcppDeepState_NumericVector();
  xs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BiDAG/inst/testfiles/takelast/inputs/xs");
  xs_stream << xs;
  std::cout << "xs values: "<< xs << std::endl;
  xs_stream.close();
  int pos  = RcppDeepState_int();
  pos_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BiDAG/inst/testfiles/takelast/inputs/pos");
  pos_stream << pos;
  std::cout << "pos values: "<< pos << std::endl;
  pos_stream.close();
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BiDAG/inst/testfiles/takelast/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    takelast(xs,pos,n);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
