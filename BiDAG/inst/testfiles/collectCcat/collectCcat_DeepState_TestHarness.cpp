#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

IntegerMatrix collectCcat(IntegerVector xs, IntegerVector ys, int n, int m);

TEST(BiDAG_deepstate_test,collectCcat_test){
  std::ofstream xs_stream;
  std::ofstream ys_stream;
  std::ofstream n_stream;
  std::ofstream m_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  IntegerVector xs  = RcppDeepState_IntegerVector();
  xs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BiDAG/inst/testfiles/collectCcat/inputs/xs");
  xs_stream << xs;
  std::cout << "xs values: "<< xs << std::endl;
  xs_stream.close();
  IntegerVector ys  = RcppDeepState_IntegerVector();
  ys_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BiDAG/inst/testfiles/collectCcat/inputs/ys");
  ys_stream << ys;
  std::cout << "ys values: "<< ys << std::endl;
  ys_stream.close();
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BiDAG/inst/testfiles/collectCcat/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  int m  = RcppDeepState_int();
  m_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BiDAG/inst/testfiles/collectCcat/inputs/m");
  m_stream << m;
  std::cout << "m values: "<< m << std::endl;
  m_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    collectCcat(xs,ys,n,m);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
