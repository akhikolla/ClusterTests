#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector takefirst(NumericVector xs, int pos);

TEST(BiDAG_deepstate_test,takefirst_test){
  std::ofstream xs_stream;
  std::ofstream pos_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector xs  = RcppDeepState_NumericVector();
  xs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BiDAG/inst/testfiles/takefirst/inputs/xs");
  xs_stream << xs;
  std::cout << "xs values: "<< xs << std::endl;
  xs_stream.close();
  int pos  = RcppDeepState_int();
  pos_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BiDAG/inst/testfiles/takefirst/inputs/pos");
  pos_stream << pos;
  std::cout << "pos values: "<< pos << std::endl;
  pos_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    takefirst(xs,pos);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
