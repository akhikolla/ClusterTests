#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector extractT(IntegerVector xs, IntegerVector ys, NumericMatrix ts);

TEST(BiDAG_deepstate_test,extractT_test){
  std::ofstream xs_stream;
  std::ofstream ys_stream;
  std::ofstream ts_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  IntegerVector xs  = RcppDeepState_IntegerVector();
  xs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BiDAG/inst/testfiles/extractT/inputs/xs");
  xs_stream << xs;
  std::cout << "xs values: "<< xs << std::endl;
  xs_stream.close();
  IntegerVector ys  = RcppDeepState_IntegerVector();
  ys_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BiDAG/inst/testfiles/extractT/inputs/ys");
  ys_stream << ys;
  std::cout << "ys values: "<< ys << std::endl;
  ys_stream.close();
  NumericMatrix ts  = RcppDeepState_NumericMatrix();
  ts_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BiDAG/inst/testfiles/extractT/inputs/ts");
  ts_stream << ts;
  std::cout << "ts values: "<< ts << std::endl;
  ts_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    extractT(xs,ys,ts);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
