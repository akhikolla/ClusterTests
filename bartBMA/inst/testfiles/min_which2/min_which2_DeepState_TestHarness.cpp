#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

List min_which2(NumericVector array, int n, double minout, int whichout);

TEST(bartBMA_deepstate_test,min_which2_test){
  std::ofstream array_stream;
  std::ofstream n_stream;
  std::ofstream minout_stream;
  std::ofstream whichout_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector array  = RcppDeepState_NumericVector();
  array_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/min_which2/inputs/array");
  array_stream << array;
  std::cout << "array values: "<< array << std::endl;
  array_stream.close();
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/min_which2/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  double minout  = RcppDeepState_double();
  minout_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/min_which2/inputs/minout");
  minout_stream << minout;
  std::cout << "minout values: "<< minout << std::endl;
  minout_stream.close();
  int whichout  = RcppDeepState_int();
  whichout_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/min_which2/inputs/whichout");
  whichout_stream << whichout;
  std::cout << "whichout values: "<< whichout << std::endl;
  whichout_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    min_which2(array,n,minout,whichout);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
