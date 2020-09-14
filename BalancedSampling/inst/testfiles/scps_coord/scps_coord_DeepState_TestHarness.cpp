#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

IntegerVector scps_coord(NumericVector prob, NumericMatrix x, NumericVector rand);

TEST(BalancedSampling_deepstate_test,scps_coord_test){
  std::ofstream prob_stream;
  std::ofstream x_stream;
  std::ofstream rand_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector prob  = RcppDeepState_NumericVector();
  prob_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/scps_coord/inputs/prob");
  prob_stream << prob;
  std::cout << "prob values: "<< prob << std::endl;
  prob_stream.close();
  NumericMatrix x  = RcppDeepState_NumericMatrix();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/scps_coord/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  NumericVector rand  = RcppDeepState_NumericVector();
  rand_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/scps_coord/inputs/rand");
  rand_stream << rand;
  std::cout << "rand values: "<< rand << std::endl;
  rand_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    scps_coord(prob,x,rand);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
