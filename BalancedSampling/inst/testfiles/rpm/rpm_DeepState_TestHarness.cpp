#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector rpm(NumericVector prob);

TEST(BalancedSampling_deepstate_test,rpm_test){
  std::ofstream prob_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector prob  = RcppDeepState_NumericVector();
  prob_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/rpm/inputs/prob");
  prob_stream << prob;
  std::cout << "prob values: "<< prob << std::endl;
  prob_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    rpm(prob);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
