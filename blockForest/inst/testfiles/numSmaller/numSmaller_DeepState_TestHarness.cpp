#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::IntegerVector numSmaller(Rcpp::NumericVector values, Rcpp::NumericVector reference);

TEST(blockForest_deepstate_test,numSmaller_test){
  std::ofstream values_stream;
  std::ofstream reference_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector values  = RcppDeepState_NumericVector();
  values_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/blockForest/inst/testfiles/numSmaller/inputs/values");
  values_stream << values;
  std::cout << "values values: "<< values << std::endl;
  values_stream.close();
  NumericVector reference  = RcppDeepState_NumericVector();
  reference_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/blockForest/inst/testfiles/numSmaller/inputs/reference");
  reference_stream << reference;
  std::cout << "reference values: "<< reference << std::endl;
  reference_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    numSmaller(values,reference);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
