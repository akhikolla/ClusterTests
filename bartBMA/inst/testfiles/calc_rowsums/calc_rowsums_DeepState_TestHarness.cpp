#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector calc_rowsums(NumericMatrix predictions);

TEST(bartBMA_deepstate_test,calc_rowsums_test){
  std::ofstream predictions_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix predictions  = RcppDeepState_NumericMatrix();
  predictions_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/calc_rowsums/inputs/predictions");
  predictions_stream << predictions;
  std::cout << "predictions values: "<< predictions << std::endl;
  predictions_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    calc_rowsums(predictions);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
