#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector calculate_resids(NumericMatrix predictions, NumericVector response);

TEST(bartBMA_deepstate_test,calculate_resids_test){
  std::ofstream predictions_stream;
  std::ofstream response_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix predictions  = RcppDeepState_NumericMatrix();
  predictions_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/calculate_resids/inputs/predictions");
  predictions_stream << predictions;
  std::cout << "predictions values: "<< predictions << std::endl;
  predictions_stream.close();
  NumericVector response  = RcppDeepState_NumericVector();
  response_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/calculate_resids/inputs/response");
  response_stream << response;
  std::cout << "response values: "<< response << std::endl;
  response_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    calculate_resids(predictions,response);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
