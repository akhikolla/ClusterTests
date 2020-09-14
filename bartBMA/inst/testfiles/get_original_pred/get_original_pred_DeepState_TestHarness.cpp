#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector get_original_pred(double low, double high, double sp_low, double sp_high, NumericVector sum_preds);

TEST(bartBMA_deepstate_test,get_original_pred_test){
  std::ofstream low_stream;
  std::ofstream high_stream;
  std::ofstream sp_low_stream;
  std::ofstream sp_high_stream;
  std::ofstream sum_preds_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double low  = RcppDeepState_double();
  low_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_original_pred/inputs/low");
  low_stream << low;
  std::cout << "low values: "<< low << std::endl;
  low_stream.close();
  double high  = RcppDeepState_double();
  high_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_original_pred/inputs/high");
  high_stream << high;
  std::cout << "high values: "<< high << std::endl;
  high_stream.close();
  double sp_low  = RcppDeepState_double();
  sp_low_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_original_pred/inputs/sp_low");
  sp_low_stream << sp_low;
  std::cout << "sp_low values: "<< sp_low << std::endl;
  sp_low_stream.close();
  double sp_high  = RcppDeepState_double();
  sp_high_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_original_pred/inputs/sp_high");
  sp_high_stream << sp_high;
  std::cout << "sp_high values: "<< sp_high << std::endl;
  sp_high_stream.close();
  NumericVector sum_preds  = RcppDeepState_NumericVector();
  sum_preds_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_original_pred/inputs/sum_preds");
  sum_preds_stream << sum_preds;
  std::cout << "sum_preds values: "<< sum_preds << std::endl;
  sum_preds_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    get_original_pred(low,high,sp_low,sp_high,sum_preds);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
