#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector get_imp_vars(NumericVector split_vars, int num_col, NumericVector current_vars);

TEST(bartBMA_deepstate_test,get_imp_vars_test){
  std::ofstream split_vars_stream;
  std::ofstream num_col_stream;
  std::ofstream current_vars_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector split_vars  = RcppDeepState_NumericVector();
  split_vars_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_imp_vars/inputs/split_vars");
  split_vars_stream << split_vars;
  std::cout << "split_vars values: "<< split_vars << std::endl;
  split_vars_stream.close();
  int num_col  = RcppDeepState_int();
  num_col_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_imp_vars/inputs/num_col");
  num_col_stream << num_col;
  std::cout << "num_col values: "<< num_col << std::endl;
  num_col_stream.close();
  NumericVector current_vars  = RcppDeepState_NumericVector();
  current_vars_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_imp_vars/inputs/current_vars");
  current_vars_stream << current_vars;
  std::cout << "current_vars values: "<< current_vars << std::endl;
  current_vars_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    get_imp_vars(split_vars,num_col,current_vars);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
