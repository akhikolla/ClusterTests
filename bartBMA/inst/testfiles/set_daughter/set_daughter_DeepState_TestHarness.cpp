#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix set_daughter(int left_daughter, int right_daughter, IntegerVector ld_obs, IntegerVector rd_obs, NumericMatrix tree_matrix_temp, double term_cols);

TEST(bartBMA_deepstate_test,set_daughter_test){
  std::ofstream left_daughter_stream;
  std::ofstream right_daughter_stream;
  std::ofstream ld_obs_stream;
  std::ofstream rd_obs_stream;
  std::ofstream tree_matrix_temp_stream;
  std::ofstream term_cols_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int left_daughter  = RcppDeepState_int();
  left_daughter_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_daughter/inputs/left_daughter");
  left_daughter_stream << left_daughter;
  std::cout << "left_daughter values: "<< left_daughter << std::endl;
  left_daughter_stream.close();
  int right_daughter  = RcppDeepState_int();
  right_daughter_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_daughter/inputs/right_daughter");
  right_daughter_stream << right_daughter;
  std::cout << "right_daughter values: "<< right_daughter << std::endl;
  right_daughter_stream.close();
  IntegerVector ld_obs  = RcppDeepState_IntegerVector();
  ld_obs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_daughter/inputs/ld_obs");
  ld_obs_stream << ld_obs;
  std::cout << "ld_obs values: "<< ld_obs << std::endl;
  ld_obs_stream.close();
  IntegerVector rd_obs  = RcppDeepState_IntegerVector();
  rd_obs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_daughter/inputs/rd_obs");
  rd_obs_stream << rd_obs;
  std::cout << "rd_obs values: "<< rd_obs << std::endl;
  rd_obs_stream.close();
  NumericMatrix tree_matrix_temp  = RcppDeepState_NumericMatrix();
  tree_matrix_temp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_daughter/inputs/tree_matrix_temp");
  tree_matrix_temp_stream << tree_matrix_temp;
  std::cout << "tree_matrix_temp values: "<< tree_matrix_temp << std::endl;
  tree_matrix_temp_stream.close();
  double term_cols  = RcppDeepState_double();
  term_cols_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_daughter/inputs/term_cols");
  term_cols_stream << term_cols;
  std::cout << "term_cols values: "<< term_cols << std::endl;
  term_cols_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    set_daughter(left_daughter,right_daughter,ld_obs,rd_obs,tree_matrix_temp,term_cols);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
