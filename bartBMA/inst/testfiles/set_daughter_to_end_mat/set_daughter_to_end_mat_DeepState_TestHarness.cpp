#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix set_daughter_to_end_mat(int d, NumericMatrix prior_tree_matrix_temp, double left_daughter, NumericVector ld_obs, NumericVector rd_obs);

TEST(bartBMA_deepstate_test,set_daughter_to_end_mat_test){
  std::ofstream d_stream;
  std::ofstream prior_tree_matrix_temp_stream;
  std::ofstream left_daughter_stream;
  std::ofstream ld_obs_stream;
  std::ofstream rd_obs_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int d  = RcppDeepState_int();
  d_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_daughter_to_end_mat/inputs/d");
  d_stream << d;
  std::cout << "d values: "<< d << std::endl;
  d_stream.close();
  NumericMatrix prior_tree_matrix_temp  = RcppDeepState_NumericMatrix();
  prior_tree_matrix_temp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_daughter_to_end_mat/inputs/prior_tree_matrix_temp");
  prior_tree_matrix_temp_stream << prior_tree_matrix_temp;
  std::cout << "prior_tree_matrix_temp values: "<< prior_tree_matrix_temp << std::endl;
  prior_tree_matrix_temp_stream.close();
  double left_daughter  = RcppDeepState_double();
  left_daughter_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_daughter_to_end_mat/inputs/left_daughter");
  left_daughter_stream << left_daughter;
  std::cout << "left_daughter values: "<< left_daughter << std::endl;
  left_daughter_stream.close();
  NumericVector ld_obs  = RcppDeepState_NumericVector();
  ld_obs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_daughter_to_end_mat/inputs/ld_obs");
  ld_obs_stream << ld_obs;
  std::cout << "ld_obs values: "<< ld_obs << std::endl;
  ld_obs_stream.close();
  NumericVector rd_obs  = RcppDeepState_NumericVector();
  rd_obs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_daughter_to_end_mat/inputs/rd_obs");
  rd_obs_stream << rd_obs;
  std::cout << "rd_obs values: "<< rd_obs << std::endl;
  rd_obs_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    set_daughter_to_end_mat(d,prior_tree_matrix_temp,left_daughter,ld_obs,rd_obs);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
