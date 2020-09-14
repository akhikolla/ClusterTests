#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix find_obs_to_update_grow(NumericMatrix prior_tree_matrix_temp, double left_daughter, int d, NumericVector ld_obs, NumericVector rd_obs);

TEST(bartBMA_deepstate_test,find_obs_to_update_grow_test){
  std::ofstream prior_tree_matrix_temp_stream;
  std::ofstream left_daughter_stream;
  std::ofstream d_stream;
  std::ofstream ld_obs_stream;
  std::ofstream rd_obs_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix prior_tree_matrix_temp  = RcppDeepState_NumericMatrix();
  prior_tree_matrix_temp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/find_obs_to_update_grow/inputs/prior_tree_matrix_temp");
  prior_tree_matrix_temp_stream << prior_tree_matrix_temp;
  std::cout << "prior_tree_matrix_temp values: "<< prior_tree_matrix_temp << std::endl;
  prior_tree_matrix_temp_stream.close();
  double left_daughter  = RcppDeepState_double();
  left_daughter_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/find_obs_to_update_grow/inputs/left_daughter");
  left_daughter_stream << left_daughter;
  std::cout << "left_daughter values: "<< left_daughter << std::endl;
  left_daughter_stream.close();
  int d  = RcppDeepState_int();
  d_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/find_obs_to_update_grow/inputs/d");
  d_stream << d;
  std::cout << "d values: "<< d << std::endl;
  d_stream.close();
  NumericVector ld_obs  = RcppDeepState_NumericVector();
  ld_obs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/find_obs_to_update_grow/inputs/ld_obs");
  ld_obs_stream << ld_obs;
  std::cout << "ld_obs values: "<< ld_obs << std::endl;
  ld_obs_stream.close();
  NumericVector rd_obs  = RcppDeepState_NumericVector();
  rd_obs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/find_obs_to_update_grow/inputs/rd_obs");
  rd_obs_stream << rd_obs;
  std::cout << "rd_obs values: "<< rd_obs << std::endl;
  rd_obs_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    find_obs_to_update_grow(prior_tree_matrix_temp,left_daughter,d,ld_obs,rd_obs);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
