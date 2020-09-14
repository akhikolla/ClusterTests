#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix addcol(NumericMatrix prior_tree_matrix_temp, int grow_node, NumericVector ld_obs, NumericVector rd_obs);

TEST(bartBMA_deepstate_test,addcol_test){
  std::ofstream prior_tree_matrix_temp_stream;
  std::ofstream grow_node_stream;
  std::ofstream ld_obs_stream;
  std::ofstream rd_obs_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix prior_tree_matrix_temp  = RcppDeepState_NumericMatrix();
  prior_tree_matrix_temp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/addcol/inputs/prior_tree_matrix_temp");
  prior_tree_matrix_temp_stream << prior_tree_matrix_temp;
  std::cout << "prior_tree_matrix_temp values: "<< prior_tree_matrix_temp << std::endl;
  prior_tree_matrix_temp_stream.close();
  int grow_node  = RcppDeepState_int();
  grow_node_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/addcol/inputs/grow_node");
  grow_node_stream << grow_node;
  std::cout << "grow_node values: "<< grow_node << std::endl;
  grow_node_stream.close();
  NumericVector ld_obs  = RcppDeepState_NumericVector();
  ld_obs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/addcol/inputs/ld_obs");
  ld_obs_stream << ld_obs;
  std::cout << "ld_obs values: "<< ld_obs << std::endl;
  ld_obs_stream.close();
  NumericVector rd_obs  = RcppDeepState_NumericVector();
  rd_obs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/addcol/inputs/rd_obs");
  rd_obs_stream << rd_obs;
  std::cout << "rd_obs values: "<< rd_obs << std::endl;
  rd_obs_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    addcol(prior_tree_matrix_temp,grow_node,ld_obs,rd_obs);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
