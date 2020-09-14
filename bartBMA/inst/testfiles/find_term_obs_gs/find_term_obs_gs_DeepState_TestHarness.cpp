#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector find_term_obs_gs(NumericMatrix tree_matrix_temp, double terminal_node);

TEST(bartBMA_deepstate_test,find_term_obs_gs_test){
  std::ofstream tree_matrix_temp_stream;
  std::ofstream terminal_node_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix tree_matrix_temp  = RcppDeepState_NumericMatrix();
  tree_matrix_temp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/find_term_obs_gs/inputs/tree_matrix_temp");
  tree_matrix_temp_stream << tree_matrix_temp;
  std::cout << "tree_matrix_temp values: "<< tree_matrix_temp << std::endl;
  tree_matrix_temp_stream.close();
  double terminal_node  = RcppDeepState_double();
  terminal_node_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/find_term_obs_gs/inputs/terminal_node");
  terminal_node_stream << terminal_node;
  std::cout << "terminal_node values: "<< terminal_node << std::endl;
  terminal_node_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    find_term_obs_gs(tree_matrix_temp,terminal_node);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
