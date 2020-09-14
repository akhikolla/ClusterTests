#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix set_tree_to_middle(NumericVector node_to_update, NumericMatrix prior_tree_table_temp, int grow_node, double left_daughter);

TEST(bartBMA_deepstate_test,set_tree_to_middle_test){
  std::ofstream node_to_update_stream;
  std::ofstream prior_tree_table_temp_stream;
  std::ofstream grow_node_stream;
  std::ofstream left_daughter_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector node_to_update  = RcppDeepState_NumericVector();
  node_to_update_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_tree_to_middle/inputs/node_to_update");
  node_to_update_stream << node_to_update;
  std::cout << "node_to_update values: "<< node_to_update << std::endl;
  node_to_update_stream.close();
  NumericMatrix prior_tree_table_temp  = RcppDeepState_NumericMatrix();
  prior_tree_table_temp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_tree_to_middle/inputs/prior_tree_table_temp");
  prior_tree_table_temp_stream << prior_tree_table_temp;
  std::cout << "prior_tree_table_temp values: "<< prior_tree_table_temp << std::endl;
  prior_tree_table_temp_stream.close();
  int grow_node  = RcppDeepState_int();
  grow_node_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_tree_to_middle/inputs/grow_node");
  grow_node_stream << grow_node;
  std::cout << "grow_node values: "<< grow_node << std::endl;
  grow_node_stream.close();
  double left_daughter  = RcppDeepState_double();
  left_daughter_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/set_tree_to_middle/inputs/left_daughter");
  left_daughter_stream << left_daughter;
  std::cout << "left_daughter values: "<< left_daughter << std::endl;
  left_daughter_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    set_tree_to_middle(node_to_update,prior_tree_table_temp,grow_node,left_daughter);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
