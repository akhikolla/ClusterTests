#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix add_rows(NumericMatrix prior_tree_table_temp, int grow_node);

TEST(bartBMA_deepstate_test,add_rows_test){
  std::ofstream prior_tree_table_temp_stream;
  std::ofstream grow_node_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix prior_tree_table_temp  = RcppDeepState_NumericMatrix();
  prior_tree_table_temp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/add_rows/inputs/prior_tree_table_temp");
  prior_tree_table_temp_stream << prior_tree_table_temp;
  std::cout << "prior_tree_table_temp values: "<< prior_tree_table_temp << std::endl;
  prior_tree_table_temp_stream.close();
  int grow_node  = RcppDeepState_int();
  grow_node_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/add_rows/inputs/grow_node");
  grow_node_stream << grow_node;
  std::cout << "grow_node values: "<< grow_node << std::endl;
  grow_node_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    add_rows(prior_tree_table_temp,grow_node);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
