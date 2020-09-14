#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector find_term_nodes_gs(NumericMatrix tree_table);

TEST(bartBMA_deepstate_test,find_term_nodes_gs_test){
  std::ofstream tree_table_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix tree_table  = RcppDeepState_NumericMatrix();
  tree_table_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/find_term_nodes_gs/inputs/tree_table");
  tree_table_stream << tree_table;
  std::cout << "tree_table values: "<< tree_table << std::endl;
  tree_table_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    find_term_nodes_gs(tree_table);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
