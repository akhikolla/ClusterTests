#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector update_mean_var(NumericMatrix tree_table, NumericMatrix tree_matrix, NumericVector resids, double a);

TEST(bartBMA_deepstate_test,update_mean_var_test){
  std::ofstream tree_table_stream;
  std::ofstream tree_matrix_stream;
  std::ofstream resids_stream;
  std::ofstream a_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix tree_table  = RcppDeepState_NumericMatrix();
  tree_table_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/update_mean_var/inputs/tree_table");
  tree_table_stream << tree_table;
  std::cout << "tree_table values: "<< tree_table << std::endl;
  tree_table_stream.close();
  NumericMatrix tree_matrix  = RcppDeepState_NumericMatrix();
  tree_matrix_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/update_mean_var/inputs/tree_matrix");
  tree_matrix_stream << tree_matrix;
  std::cout << "tree_matrix values: "<< tree_matrix << std::endl;
  tree_matrix_stream.close();
  NumericVector resids  = RcppDeepState_NumericVector();
  resids_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/update_mean_var/inputs/resids");
  resids_stream << resids;
  std::cout << "resids values: "<< resids << std::endl;
  resids_stream.close();
  double a  = RcppDeepState_double();
  a_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/update_mean_var/inputs/a");
  a_stream << a;
  std::cout << "a values: "<< a << std::endl;
  a_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    update_mean_var(tree_table,tree_matrix,resids,a);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
