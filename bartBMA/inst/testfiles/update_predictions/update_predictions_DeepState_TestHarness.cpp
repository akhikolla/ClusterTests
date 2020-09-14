#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

List update_predictions(NumericMatrix tree_table, NumericMatrix tree_matrix, NumericVector new_mean, int n);

TEST(bartBMA_deepstate_test,update_predictions_test){
  std::ofstream tree_table_stream;
  std::ofstream tree_matrix_stream;
  std::ofstream new_mean_stream;
  std::ofstream n_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix tree_table  = RcppDeepState_NumericMatrix();
  tree_table_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/update_predictions/inputs/tree_table");
  tree_table_stream << tree_table;
  std::cout << "tree_table values: "<< tree_table << std::endl;
  tree_table_stream.close();
  NumericMatrix tree_matrix  = RcppDeepState_NumericMatrix();
  tree_matrix_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/update_predictions/inputs/tree_matrix");
  tree_matrix_stream << tree_matrix;
  std::cout << "tree_matrix values: "<< tree_matrix << std::endl;
  tree_matrix_stream.close();
  NumericVector new_mean  = RcppDeepState_NumericVector();
  new_mean_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/update_predictions/inputs/new_mean");
  new_mean_stream << new_mean;
  std::cout << "new_mean values: "<< new_mean << std::endl;
  new_mean_stream.close();
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/update_predictions/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    update_predictions(tree_table,tree_matrix,new_mean,n);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
