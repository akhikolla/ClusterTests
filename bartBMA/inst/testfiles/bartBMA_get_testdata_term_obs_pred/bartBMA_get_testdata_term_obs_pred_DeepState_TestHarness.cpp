#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector bartBMA_get_testdata_term_obs_pred(NumericMatrix test_data, NumericMatrix tree_data, NumericVector term_node_means);

TEST(bartBMA_deepstate_test,bartBMA_get_testdata_term_obs_pred_test){
  std::ofstream test_data_stream;
  std::ofstream tree_data_stream;
  std::ofstream term_node_means_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix test_data  = RcppDeepState_NumericMatrix();
  test_data_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/bartBMA_get_testdata_term_obs_pred/inputs/test_data");
  test_data_stream << test_data;
  std::cout << "test_data values: "<< test_data << std::endl;
  test_data_stream.close();
  NumericMatrix tree_data  = RcppDeepState_NumericMatrix();
  tree_data_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/bartBMA_get_testdata_term_obs_pred/inputs/tree_data");
  tree_data_stream << tree_data;
  std::cout << "tree_data values: "<< tree_data << std::endl;
  tree_data_stream.close();
  NumericVector term_node_means  = RcppDeepState_NumericVector();
  term_node_means_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/bartBMA_get_testdata_term_obs_pred/inputs/term_node_means");
  term_node_means_stream << term_node_means;
  std::cout << "term_node_means values: "<< term_node_means << std::endl;
  term_node_means_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bartBMA_get_testdata_term_obs_pred(test_data,tree_data,term_node_means);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
