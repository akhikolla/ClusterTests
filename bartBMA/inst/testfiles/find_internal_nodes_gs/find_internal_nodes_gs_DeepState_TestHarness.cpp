#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector find_internal_nodes_gs(NumericMatrix treetable);

TEST(bartBMA_deepstate_test,find_internal_nodes_gs_test){
  std::ofstream treetable_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix treetable  = RcppDeepState_NumericMatrix();
  treetable_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/find_internal_nodes_gs/inputs/treetable");
  treetable_stream << treetable;
  std::cout << "treetable values: "<< treetable << std::endl;
  treetable_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    find_internal_nodes_gs(treetable);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
