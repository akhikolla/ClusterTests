#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::DataFrame find_min_dists_cpp(NumericMatrix mOne, NumericMatrix mTwo);

TEST(Bioi_deepstate_test,find_min_dists_cpp_test){
  std::ofstream mOne_stream;
  std::ofstream mTwo_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix mOne  = RcppDeepState_NumericMatrix();
  mOne_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Bioi/inst/testfiles/find_min_dists_cpp/inputs/mOne");
  mOne_stream << mOne;
  std::cout << "mOne values: "<< mOne << std::endl;
  mOne_stream.close();
  NumericMatrix mTwo  = RcppDeepState_NumericMatrix();
  mTwo_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Bioi/inst/testfiles/find_min_dists_cpp/inputs/mTwo");
  mTwo_stream << mTwo;
  std::cout << "mTwo values: "<< mTwo << std::endl;
  mTwo_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    find_min_dists_cpp(mOne,mTwo);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
