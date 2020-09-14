#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

DataFrame blr_pairs_cpp(NumericVector x, NumericVector y);

TEST(blorr_deepstate_test,blr_pairs_cpp_test){
  std::ofstream x_stream;
  std::ofstream y_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector x  = RcppDeepState_NumericVector();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/blorr/inst/testfiles/blr_pairs_cpp/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  NumericVector y  = RcppDeepState_NumericVector();
  y_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/blorr/inst/testfiles/blr_pairs_cpp/inputs/y");
  y_stream << y;
  std::cout << "y values: "<< y << std::endl;
  y_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    blr_pairs_cpp(x,y);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
