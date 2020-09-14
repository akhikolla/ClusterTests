#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

bool myIn(int x, NumericVector y);

TEST(bigReg_deepstate_test,myIn_test){
  std::ofstream x_stream;
  std::ofstream y_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int x  = RcppDeepState_int();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bigReg/inst/testfiles/myIn/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  NumericVector y  = RcppDeepState_NumericVector();
  y_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bigReg/inst/testfiles/myIn/inputs/y");
  y_stream << y;
  std::cout << "y values: "<< y << std::endl;
  y_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    myIn(x,y);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
