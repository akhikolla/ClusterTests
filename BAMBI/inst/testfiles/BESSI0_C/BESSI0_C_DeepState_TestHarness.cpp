#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double BESSI0_C(double x);

TEST(BAMBI_deepstate_test,BESSI0_C_test){
  std::ofstream x_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double x  = RcppDeepState_double();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/BESSI0_C/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    BESSI0_C(x);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
