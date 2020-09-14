#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double pnorm_cpp(double x);

TEST(binnednp_deepstate_test,pnorm_cpp_test){
  std::ofstream x_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double x  = RcppDeepState_double();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/pnorm_cpp/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    pnorm_cpp(x);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
