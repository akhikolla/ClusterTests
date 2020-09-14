#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

SEXP r_bind(arma::mat x, arma::mat y);

TEST(bigReg_deepstate_test,r_bind_test){
  std::ofstream x_stream;
  std::ofstream y_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat x  = RcppDeepState_mat();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bigReg/inst/testfiles/r_bind/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  arma::mat y  = RcppDeepState_mat();
  y_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bigReg/inst/testfiles/r_bind/inputs/y");
  y_stream << y;
  std::cout << "y values: "<< y << std::endl;
  y_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    r_bind(x,y);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
