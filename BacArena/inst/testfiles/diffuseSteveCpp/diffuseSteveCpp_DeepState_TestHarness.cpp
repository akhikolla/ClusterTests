#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

void diffuseSteveCpp(Rcpp::NumericMatrix y, double D, double h, double tstep);

TEST(BacArena_deepstate_test,diffuseSteveCpp_test){
  std::ofstream y_stream;
  std::ofstream D_stream;
  std::ofstream h_stream;
  std::ofstream tstep_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix y  = RcppDeepState_NumericMatrix();
  y_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BacArena/inst/testfiles/diffuseSteveCpp/inputs/y");
  y_stream << y;
  std::cout << "y values: "<< y << std::endl;
  y_stream.close();
  double D  = RcppDeepState_double();
  D_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BacArena/inst/testfiles/diffuseSteveCpp/inputs/D");
  D_stream << D;
  std::cout << "D values: "<< D << std::endl;
  D_stream.close();
  double h  = RcppDeepState_double();
  h_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BacArena/inst/testfiles/diffuseSteveCpp/inputs/h");
  h_stream << h;
  std::cout << "h values: "<< h << std::endl;
  h_stream.close();
  double tstep  = RcppDeepState_double();
  tstep_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BacArena/inst/testfiles/diffuseSteveCpp/inputs/tstep");
  tstep_stream << tstep;
  std::cout << "tstep values: "<< tstep << std::endl;
  tstep_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    diffuseSteveCpp(y,D,h,tstep);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
