#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

IntegerVector myc(IntegerVector x, IntegerVector y);

TEST(BClustLonG_deepstate_test,myc_test){
  std::ofstream x_stream;
  std::ofstream y_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  IntegerVector x  = RcppDeepState_IntegerVector();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BClustLonG/inst/testfiles/myc/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  IntegerVector y  = RcppDeepState_IntegerVector();
  y_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BClustLonG/inst/testfiles/myc/inputs/y");
  y_stream << y;
  std::cout << "y values: "<< y << std::endl;
  y_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    myc(x,y);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
