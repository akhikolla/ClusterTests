#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double randIGamma(double shape, double scale);

TEST(BayesSUR_deepstate_test,randIGamma_test){
  std::ofstream shape_stream;
  std::ofstream scale_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double shape  = RcppDeepState_double();
  shape_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesSUR/inst/testfiles/randIGamma/inputs/shape");
  shape_stream << shape;
  std::cout << "shape values: "<< shape << std::endl;
  shape_stream.close();
  double scale  = RcppDeepState_double();
  scale_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesSUR/inst/testfiles/randIGamma/inputs/scale");
  scale_stream << scale;
  std::cout << "scale values: "<< scale << std::endl;
  scale_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    randIGamma(shape,scale);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
