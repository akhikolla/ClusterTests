#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector vFromP(NumericVector p, const double eps);

TEST(beyondWhittle_deepstate_test,vFromP_test){
  std::ofstream p_stream;
  std::ofstream eps_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector p  = RcppDeepState_NumericVector();
  p_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/vFromP/inputs/p");
  p_stream << p;
  std::cout << "p values: "<< p << std::endl;
  p_stream.close();
  double eps  = RcppDeepState_double();
  eps_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/vFromP/inputs/eps");
  eps_stream << eps;
  std::cout << "eps values: "<< eps << std::endl;
  eps_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    vFromP(p,eps);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
