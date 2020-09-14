#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double randT(const double nu);

TEST(BayesSUR_deepstate_test,randT_test){
  std::ofstream nu_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double nu  = RcppDeepState_double();
  nu_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesSUR/inst/testfiles/randT/inputs/nu");
  nu_stream << nu;
  std::cout << "nu values: "<< nu << std::endl;
  nu_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    randT(nu);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
