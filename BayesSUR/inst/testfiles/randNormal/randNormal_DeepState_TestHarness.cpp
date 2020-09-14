#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double randNormal(const double m, const double sigmaSquare);

TEST(BayesSUR_deepstate_test,randNormal_test){
  std::ofstream m_stream;
  std::ofstream sigmaSquare_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double m  = RcppDeepState_double();
  m_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesSUR/inst/testfiles/randNormal/inputs/m");
  m_stream << m;
  std::cout << "m values: "<< m << std::endl;
  m_stream.close();
  double sigmaSquare  = RcppDeepState_double();
  sigmaSquare_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesSUR/inst/testfiles/randNormal/inputs/sigmaSquare");
  sigmaSquare_stream << sigmaSquare;
  std::cout << "sigmaSquare values: "<< sigmaSquare << std::endl;
  sigmaSquare_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    randNormal(m,sigmaSquare);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
