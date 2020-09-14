#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

List ARFitVARXR(NumericMatrix K21, const int k, const int p, int m, int s);

TEST(BigVAR_deepstate_test,ARFitVARXR_test){
  std::ofstream K21_stream;
  std::ofstream k_stream;
  std::ofstream p_stream;
  std::ofstream m_stream;
  std::ofstream s_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix K21  = RcppDeepState_NumericMatrix();
  K21_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BigVAR/inst/testfiles/ARFitVARXR/inputs/K21");
  K21_stream << K21;
  std::cout << "K21 values: "<< K21 << std::endl;
  K21_stream.close();
  int k  = RcppDeepState_int();
  k_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BigVAR/inst/testfiles/ARFitVARXR/inputs/k");
  k_stream << k;
  std::cout << "k values: "<< k << std::endl;
  k_stream.close();
  int p  = RcppDeepState_int();
  p_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BigVAR/inst/testfiles/ARFitVARXR/inputs/p");
  p_stream << p;
  std::cout << "p values: "<< p << std::endl;
  p_stream.close();
  int m  = RcppDeepState_int();
  m_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BigVAR/inst/testfiles/ARFitVARXR/inputs/m");
  m_stream << m;
  std::cout << "m values: "<< m << std::endl;
  m_stream.close();
  int s  = RcppDeepState_int();
  s_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BigVAR/inst/testfiles/ARFitVARXR/inputs/s");
  s_stream << s;
  std::cout << "s values: "<< s << std::endl;
  s_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    ARFitVARXR(K21,k,p,m,s);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
