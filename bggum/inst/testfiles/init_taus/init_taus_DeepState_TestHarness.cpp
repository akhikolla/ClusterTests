#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

List init_taus(int m, double shape1, double shape2, double a, double b, IntegerVector K);

TEST(bggum_deepstate_test,init_taus_test){
  std::ofstream m_stream;
  std::ofstream shape1_stream;
  std::ofstream shape2_stream;
  std::ofstream a_stream;
  std::ofstream b_stream;
  std::ofstream K_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int m  = RcppDeepState_int();
  m_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bggum/inst/testfiles/init_taus/inputs/m");
  m_stream << m;
  std::cout << "m values: "<< m << std::endl;
  m_stream.close();
  double shape1  = RcppDeepState_double();
  shape1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bggum/inst/testfiles/init_taus/inputs/shape1");
  shape1_stream << shape1;
  std::cout << "shape1 values: "<< shape1 << std::endl;
  shape1_stream.close();
  double shape2  = RcppDeepState_double();
  shape2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bggum/inst/testfiles/init_taus/inputs/shape2");
  shape2_stream << shape2;
  std::cout << "shape2 values: "<< shape2 << std::endl;
  shape2_stream.close();
  double a  = RcppDeepState_double();
  a_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bggum/inst/testfiles/init_taus/inputs/a");
  a_stream << a;
  std::cout << "a values: "<< a << std::endl;
  a_stream.close();
  double b  = RcppDeepState_double();
  b_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bggum/inst/testfiles/init_taus/inputs/b");
  b_stream << b;
  std::cout << "b values: "<< b << std::endl;
  b_stream.close();
  IntegerVector K  = RcppDeepState_IntegerVector();
  K_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bggum/inst/testfiles/init_taus/inputs/K");
  K_stream << K;
  std::cout << "K values: "<< K << std::endl;
  K_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    init_taus(m,shape1,shape2,a,b,K);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
