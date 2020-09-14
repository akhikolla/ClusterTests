#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double PoissonClosed_(double alpha_1, double beta_1, double alpha_2, double beta_2);

TEST(bayesAB_deepstate_test,PoissonClosed__test){
  std::ofstream alpha_1_stream;
  std::ofstream beta_1_stream;
  std::ofstream alpha_2_stream;
  std::ofstream beta_2_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double alpha_1  = RcppDeepState_double();
  alpha_1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bayesAB/inst/testfiles/PoissonClosed_/inputs/alpha_1");
  alpha_1_stream << alpha_1;
  std::cout << "alpha_1 values: "<< alpha_1 << std::endl;
  alpha_1_stream.close();
  double beta_1  = RcppDeepState_double();
  beta_1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bayesAB/inst/testfiles/PoissonClosed_/inputs/beta_1");
  beta_1_stream << beta_1;
  std::cout << "beta_1 values: "<< beta_1 << std::endl;
  beta_1_stream.close();
  double alpha_2  = RcppDeepState_double();
  alpha_2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bayesAB/inst/testfiles/PoissonClosed_/inputs/alpha_2");
  alpha_2_stream << alpha_2;
  std::cout << "alpha_2 values: "<< alpha_2 << std::endl;
  alpha_2_stream.close();
  double beta_2  = RcppDeepState_double();
  beta_2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bayesAB/inst/testfiles/PoissonClosed_/inputs/beta_2");
  beta_2_stream << beta_2;
  std::cout << "beta_2 values: "<< beta_2 << std::endl;
  beta_2_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    PoissonClosed_(alpha_1,beta_1,alpha_2,beta_2);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
