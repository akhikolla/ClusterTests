#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

void gaussian_mise_initialize(int n, int k, NumericVector wg, NumericVector w, NumericVector t, double gboot, double AK);

TEST(binnednp_deepstate_test,gaussian_mise_initialize_test){
  std::ofstream n_stream;
  std::ofstream k_stream;
  std::ofstream wg_stream;
  std::ofstream w_stream;
  std::ofstream t_stream;
  std::ofstream gboot_stream;
  std::ofstream AK_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/gaussian_mise_initialize/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  int k  = RcppDeepState_int();
  k_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/gaussian_mise_initialize/inputs/k");
  k_stream << k;
  std::cout << "k values: "<< k << std::endl;
  k_stream.close();
  NumericVector wg  = RcppDeepState_NumericVector();
  wg_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/gaussian_mise_initialize/inputs/wg");
  wg_stream << wg;
  std::cout << "wg values: "<< wg << std::endl;
  wg_stream.close();
  NumericVector w  = RcppDeepState_NumericVector();
  w_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/gaussian_mise_initialize/inputs/w");
  w_stream << w;
  std::cout << "w values: "<< w << std::endl;
  w_stream.close();
  NumericVector t  = RcppDeepState_NumericVector();
  t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/gaussian_mise_initialize/inputs/t");
  t_stream << t;
  std::cout << "t values: "<< t << std::endl;
  t_stream.close();
  double gboot  = RcppDeepState_double();
  gboot_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/gaussian_mise_initialize/inputs/gboot");
  gboot_stream << gboot;
  std::cout << "gboot values: "<< gboot << std::endl;
  gboot_stream.close();
  double AK  = RcppDeepState_double();
  AK_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/gaussian_mise_initialize/inputs/AK");
  AK_stream << AK;
  std::cout << "AK values: "<< AK << std::endl;
  AK_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    gaussian_mise_initialize(n,k,wg,w,t,gboot,AK);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
