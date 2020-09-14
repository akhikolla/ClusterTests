#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::mat rsin_unimodal(int n, double k1, double k2, double k3, double mu1, double mu2, double kappa_opt, double log_I0_kappa_opt, double logK, double log_const_vmsin);

TEST(BAMBI_deepstate_test,rsin_unimodal_test){
  std::ofstream n_stream;
  std::ofstream k1_stream;
  std::ofstream k2_stream;
  std::ofstream k3_stream;
  std::ofstream mu1_stream;
  std::ofstream mu2_stream;
  std::ofstream kappa_opt_stream;
  std::ofstream log_I0_kappa_opt_stream;
  std::ofstream logK_stream;
  std::ofstream log_const_vmsin_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/rsin_unimodal/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  double k1  = RcppDeepState_double();
  k1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/rsin_unimodal/inputs/k1");
  k1_stream << k1;
  std::cout << "k1 values: "<< k1 << std::endl;
  k1_stream.close();
  double k2  = RcppDeepState_double();
  k2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/rsin_unimodal/inputs/k2");
  k2_stream << k2;
  std::cout << "k2 values: "<< k2 << std::endl;
  k2_stream.close();
  double k3  = RcppDeepState_double();
  k3_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/rsin_unimodal/inputs/k3");
  k3_stream << k3;
  std::cout << "k3 values: "<< k3 << std::endl;
  k3_stream.close();
  double mu1  = RcppDeepState_double();
  mu1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/rsin_unimodal/inputs/mu1");
  mu1_stream << mu1;
  std::cout << "mu1 values: "<< mu1 << std::endl;
  mu1_stream.close();
  double mu2  = RcppDeepState_double();
  mu2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/rsin_unimodal/inputs/mu2");
  mu2_stream << mu2;
  std::cout << "mu2 values: "<< mu2 << std::endl;
  mu2_stream.close();
  double kappa_opt  = RcppDeepState_double();
  kappa_opt_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/rsin_unimodal/inputs/kappa_opt");
  kappa_opt_stream << kappa_opt;
  std::cout << "kappa_opt values: "<< kappa_opt << std::endl;
  kappa_opt_stream.close();
  double log_I0_kappa_opt  = RcppDeepState_double();
  log_I0_kappa_opt_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/rsin_unimodal/inputs/log_I0_kappa_opt");
  log_I0_kappa_opt_stream << log_I0_kappa_opt;
  std::cout << "log_I0_kappa_opt values: "<< log_I0_kappa_opt << std::endl;
  log_I0_kappa_opt_stream.close();
  double logK  = RcppDeepState_double();
  logK_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/rsin_unimodal/inputs/logK");
  logK_stream << logK;
  std::cout << "logK values: "<< logK << std::endl;
  logK_stream.close();
  double log__vmsin  = RcppDeepState_double();
  log_const_vmsin_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAMBI/inst/testfiles/rsin_unimodal/inputs/log_const_vmsin");
  log_const_vmsin_stream << log_const_vmsin;
  std::cout << "log_const_vmsin values: "<< log_const_vmsin << std::endl;
  log_const_vmsin_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    rsin_unimodal(n,k1,k2,k3,mu1,mu2,kappa_opt,log_I0_kappa_opt,logK,log_const_vmsin);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
