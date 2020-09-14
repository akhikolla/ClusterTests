#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double gaussian_mise_loop(int hn, NumericVector seq, double rho);

TEST(binnednp_deepstate_test,gaussian_mise_loop_test){
  std::ofstream hn_stream;
  std::ofstream seq_stream;
  std::ofstream rho_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int hn  = RcppDeepState_int();
  hn_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/gaussian_mise_loop/inputs/hn");
  hn_stream << hn;
  std::cout << "hn values: "<< hn << std::endl;
  hn_stream.close();
  NumericVector seq  = RcppDeepState_NumericVector();
  seq_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/gaussian_mise_loop/inputs/seq");
  seq_stream << seq;
  std::cout << "seq values: "<< seq << std::endl;
  seq_stream.close();
  double rho  = RcppDeepState_double();
  rho_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/gaussian_mise_loop/inputs/rho");
  rho_stream << rho;
  std::cout << "rho values: "<< rho << std::endl;
  rho_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    gaussian_mise_loop(hn,seq,rho);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
