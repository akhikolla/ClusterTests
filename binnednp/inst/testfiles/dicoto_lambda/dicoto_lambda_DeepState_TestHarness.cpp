#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector dicoto_lambda(double lambda, int nith, double h0, double h1, double rho, NumericVector hist, NumericVector combt, NumericVector combw, double lim1, double lim2);

TEST(binnednp_deepstate_test,dicoto_lambda_test){
  std::ofstream lambda_stream;
  std::ofstream nith_stream;
  std::ofstream h0_stream;
  std::ofstream h1_stream;
  std::ofstream rho_stream;
  std::ofstream hist_stream;
  std::ofstream combt_stream;
  std::ofstream combw_stream;
  std::ofstream lim1_stream;
  std::ofstream lim2_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double lambda  = RcppDeepState_double();
  lambda_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/dicoto_lambda/inputs/lambda");
  lambda_stream << lambda;
  std::cout << "lambda values: "<< lambda << std::endl;
  lambda_stream.close();
  int nith  = RcppDeepState_int();
  nith_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/dicoto_lambda/inputs/nith");
  nith_stream << nith;
  std::cout << "nith values: "<< nith << std::endl;
  nith_stream.close();
  double h0  = RcppDeepState_double();
  h0_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/dicoto_lambda/inputs/h0");
  h0_stream << h0;
  std::cout << "h0 values: "<< h0 << std::endl;
  h0_stream.close();
  double h1  = RcppDeepState_double();
  h1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/dicoto_lambda/inputs/h1");
  h1_stream << h1;
  std::cout << "h1 values: "<< h1 << std::endl;
  h1_stream.close();
  double rho  = RcppDeepState_double();
  rho_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/dicoto_lambda/inputs/rho");
  rho_stream << rho;
  std::cout << "rho values: "<< rho << std::endl;
  rho_stream.close();
  NumericVector hist  = RcppDeepState_NumericVector();
  hist_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/dicoto_lambda/inputs/hist");
  hist_stream << hist;
  std::cout << "hist values: "<< hist << std::endl;
  hist_stream.close();
  NumericVector combt  = RcppDeepState_NumericVector();
  combt_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/dicoto_lambda/inputs/combt");
  combt_stream << combt;
  std::cout << "combt values: "<< combt << std::endl;
  combt_stream.close();
  NumericVector combw  = RcppDeepState_NumericVector();
  combw_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/dicoto_lambda/inputs/combw");
  combw_stream << combw;
  std::cout << "combw values: "<< combw << std::endl;
  combw_stream.close();
  double lim1  = RcppDeepState_double();
  lim1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/dicoto_lambda/inputs/lim1");
  lim1_stream << lim1;
  std::cout << "lim1 values: "<< lim1 << std::endl;
  lim1_stream.close();
  double lim2  = RcppDeepState_double();
  lim2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/dicoto_lambda/inputs/lim2");
  lim2_stream << lim2;
  std::cout << "lim2 values: "<< lim2 << std::endl;
  lim2_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    dicoto_lambda(lambda,nith,h0,h1,rho,hist,combt,combw,lim1,lim2);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
