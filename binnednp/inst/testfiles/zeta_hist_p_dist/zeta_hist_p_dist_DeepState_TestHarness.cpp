#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double zeta_hist_p_dist(NumericVector emp, NumericVector combt, NumericVector comby, NumericVector combw, double Af1_mixt, double l0, double l1, double h0, double h1, double lrho, double rho, int nitlambda, int nith, double lim1, double lim2);

TEST(binnednp_deepstate_test,zeta_hist_p_dist_test){
  std::ofstream emp_stream;
  std::ofstream combt_stream;
  std::ofstream comby_stream;
  std::ofstream combw_stream;
  std::ofstream Af1_mixt_stream;
  std::ofstream l0_stream;
  std::ofstream l1_stream;
  std::ofstream h0_stream;
  std::ofstream h1_stream;
  std::ofstream lrho_stream;
  std::ofstream rho_stream;
  std::ofstream nitlambda_stream;
  std::ofstream nith_stream;
  std::ofstream lim1_stream;
  std::ofstream lim2_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector emp  = RcppDeepState_NumericVector();
  emp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/emp");
  emp_stream << emp;
  std::cout << "emp values: "<< emp << std::endl;
  emp_stream.close();
  NumericVector combt  = RcppDeepState_NumericVector();
  combt_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/combt");
  combt_stream << combt;
  std::cout << "combt values: "<< combt << std::endl;
  combt_stream.close();
  NumericVector comby  = RcppDeepState_NumericVector();
  comby_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/comby");
  comby_stream << comby;
  std::cout << "comby values: "<< comby << std::endl;
  comby_stream.close();
  NumericVector combw  = RcppDeepState_NumericVector();
  combw_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/combw");
  combw_stream << combw;
  std::cout << "combw values: "<< combw << std::endl;
  combw_stream.close();
  double Af1_mixt  = RcppDeepState_double();
  Af1_mixt_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/Af1_mixt");
  Af1_mixt_stream << Af1_mixt;
  std::cout << "Af1_mixt values: "<< Af1_mixt << std::endl;
  Af1_mixt_stream.close();
  double l0  = RcppDeepState_double();
  l0_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/l0");
  l0_stream << l0;
  std::cout << "l0 values: "<< l0 << std::endl;
  l0_stream.close();
  double l1  = RcppDeepState_double();
  l1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/l1");
  l1_stream << l1;
  std::cout << "l1 values: "<< l1 << std::endl;
  l1_stream.close();
  double h0  = RcppDeepState_double();
  h0_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/h0");
  h0_stream << h0;
  std::cout << "h0 values: "<< h0 << std::endl;
  h0_stream.close();
  double h1  = RcppDeepState_double();
  h1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/h1");
  h1_stream << h1;
  std::cout << "h1 values: "<< h1 << std::endl;
  h1_stream.close();
  double lrho  = RcppDeepState_double();
  lrho_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/lrho");
  lrho_stream << lrho;
  std::cout << "lrho values: "<< lrho << std::endl;
  lrho_stream.close();
  double rho  = RcppDeepState_double();
  rho_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/rho");
  rho_stream << rho;
  std::cout << "rho values: "<< rho << std::endl;
  rho_stream.close();
  int nitlambda  = RcppDeepState_int();
  nitlambda_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/nitlambda");
  nitlambda_stream << nitlambda;
  std::cout << "nitlambda values: "<< nitlambda << std::endl;
  nitlambda_stream.close();
  int nith  = RcppDeepState_int();
  nith_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/nith");
  nith_stream << nith;
  std::cout << "nith values: "<< nith << std::endl;
  nith_stream.close();
  double lim1  = RcppDeepState_double();
  lim1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/lim1");
  lim1_stream << lim1;
  std::cout << "lim1 values: "<< lim1 << std::endl;
  lim1_stream.close();
  double lim2  = RcppDeepState_double();
  lim2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/zeta_hist_p_dist/inputs/lim2");
  lim2_stream << lim2;
  std::cout << "lim2 values: "<< lim2 << std::endl;
  lim2_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    zeta_hist_p_dist(emp,combt,comby,combw,Af1_mixt,l0,l1,h0,h1,lrho,rho,nitlambda,nith,lim1,lim2);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
