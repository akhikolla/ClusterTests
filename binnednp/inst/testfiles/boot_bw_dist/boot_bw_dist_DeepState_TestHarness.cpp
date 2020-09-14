#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double boot_bw_dist(int nit, double h0, double h1, double rho, int n, NumericVector t, NumericVector w, NumericVector p, double g, int lgrid, double lim1, double lim2);

TEST(binnednp_deepstate_test,boot_bw_dist_test){
  std::ofstream nit_stream;
  std::ofstream h0_stream;
  std::ofstream h1_stream;
  std::ofstream rho_stream;
  std::ofstream n_stream;
  std::ofstream t_stream;
  std::ofstream w_stream;
  std::ofstream p_stream;
  std::ofstream g_stream;
  std::ofstream lgrid_stream;
  std::ofstream lim1_stream;
  std::ofstream lim2_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int nit  = RcppDeepState_int();
  nit_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/boot_bw_dist/inputs/nit");
  nit_stream << nit;
  std::cout << "nit values: "<< nit << std::endl;
  nit_stream.close();
  double h0  = RcppDeepState_double();
  h0_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/boot_bw_dist/inputs/h0");
  h0_stream << h0;
  std::cout << "h0 values: "<< h0 << std::endl;
  h0_stream.close();
  double h1  = RcppDeepState_double();
  h1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/boot_bw_dist/inputs/h1");
  h1_stream << h1;
  std::cout << "h1 values: "<< h1 << std::endl;
  h1_stream.close();
  double rho  = RcppDeepState_double();
  rho_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/boot_bw_dist/inputs/rho");
  rho_stream << rho;
  std::cout << "rho values: "<< rho << std::endl;
  rho_stream.close();
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/boot_bw_dist/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  NumericVector t  = RcppDeepState_NumericVector();
  t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/boot_bw_dist/inputs/t");
  t_stream << t;
  std::cout << "t values: "<< t << std::endl;
  t_stream.close();
  NumericVector w  = RcppDeepState_NumericVector();
  w_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/boot_bw_dist/inputs/w");
  w_stream << w;
  std::cout << "w values: "<< w << std::endl;
  w_stream.close();
  NumericVector p  = RcppDeepState_NumericVector();
  p_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/boot_bw_dist/inputs/p");
  p_stream << p;
  std::cout << "p values: "<< p << std::endl;
  p_stream.close();
  double g  = RcppDeepState_double();
  g_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/boot_bw_dist/inputs/g");
  g_stream << g;
  std::cout << "g values: "<< g << std::endl;
  g_stream.close();
  int lgrid  = RcppDeepState_int();
  lgrid_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/boot_bw_dist/inputs/lgrid");
  lgrid_stream << lgrid;
  std::cout << "lgrid values: "<< lgrid << std::endl;
  lgrid_stream.close();
  double lim1  = RcppDeepState_double();
  lim1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/boot_bw_dist/inputs/lim1");
  lim1_stream << lim1;
  std::cout << "lim1 values: "<< lim1 << std::endl;
  lim1_stream.close();
  double lim2  = RcppDeepState_double();
  lim2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/boot_bw_dist/inputs/lim2");
  lim2_stream << lim2;
  std::cout << "lim2 values: "<< lim2 << std::endl;
  lim2_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    boot_bw_dist(nit,h0,h1,rho,n,t,w,p,g,lgrid,lim1,lim2);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
