#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double mise_Fh(double h, int n, NumericVector t, NumericVector w, NumericVector p, double g, int lgrid, double lim1, double lim2);

TEST(binnednp_deepstate_test,mise_Fh_test){
  std::ofstream h_stream;
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
  double h  = RcppDeepState_double();
  h_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/mise_Fh/inputs/h");
  h_stream << h;
  std::cout << "h values: "<< h << std::endl;
  h_stream.close();
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/mise_Fh/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  NumericVector t  = RcppDeepState_NumericVector();
  t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/mise_Fh/inputs/t");
  t_stream << t;
  std::cout << "t values: "<< t << std::endl;
  t_stream.close();
  NumericVector w  = RcppDeepState_NumericVector();
  w_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/mise_Fh/inputs/w");
  w_stream << w;
  std::cout << "w values: "<< w << std::endl;
  w_stream.close();
  NumericVector p  = RcppDeepState_NumericVector();
  p_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/mise_Fh/inputs/p");
  p_stream << p;
  std::cout << "p values: "<< p << std::endl;
  p_stream.close();
  double g  = RcppDeepState_double();
  g_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/mise_Fh/inputs/g");
  g_stream << g;
  std::cout << "g values: "<< g << std::endl;
  g_stream.close();
  int lgrid  = RcppDeepState_int();
  lgrid_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/mise_Fh/inputs/lgrid");
  lgrid_stream << lgrid;
  std::cout << "lgrid values: "<< lgrid << std::endl;
  lgrid_stream.close();
  double lim1  = RcppDeepState_double();
  lim1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/mise_Fh/inputs/lim1");
  lim1_stream << lim1;
  std::cout << "lim1 values: "<< lim1 << std::endl;
  lim1_stream.close();
  double lim2  = RcppDeepState_double();
  lim2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/mise_Fh/inputs/lim2");
  lim2_stream << lim2;
  std::cout << "lim2 values: "<< lim2 << std::endl;
  lim2_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    mise_Fh(h,n,t,w,p,g,lgrid,lim1,lim2);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
