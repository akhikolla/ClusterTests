#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double Fg(double x, NumericVector w, NumericVector t, double g);

TEST(binnednp_deepstate_test,Fg_test){
  std::ofstream x_stream;
  std::ofstream w_stream;
  std::ofstream t_stream;
  std::ofstream g_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double x  = RcppDeepState_double();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Fg/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  NumericVector w  = RcppDeepState_NumericVector();
  w_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Fg/inputs/w");
  w_stream << w;
  std::cout << "w values: "<< w << std::endl;
  w_stream.close();
  NumericVector t  = RcppDeepState_NumericVector();
  t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Fg/inputs/t");
  t_stream << t;
  std::cout << "t values: "<< t << std::endl;
  t_stream.close();
  double g  = RcppDeepState_double();
  g_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Fg/inputs/g");
  g_stream << g;
  std::cout << "g values: "<< g << std::endl;
  g_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    Fg(x,w,t,g);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
