#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double Fh_combt(double x, NumericVector t, NumericVector w, double h);

TEST(binnednp_deepstate_test,Fh_combt_test){
  std::ofstream x_stream;
  std::ofstream t_stream;
  std::ofstream w_stream;
  std::ofstream h_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double x  = RcppDeepState_double();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Fh_combt/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  NumericVector t  = RcppDeepState_NumericVector();
  t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Fh_combt/inputs/t");
  t_stream << t;
  std::cout << "t values: "<< t << std::endl;
  t_stream.close();
  NumericVector w  = RcppDeepState_NumericVector();
  w_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Fh_combt/inputs/w");
  w_stream << w;
  std::cout << "w values: "<< w << std::endl;
  w_stream.close();
  double h  = RcppDeepState_double();
  h_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Fh_combt/inputs/h");
  h_stream << h;
  std::cout << "h values: "<< h << std::endl;
  h_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    Fh_combt(x,t,w,h);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
