#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double slope_cpp(NumericVector t, NumericVector w, double h, double lim1, double lim2, int lgrid);

TEST(binnednp_deepstate_test,slope_cpp_test){
  std::ofstream t_stream;
  std::ofstream w_stream;
  std::ofstream h_stream;
  std::ofstream lim1_stream;
  std::ofstream lim2_stream;
  std::ofstream lgrid_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector t  = RcppDeepState_NumericVector();
  t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/slope_cpp/inputs/t");
  t_stream << t;
  std::cout << "t values: "<< t << std::endl;
  t_stream.close();
  NumericVector w  = RcppDeepState_NumericVector();
  w_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/slope_cpp/inputs/w");
  w_stream << w;
  std::cout << "w values: "<< w << std::endl;
  w_stream.close();
  double h  = RcppDeepState_double();
  h_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/slope_cpp/inputs/h");
  h_stream << h;
  std::cout << "h values: "<< h << std::endl;
  h_stream.close();
  double lim1  = RcppDeepState_double();
  lim1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/slope_cpp/inputs/lim1");
  lim1_stream << lim1;
  std::cout << "lim1 values: "<< lim1 << std::endl;
  lim1_stream.close();
  double lim2  = RcppDeepState_double();
  lim2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/slope_cpp/inputs/lim2");
  lim2_stream << lim2;
  std::cout << "lim2 values: "<< lim2 << std::endl;
  lim2_stream.close();
  int lgrid  = RcppDeepState_int();
  lgrid_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/slope_cpp/inputs/lgrid");
  lgrid_stream << lgrid;
  std::cout << "lgrid values: "<< lgrid << std::endl;
  lgrid_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    slope_cpp(t,w,h,lim1,lim2,lgrid);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
