#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector buildLFactors(NumericVector t, double h1, double h2);

TEST(binnednp_deepstate_test,buildLFactors_test){
  std::ofstream t_stream;
  std::ofstream h1_stream;
  std::ofstream h2_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector t  = RcppDeepState_NumericVector();
  t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/buildLFactors/inputs/t");
  t_stream << t;
  std::cout << "t values: "<< t << std::endl;
  t_stream.close();
  double h1  = RcppDeepState_double();
  h1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/buildLFactors/inputs/h1");
  h1_stream << h1;
  std::cout << "h1 values: "<< h1 << std::endl;
  h1_stream.close();
  double h2  = RcppDeepState_double();
  h2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/buildLFactors/inputs/h2");
  h2_stream << h2;
  std::cout << "h2 values: "<< h2 << std::endl;
  h2_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    buildLFactors(t,h1,h2);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
