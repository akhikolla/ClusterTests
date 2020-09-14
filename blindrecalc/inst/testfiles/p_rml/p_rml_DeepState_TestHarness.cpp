#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector p_rml(double p_c, double p_e, double r, double margin);

TEST(blindrecalc_deepstate_test,p_rml_test){
  std::ofstream p_c_stream;
  std::ofstream p_e_stream;
  std::ofstream r_stream;
  std::ofstream margin_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double p_c  = RcppDeepState_double();
  p_c_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/blindrecalc/inst/testfiles/p_rml/inputs/p_c");
  p_c_stream << p_c;
  std::cout << "p_c values: "<< p_c << std::endl;
  p_c_stream.close();
  double p_e  = RcppDeepState_double();
  p_e_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/blindrecalc/inst/testfiles/p_rml/inputs/p_e");
  p_e_stream << p_e;
  std::cout << "p_e values: "<< p_e << std::endl;
  p_e_stream.close();
  double r  = RcppDeepState_double();
  r_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/blindrecalc/inst/testfiles/p_rml/inputs/r");
  r_stream << r;
  std::cout << "r values: "<< r << std::endl;
  r_stream.close();
  double margin  = RcppDeepState_double();
  margin_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/blindrecalc/inst/testfiles/p_rml/inputs/margin");
  margin_stream << margin;
  std::cout << "margin values: "<< margin << std::endl;
  margin_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    p_rml(p_c,p_e,r,margin);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
