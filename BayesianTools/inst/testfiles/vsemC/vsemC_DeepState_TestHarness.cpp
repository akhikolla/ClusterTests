#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix vsemC(NumericVector par, NumericVector PAR);

TEST(BayesianTools_deepstate_test,vsemC_test){
  std::ofstream par_stream;
  std::ofstream PAR_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector par  = RcppDeepState_NumericVector();
  par_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesianTools/inst/testfiles/vsemC/inputs/par");
  par_stream << par;
  std::cout << "par values: "<< par << std::endl;
  par_stream.close();
  NumericVector PAR  = RcppDeepState_NumericVector();
  PAR_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BayesianTools/inst/testfiles/vsemC/inputs/PAR");
  PAR_stream << PAR;
  std::cout << "PAR values: "<< PAR << std::endl;
  PAR_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    vsemC(par,PAR);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
