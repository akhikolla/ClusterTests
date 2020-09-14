#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

List mfwrapper(NumericMatrix xi, IntegerVector weights, NumericMatrix qJ, IntegerVector Lv, NumericVector Eps, NumericVector priorCount);

TEST(bbl_deepstate_test,mfwrapper_test){
  std::ofstream xi_stream;
  std::ofstream weights_stream;
  std::ofstream qJ_stream;
  std::ofstream Lv_stream;
  std::ofstream Eps_stream;
  std::ofstream priorCount_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix xi  = RcppDeepState_NumericMatrix();
  xi_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bbl/inst/testfiles/mfwrapper/inputs/xi");
  xi_stream << xi;
  std::cout << "xi values: "<< xi << std::endl;
  xi_stream.close();
  IntegerVector weights  = RcppDeepState_IntegerVector();
  weights_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bbl/inst/testfiles/mfwrapper/inputs/weights");
  weights_stream << weights;
  std::cout << "weights values: "<< weights << std::endl;
  weights_stream.close();
  NumericMatrix qJ  = RcppDeepState_NumericMatrix();
  qJ_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bbl/inst/testfiles/mfwrapper/inputs/qJ");
  qJ_stream << qJ;
  std::cout << "qJ values: "<< qJ << std::endl;
  qJ_stream.close();
  IntegerVector Lv  = RcppDeepState_IntegerVector();
  Lv_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bbl/inst/testfiles/mfwrapper/inputs/Lv");
  Lv_stream << Lv;
  std::cout << "Lv values: "<< Lv << std::endl;
  Lv_stream.close();
  NumericVector Eps  = RcppDeepState_NumericVector();
  Eps_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bbl/inst/testfiles/mfwrapper/inputs/Eps");
  Eps_stream << Eps;
  std::cout << "Eps values: "<< Eps << std::endl;
  Eps_stream.close();
  NumericVector priorCount  = RcppDeepState_NumericVector();
  priorCount_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bbl/inst/testfiles/mfwrapper/inputs/priorCount");
  priorCount_stream << priorCount;
  std::cout << "priorCount values: "<< priorCount << std::endl;
  priorCount_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    mfwrapper(xi,weights,qJ,Lv,Eps,priorCount);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
