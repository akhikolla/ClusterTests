#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector psd_arma(NumericVector freq, NumericVector ar, NumericVector ma, double sigma2);

TEST(beyondWhittle_deepstate_test,psd_arma_test){
  std::ofstream freq_stream;
  std::ofstream ar_stream;
  std::ofstream ma_stream;
  std::ofstream sigma2_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector freq  = RcppDeepState_NumericVector();
  freq_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/psd_arma/inputs/freq");
  freq_stream << freq;
  std::cout << "freq values: "<< freq << std::endl;
  freq_stream.close();
  NumericVector ar  = RcppDeepState_NumericVector();
  ar_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/psd_arma/inputs/ar");
  ar_stream << ar;
  std::cout << "ar values: "<< ar << std::endl;
  ar_stream.close();
  NumericVector ma  = RcppDeepState_NumericVector();
  ma_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/psd_arma/inputs/ma");
  ma_stream << ma;
  std::cout << "ma values: "<< ma << std::endl;
  ma_stream.close();
  double sigma2  = RcppDeepState_double();
  sigma2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/psd_arma/inputs/sigma2");
  sigma2_stream << sigma2;
  std::cout << "sigma2 values: "<< sigma2 << std::endl;
  sigma2_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    psd_arma(freq,ar,ma,sigma2);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
