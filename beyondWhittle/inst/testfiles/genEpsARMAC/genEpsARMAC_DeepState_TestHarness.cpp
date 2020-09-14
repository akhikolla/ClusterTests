#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector genEpsARMAC(NumericVector zt, NumericVector ar, NumericVector ma);

TEST(beyondWhittle_deepstate_test,genEpsARMAC_test){
  std::ofstream zt_stream;
  std::ofstream ar_stream;
  std::ofstream ma_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector zt  = RcppDeepState_NumericVector();
  zt_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/genEpsARMAC/inputs/zt");
  zt_stream << zt;
  std::cout << "zt values: "<< zt << std::endl;
  zt_stream.close();
  NumericVector ar  = RcppDeepState_NumericVector();
  ar_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/genEpsARMAC/inputs/ar");
  ar_stream << ar;
  std::cout << "ar values: "<< ar << std::endl;
  ar_stream.close();
  NumericVector ma  = RcppDeepState_NumericVector();
  ma_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/genEpsARMAC/inputs/ma");
  ma_stream << ma;
  std::cout << "ma values: "<< ma << std::endl;
  ma_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    genEpsARMAC(zt,ar,ma);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
