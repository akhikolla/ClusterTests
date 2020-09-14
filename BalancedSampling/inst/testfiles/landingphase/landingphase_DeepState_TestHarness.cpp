#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

IntegerVector landingphase(NumericVector prob, NumericVector probflight, NumericMatrix Xbal);

TEST(BalancedSampling_deepstate_test,landingphase_test){
  std::ofstream prob_stream;
  std::ofstream probflight_stream;
  std::ofstream Xbal_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector prob  = RcppDeepState_NumericVector();
  prob_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/landingphase/inputs/prob");
  prob_stream << prob;
  std::cout << "prob values: "<< prob << std::endl;
  prob_stream.close();
  NumericVector probflight  = RcppDeepState_NumericVector();
  probflight_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/landingphase/inputs/probflight");
  probflight_stream << probflight;
  std::cout << "probflight values: "<< probflight << std::endl;
  probflight_stream.close();
  NumericMatrix Xbal  = RcppDeepState_NumericMatrix();
  Xbal_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/landingphase/inputs/Xbal");
  Xbal_stream << Xbal;
  std::cout << "Xbal values: "<< Xbal << std::endl;
  Xbal_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    landingphase(prob,probflight,Xbal);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
