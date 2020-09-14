#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

IntegerVector lcubelandingphase(NumericVector prob, NumericVector probflight, NumericMatrix Xspread, NumericMatrix Xbal);

TEST(BalancedSampling_deepstate_test,lcubelandingphase_test){
  std::ofstream prob_stream;
  std::ofstream probflight_stream;
  std::ofstream Xspread_stream;
  std::ofstream Xbal_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector prob  = RcppDeepState_NumericVector();
  prob_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/lcubelandingphase/inputs/prob");
  prob_stream << prob;
  std::cout << "prob values: "<< prob << std::endl;
  prob_stream.close();
  NumericVector probflight  = RcppDeepState_NumericVector();
  probflight_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/lcubelandingphase/inputs/probflight");
  probflight_stream << probflight;
  std::cout << "probflight values: "<< probflight << std::endl;
  probflight_stream.close();
  NumericMatrix Xspread  = RcppDeepState_NumericMatrix();
  Xspread_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/lcubelandingphase/inputs/Xspread");
  Xspread_stream << Xspread;
  std::cout << "Xspread values: "<< Xspread << std::endl;
  Xspread_stream.close();
  NumericMatrix Xbal  = RcppDeepState_NumericMatrix();
  Xbal_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/lcubelandingphase/inputs/Xbal");
  Xbal_stream << Xbal;
  std::cout << "Xbal values: "<< Xbal << std::endl;
  Xbal_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    lcubelandingphase(prob,probflight,Xspread,Xbal);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
