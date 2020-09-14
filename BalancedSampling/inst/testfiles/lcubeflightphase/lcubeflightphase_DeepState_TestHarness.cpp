#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector lcubeflightphase(NumericVector prob, NumericMatrix Xspread, NumericMatrix Xbal);

TEST(BalancedSampling_deepstate_test,lcubeflightphase_test){
  std::ofstream prob_stream;
  std::ofstream Xspread_stream;
  std::ofstream Xbal_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector prob  = RcppDeepState_NumericVector();
  prob_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/lcubeflightphase/inputs/prob");
  prob_stream << prob;
  std::cout << "prob values: "<< prob << std::endl;
  prob_stream.close();
  NumericMatrix Xspread  = RcppDeepState_NumericMatrix();
  Xspread_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/lcubeflightphase/inputs/Xspread");
  Xspread_stream << Xspread;
  std::cout << "Xspread values: "<< Xspread << std::endl;
  Xspread_stream.close();
  NumericMatrix Xbal  = RcppDeepState_NumericMatrix();
  Xbal_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/lcubeflightphase/inputs/Xbal");
  Xbal_stream << Xbal;
  std::cout << "Xbal values: "<< Xbal << std::endl;
  Xbal_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    lcubeflightphase(prob,Xspread,Xbal);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
