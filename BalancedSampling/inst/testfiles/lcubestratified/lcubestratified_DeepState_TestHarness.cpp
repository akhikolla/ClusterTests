#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

IntegerVector lcubestratified(NumericVector prob, NumericMatrix Xspread, NumericMatrix Xbal, IntegerVector integerStrata);

TEST(BalancedSampling_deepstate_test,lcubestratified_test){
  std::ofstream prob_stream;
  std::ofstream Xspread_stream;
  std::ofstream Xbal_stream;
  std::ofstream integerStrata_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector prob  = RcppDeepState_NumericVector();
  prob_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/lcubestratified/inputs/prob");
  prob_stream << prob;
  std::cout << "prob values: "<< prob << std::endl;
  prob_stream.close();
  NumericMatrix Xspread  = RcppDeepState_NumericMatrix();
  Xspread_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/lcubestratified/inputs/Xspread");
  Xspread_stream << Xspread;
  std::cout << "Xspread values: "<< Xspread << std::endl;
  Xspread_stream.close();
  NumericMatrix Xbal  = RcppDeepState_NumericMatrix();
  Xbal_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/lcubestratified/inputs/Xbal");
  Xbal_stream << Xbal;
  std::cout << "Xbal values: "<< Xbal << std::endl;
  Xbal_stream.close();
  IntegerVector integerStrata  = RcppDeepState_IntegerVector();
  integerStrata_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BalancedSampling/inst/testfiles/lcubestratified/inputs/integerStrata");
  integerStrata_stream << integerStrata;
  std::cout << "integerStrata values: "<< integerStrata << std::endl;
  integerStrata_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    lcubestratified(prob,Xspread,Xbal,integerStrata);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
