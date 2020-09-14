#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifiesurvey_rcpp_rubin_rules(Rcpp::NumericMatrix estimates, Rcpp::NumericMatrix variances);

TEST(BIFIEsurvey_deepstate_test,bifiesurvey_rcpp_rubin_rules_test){
  std::ofstream estimates_stream;
  std::ofstream variances_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix estimates  = RcppDeepState_NumericMatrix();
  estimates_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_rubin_rules/inputs/estimates");
  estimates_stream << estimates;
  std::cout << "estimates values: "<< estimates << std::endl;
  estimates_stream.close();
  NumericMatrix variances  = RcppDeepState_NumericMatrix();
  variances_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_rubin_rules/inputs/variances");
  variances_stream << variances;
  std::cout << "variances values: "<< variances << std::endl;
  variances_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifiesurvey_rcpp_rubin_rules(estimates,variances);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
