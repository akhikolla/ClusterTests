#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifiesurvey_rcpp_bootstrap(Rcpp::NumericVector cumwgt, Rcpp::NumericMatrix rand_wgt);

TEST(BIFIEsurvey_deepstate_test,bifiesurvey_rcpp_bootstrap_test){
  std::ofstream cumwgt_stream;
  std::ofstream rand_wgt_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector cumwgt  = RcppDeepState_NumericVector();
  cumwgt_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_bootstrap/inputs/cumwgt");
  cumwgt_stream << cumwgt;
  std::cout << "cumwgt values: "<< cumwgt << std::endl;
  cumwgt_stream.close();
  NumericMatrix rand_wgt  = RcppDeepState_NumericMatrix();
  rand_wgt_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_bootstrap/inputs/rand_wgt");
  rand_wgt_stream << rand_wgt;
  std::cout << "rand_wgt values: "<< rand_wgt << std::endl;
  rand_wgt_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifiesurvey_rcpp_bootstrap(cumwgt,rand_wgt);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
