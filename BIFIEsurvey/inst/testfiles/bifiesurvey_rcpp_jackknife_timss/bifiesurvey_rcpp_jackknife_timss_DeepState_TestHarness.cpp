#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::NumericMatrix bifiesurvey_rcpp_jackknife_timss(Rcpp::NumericVector wgt, Rcpp::NumericVector jkzone, Rcpp::NumericVector jkrep, int RR, double jkfac, Rcpp::NumericVector prbar);

TEST(BIFIEsurvey_deepstate_test,bifiesurvey_rcpp_jackknife_timss_test){
  std::ofstream wgt_stream;
  std::ofstream jkzone_stream;
  std::ofstream jkrep_stream;
  std::ofstream RR_stream;
  std::ofstream jkfac_stream;
  std::ofstream prbar_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector wgt  = RcppDeepState_NumericVector();
  wgt_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_jackknife_timss/inputs/wgt");
  wgt_stream << wgt;
  std::cout << "wgt values: "<< wgt << std::endl;
  wgt_stream.close();
  NumericVector jkzone  = RcppDeepState_NumericVector();
  jkzone_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_jackknife_timss/inputs/jkzone");
  jkzone_stream << jkzone;
  std::cout << "jkzone values: "<< jkzone << std::endl;
  jkzone_stream.close();
  NumericVector jkrep  = RcppDeepState_NumericVector();
  jkrep_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_jackknife_timss/inputs/jkrep");
  jkrep_stream << jkrep;
  std::cout << "jkrep values: "<< jkrep << std::endl;
  jkrep_stream.close();
  int RR  = RcppDeepState_int();
  RR_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_jackknife_timss/inputs/RR");
  RR_stream << RR;
  std::cout << "RR values: "<< RR << std::endl;
  RR_stream.close();
  double jkfac  = RcppDeepState_double();
  jkfac_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_jackknife_timss/inputs/jkfac");
  jkfac_stream << jkfac;
  std::cout << "jkfac values: "<< jkfac << std::endl;
  jkfac_stream.close();
  NumericVector prbar  = RcppDeepState_NumericVector();
  prbar_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_jackknife_timss/inputs/prbar");
  prbar_stream << prbar;
  std::cout << "prbar values: "<< prbar << std::endl;
  prbar_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifiesurvey_rcpp_jackknife_timss(wgt,jkzone,jkrep,RR,jkfac,prbar);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
