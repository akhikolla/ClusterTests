#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifie_comp_vcov_within(Rcpp::NumericMatrix parsM, Rcpp::NumericMatrix parsrepM, Rcpp::NumericVector fayfac, int RR, int Nimp);

TEST(BIFIEsurvey_deepstate_test,bifie_comp_vcov_within_test){
  std::ofstream parsM_stream;
  std::ofstream parsrepM_stream;
  std::ofstream fayfac_stream;
  std::ofstream RR_stream;
  std::ofstream Nimp_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix parsM  = RcppDeepState_NumericMatrix();
  parsM_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_comp_vcov_within/inputs/parsM");
  parsM_stream << parsM;
  std::cout << "parsM values: "<< parsM << std::endl;
  parsM_stream.close();
  NumericMatrix parsrepM  = RcppDeepState_NumericMatrix();
  parsrepM_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_comp_vcov_within/inputs/parsrepM");
  parsrepM_stream << parsrepM;
  std::cout << "parsrepM values: "<< parsrepM << std::endl;
  parsrepM_stream.close();
  NumericVector fayfac  = RcppDeepState_NumericVector();
  fayfac_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_comp_vcov_within/inputs/fayfac");
  fayfac_stream << fayfac;
  std::cout << "fayfac values: "<< fayfac << std::endl;
  fayfac_stream.close();
  int RR  = RcppDeepState_int();
  RR_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_comp_vcov_within/inputs/RR");
  RR_stream << RR;
  std::cout << "RR values: "<< RR << std::endl;
  RR_stream.close();
  int Nimp  = RcppDeepState_int();
  Nimp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_comp_vcov_within/inputs/Nimp");
  Nimp_stream << Nimp;
  std::cout << "Nimp values: "<< Nimp << std::endl;
  Nimp_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifie_comp_vcov_within(parsM,parsrepM,fayfac,RR,Nimp);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
