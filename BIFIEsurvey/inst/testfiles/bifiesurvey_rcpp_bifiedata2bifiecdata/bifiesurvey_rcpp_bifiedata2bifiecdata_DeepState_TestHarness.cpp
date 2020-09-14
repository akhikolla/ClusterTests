#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifiesurvey_rcpp_bifiedata2bifiecdata(Rcpp::NumericMatrix datalistM, int Nimp);

TEST(BIFIEsurvey_deepstate_test,bifiesurvey_rcpp_bifiedata2bifiecdata_test){
  std::ofstream datalistM_stream;
  std::ofstream Nimp_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix datalistM  = RcppDeepState_NumericMatrix();
  datalistM_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_bifiedata2bifiecdata/inputs/datalistM");
  datalistM_stream << datalistM;
  std::cout << "datalistM values: "<< datalistM << std::endl;
  datalistM_stream.close();
  int Nimp  = RcppDeepState_int();
  Nimp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_bifiedata2bifiecdata/inputs/Nimp");
  Nimp_stream << Nimp;
  std::cout << "Nimp values: "<< Nimp << std::endl;
  Nimp_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifiesurvey_rcpp_bifiedata2bifiecdata(datalistM,Nimp);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
