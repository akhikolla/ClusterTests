#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifie_comp_vcov(Rcpp::NumericMatrix parsM, Rcpp::NumericMatrix parsrepM, Rcpp::NumericMatrix Cdes, Rcpp::NumericVector rdes, Rcpp::NumericVector Ccols, Rcpp::NumericVector fayfac);

TEST(BIFIEsurvey_deepstate_test,bifie_comp_vcov_test){
  std::ofstream parsM_stream;
  std::ofstream parsrepM_stream;
  std::ofstream Cdes_stream;
  std::ofstream rdes_stream;
  std::ofstream Ccols_stream;
  std::ofstream fayfac_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix parsM  = RcppDeepState_NumericMatrix();
  parsM_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_comp_vcov/inputs/parsM");
  parsM_stream << parsM;
  std::cout << "parsM values: "<< parsM << std::endl;
  parsM_stream.close();
  NumericMatrix parsrepM  = RcppDeepState_NumericMatrix();
  parsrepM_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_comp_vcov/inputs/parsrepM");
  parsrepM_stream << parsrepM;
  std::cout << "parsrepM values: "<< parsrepM << std::endl;
  parsrepM_stream.close();
  NumericMatrix Cdes  = RcppDeepState_NumericMatrix();
  Cdes_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_comp_vcov/inputs/Cdes");
  Cdes_stream << Cdes;
  std::cout << "Cdes values: "<< Cdes << std::endl;
  Cdes_stream.close();
  NumericVector rdes  = RcppDeepState_NumericVector();
  rdes_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_comp_vcov/inputs/rdes");
  rdes_stream << rdes;
  std::cout << "rdes values: "<< rdes << std::endl;
  rdes_stream.close();
  NumericVector Ccols  = RcppDeepState_NumericVector();
  Ccols_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_comp_vcov/inputs/Ccols");
  Ccols_stream << Ccols;
  std::cout << "Ccols values: "<< Ccols << std::endl;
  Ccols_stream.close();
  NumericVector fayfac  = RcppDeepState_NumericVector();
  fayfac_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_comp_vcov/inputs/fayfac");
  fayfac_stream << fayfac;
  std::cout << "fayfac values: "<< fayfac << std::endl;
  fayfac_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifie_comp_vcov(parsM,parsrepM,Cdes,rdes,Ccols,fayfac);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
