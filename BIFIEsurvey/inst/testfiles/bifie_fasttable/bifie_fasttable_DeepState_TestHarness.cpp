#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifie_fasttable(Rcpp::NumericMatrix datavec);

TEST(BIFIEsurvey_deepstate_test,bifie_fasttable_test){
  std::ofstream datavec_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix datavec  = RcppDeepState_NumericMatrix();
  datavec_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_fasttable/inputs/datavec");
  datavec_stream << datavec;
  std::cout << "datavec values: "<< datavec << std::endl;
  datavec_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifie_fasttable(datavec);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
