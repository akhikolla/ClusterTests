#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifiesurvey_rcpp_bifiedata_stepwise(Rcpp::NumericMatrix dat1, Rcpp::NumericMatrix dat_ind, int Nmiss);

TEST(BIFIEsurvey_deepstate_test,bifiesurvey_rcpp_bifiedata_stepwise_test){
  std::ofstream dat1_stream;
  std::ofstream dat_ind_stream;
  std::ofstream Nmiss_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix dat1  = RcppDeepState_NumericMatrix();
  dat1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_bifiedata_stepwise/inputs/dat1");
  dat1_stream << dat1;
  std::cout << "dat1 values: "<< dat1 << std::endl;
  dat1_stream.close();
  NumericMatrix dat_ind  = RcppDeepState_NumericMatrix();
  dat_ind_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_bifiedata_stepwise/inputs/dat_ind");
  dat_ind_stream << dat_ind;
  std::cout << "dat_ind values: "<< dat_ind << std::endl;
  dat_ind_stream.close();
  int Nmiss  = RcppDeepState_int();
  Nmiss_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_bifiedata_stepwise/inputs/Nmiss");
  Nmiss_stream << Nmiss;
  std::cout << "Nmiss values: "<< Nmiss << std::endl;
  Nmiss_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifiesurvey_rcpp_bifiedata_stepwise(dat1,dat_ind,Nmiss);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
