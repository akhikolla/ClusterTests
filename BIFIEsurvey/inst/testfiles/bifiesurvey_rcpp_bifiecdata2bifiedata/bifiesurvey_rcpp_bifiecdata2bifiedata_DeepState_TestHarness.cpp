#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifiesurvey_rcpp_bifiecdata2bifiedata(Rcpp::NumericMatrix datalistM_ind, Rcpp::NumericMatrix datalistM_imputed, int Nimp, Rcpp::NumericMatrix dat1, Rcpp::NumericMatrix datalistM_impindex);

TEST(BIFIEsurvey_deepstate_test,bifiesurvey_rcpp_bifiecdata2bifiedata_test){
  std::ofstream datalistM_ind_stream;
  std::ofstream datalistM_imputed_stream;
  std::ofstream Nimp_stream;
  std::ofstream dat1_stream;
  std::ofstream datalistM_impindex_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix datalistM_ind  = RcppDeepState_NumericMatrix();
  datalistM_ind_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_bifiecdata2bifiedata/inputs/datalistM_ind");
  datalistM_ind_stream << datalistM_ind;
  std::cout << "datalistM_ind values: "<< datalistM_ind << std::endl;
  datalistM_ind_stream.close();
  NumericMatrix datalistM_imputed  = RcppDeepState_NumericMatrix();
  datalistM_imputed_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_bifiecdata2bifiedata/inputs/datalistM_imputed");
  datalistM_imputed_stream << datalistM_imputed;
  std::cout << "datalistM_imputed values: "<< datalistM_imputed << std::endl;
  datalistM_imputed_stream.close();
  int Nimp  = RcppDeepState_int();
  Nimp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_bifiecdata2bifiedata/inputs/Nimp");
  Nimp_stream << Nimp;
  std::cout << "Nimp values: "<< Nimp << std::endl;
  Nimp_stream.close();
  NumericMatrix dat1  = RcppDeepState_NumericMatrix();
  dat1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_bifiecdata2bifiedata/inputs/dat1");
  dat1_stream << dat1;
  std::cout << "dat1 values: "<< dat1 << std::endl;
  dat1_stream.close();
  NumericMatrix datalistM_impindex  = RcppDeepState_NumericMatrix();
  datalistM_impindex_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_bifiecdata2bifiedata/inputs/datalistM_impindex");
  datalistM_impindex_stream << datalistM_impindex;
  std::cout << "datalistM_impindex values: "<< datalistM_impindex << std::endl;
  datalistM_impindex_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifiesurvey_rcpp_bifiecdata2bifiedata(datalistM_ind,datalistM_imputed,Nimp,dat1,datalistM_impindex);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
