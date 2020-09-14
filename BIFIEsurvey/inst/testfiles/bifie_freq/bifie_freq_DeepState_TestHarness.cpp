#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifie_freq(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values, Rcpp::NumericMatrix vars_values, Rcpp::NumericVector vars_values_numb);

TEST(BIFIEsurvey_deepstate_test,bifie_freq_test){
  std::ofstream datalist_stream;
  std::ofstream wgt1_stream;
  std::ofstream wgtrep_stream;
  std::ofstream vars_index_stream;
  std::ofstream fayfac_stream;
  std::ofstream NI_stream;
  std::ofstream group_index1_stream;
  std::ofstream group_values_stream;
  std::ofstream vars_values_stream;
  std::ofstream vars_values_numb_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix datalist  = RcppDeepState_NumericMatrix();
  datalist_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_freq/inputs/datalist");
  datalist_stream << datalist;
  std::cout << "datalist values: "<< datalist << std::endl;
  datalist_stream.close();
  NumericMatrix wgt1  = RcppDeepState_NumericMatrix();
  wgt1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_freq/inputs/wgt1");
  wgt1_stream << wgt1;
  std::cout << "wgt1 values: "<< wgt1 << std::endl;
  wgt1_stream.close();
  NumericMatrix wgtrep  = RcppDeepState_NumericMatrix();
  wgtrep_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_freq/inputs/wgtrep");
  wgtrep_stream << wgtrep;
  std::cout << "wgtrep values: "<< wgtrep << std::endl;
  wgtrep_stream.close();
  NumericVector vars_index  = RcppDeepState_NumericVector();
  vars_index_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_freq/inputs/vars_index");
  vars_index_stream << vars_index;
  std::cout << "vars_index values: "<< vars_index << std::endl;
  vars_index_stream.close();
  NumericVector fayfac  = RcppDeepState_NumericVector();
  fayfac_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_freq/inputs/fayfac");
  fayfac_stream << fayfac;
  std::cout << "fayfac values: "<< fayfac << std::endl;
  fayfac_stream.close();
  NumericVector NI  = RcppDeepState_NumericVector();
  NI_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_freq/inputs/NI");
  NI_stream << NI;
  std::cout << "NI values: "<< NI << std::endl;
  NI_stream.close();
  NumericVector group_index1  = RcppDeepState_NumericVector();
  group_index1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_freq/inputs/group_index1");
  group_index1_stream << group_index1;
  std::cout << "group_index1 values: "<< group_index1 << std::endl;
  group_index1_stream.close();
  NumericVector group_values  = RcppDeepState_NumericVector();
  group_values_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_freq/inputs/group_values");
  group_values_stream << group_values;
  std::cout << "group_values values: "<< group_values << std::endl;
  group_values_stream.close();
  NumericMatrix vars_values  = RcppDeepState_NumericMatrix();
  vars_values_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_freq/inputs/vars_values");
  vars_values_stream << vars_values;
  std::cout << "vars_values values: "<< vars_values << std::endl;
  vars_values_stream.close();
  NumericVector vars_values_numb  = RcppDeepState_NumericVector();
  vars_values_numb_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_freq/inputs/vars_values_numb");
  vars_values_numb_stream << vars_values_numb;
  std::cout << "vars_values_numb values: "<< vars_values_numb << std::endl;
  vars_values_numb_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifie_freq(datalist,wgt1,wgtrep,vars_index,fayfac,NI,group_index1,group_values,vars_values,vars_values_numb);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
