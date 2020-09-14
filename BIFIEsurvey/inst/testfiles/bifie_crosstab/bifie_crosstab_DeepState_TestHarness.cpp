#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifie_crosstab(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_values1, Rcpp::NumericVector vars_index1, Rcpp::NumericVector vars_values2, Rcpp::NumericVector vars_index2, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values);

TEST(BIFIEsurvey_deepstate_test,bifie_crosstab_test){
  std::ofstream datalist_stream;
  std::ofstream wgt1_stream;
  std::ofstream wgtrep_stream;
  std::ofstream vars_values1_stream;
  std::ofstream vars_index1_stream;
  std::ofstream vars_values2_stream;
  std::ofstream vars_index2_stream;
  std::ofstream fayfac_stream;
  std::ofstream NI_stream;
  std::ofstream group_index1_stream;
  std::ofstream group_values_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix datalist  = RcppDeepState_NumericMatrix();
  datalist_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_crosstab/inputs/datalist");
  datalist_stream << datalist;
  std::cout << "datalist values: "<< datalist << std::endl;
  datalist_stream.close();
  NumericMatrix wgt1  = RcppDeepState_NumericMatrix();
  wgt1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_crosstab/inputs/wgt1");
  wgt1_stream << wgt1;
  std::cout << "wgt1 values: "<< wgt1 << std::endl;
  wgt1_stream.close();
  NumericMatrix wgtrep  = RcppDeepState_NumericMatrix();
  wgtrep_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_crosstab/inputs/wgtrep");
  wgtrep_stream << wgtrep;
  std::cout << "wgtrep values: "<< wgtrep << std::endl;
  wgtrep_stream.close();
  NumericVector vars_values1  = RcppDeepState_NumericVector();
  vars_values1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_crosstab/inputs/vars_values1");
  vars_values1_stream << vars_values1;
  std::cout << "vars_values1 values: "<< vars_values1 << std::endl;
  vars_values1_stream.close();
  NumericVector vars_index1  = RcppDeepState_NumericVector();
  vars_index1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_crosstab/inputs/vars_index1");
  vars_index1_stream << vars_index1;
  std::cout << "vars_index1 values: "<< vars_index1 << std::endl;
  vars_index1_stream.close();
  NumericVector vars_values2  = RcppDeepState_NumericVector();
  vars_values2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_crosstab/inputs/vars_values2");
  vars_values2_stream << vars_values2;
  std::cout << "vars_values2 values: "<< vars_values2 << std::endl;
  vars_values2_stream.close();
  NumericVector vars_index2  = RcppDeepState_NumericVector();
  vars_index2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_crosstab/inputs/vars_index2");
  vars_index2_stream << vars_index2;
  std::cout << "vars_index2 values: "<< vars_index2 << std::endl;
  vars_index2_stream.close();
  NumericVector fayfac  = RcppDeepState_NumericVector();
  fayfac_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_crosstab/inputs/fayfac");
  fayfac_stream << fayfac;
  std::cout << "fayfac values: "<< fayfac << std::endl;
  fayfac_stream.close();
  NumericVector NI  = RcppDeepState_NumericVector();
  NI_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_crosstab/inputs/NI");
  NI_stream << NI;
  std::cout << "NI values: "<< NI << std::endl;
  NI_stream.close();
  NumericVector group_index1  = RcppDeepState_NumericVector();
  group_index1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_crosstab/inputs/group_index1");
  group_index1_stream << group_index1;
  std::cout << "group_index1 values: "<< group_index1 << std::endl;
  group_index1_stream.close();
  NumericVector group_values  = RcppDeepState_NumericVector();
  group_values_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_crosstab/inputs/group_values");
  group_values_stream << group_values;
  std::cout << "group_values values: "<< group_values << std::endl;
  group_values_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifie_crosstab(datalist,wgt1,wgtrep,vars_values1,vars_index1,vars_values2,vars_index2,fayfac,NI,group_index1,group_values);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
