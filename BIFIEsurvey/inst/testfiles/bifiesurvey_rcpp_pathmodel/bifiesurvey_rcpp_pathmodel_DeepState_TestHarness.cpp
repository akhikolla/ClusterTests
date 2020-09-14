#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifiesurvey_rcpp_pathmodel(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values, Rcpp::NumericMatrix L, Rcpp::NumericVector L_row_index, int NL, Rcpp::NumericMatrix E, Rcpp::NumericMatrix R, Rcpp::NumericVector R_row_index, Rcpp::NumericMatrix coeff_index, int NP0, Rcpp::NumericVector unreliability);

TEST(BIFIEsurvey_deepstate_test,bifiesurvey_rcpp_pathmodel_test){
  std::ofstream datalist_stream;
  std::ofstream wgt1_stream;
  std::ofstream wgtrep_stream;
  std::ofstream vars_index_stream;
  std::ofstream fayfac_stream;
  std::ofstream NI_stream;
  std::ofstream group_index1_stream;
  std::ofstream group_values_stream;
  std::ofstream L_stream;
  std::ofstream L_row_index_stream;
  std::ofstream NL_stream;
  std::ofstream E_stream;
  std::ofstream R_stream;
  std::ofstream R_row_index_stream;
  std::ofstream coeff_index_stream;
  std::ofstream NP0_stream;
  std::ofstream unreliability_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix datalist  = RcppDeepState_NumericMatrix();
  datalist_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/datalist");
  datalist_stream << datalist;
  std::cout << "datalist values: "<< datalist << std::endl;
  datalist_stream.close();
  NumericMatrix wgt1  = RcppDeepState_NumericMatrix();
  wgt1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/wgt1");
  wgt1_stream << wgt1;
  std::cout << "wgt1 values: "<< wgt1 << std::endl;
  wgt1_stream.close();
  NumericMatrix wgtrep  = RcppDeepState_NumericMatrix();
  wgtrep_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/wgtrep");
  wgtrep_stream << wgtrep;
  std::cout << "wgtrep values: "<< wgtrep << std::endl;
  wgtrep_stream.close();
  NumericVector vars_index  = RcppDeepState_NumericVector();
  vars_index_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/vars_index");
  vars_index_stream << vars_index;
  std::cout << "vars_index values: "<< vars_index << std::endl;
  vars_index_stream.close();
  NumericVector fayfac  = RcppDeepState_NumericVector();
  fayfac_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/fayfac");
  fayfac_stream << fayfac;
  std::cout << "fayfac values: "<< fayfac << std::endl;
  fayfac_stream.close();
  NumericVector NI  = RcppDeepState_NumericVector();
  NI_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/NI");
  NI_stream << NI;
  std::cout << "NI values: "<< NI << std::endl;
  NI_stream.close();
  NumericVector group_index1  = RcppDeepState_NumericVector();
  group_index1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/group_index1");
  group_index1_stream << group_index1;
  std::cout << "group_index1 values: "<< group_index1 << std::endl;
  group_index1_stream.close();
  NumericVector group_values  = RcppDeepState_NumericVector();
  group_values_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/group_values");
  group_values_stream << group_values;
  std::cout << "group_values values: "<< group_values << std::endl;
  group_values_stream.close();
  NumericMatrix L  = RcppDeepState_NumericMatrix();
  L_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/L");
  L_stream << L;
  std::cout << "L values: "<< L << std::endl;
  L_stream.close();
  NumericVector L_row_index  = RcppDeepState_NumericVector();
  L_row_index_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/L_row_index");
  L_row_index_stream << L_row_index;
  std::cout << "L_row_index values: "<< L_row_index << std::endl;
  L_row_index_stream.close();
  int NL  = RcppDeepState_int();
  NL_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/NL");
  NL_stream << NL;
  std::cout << "NL values: "<< NL << std::endl;
  NL_stream.close();
  NumericMatrix E  = RcppDeepState_NumericMatrix();
  E_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/E");
  E_stream << E;
  std::cout << "E values: "<< E << std::endl;
  E_stream.close();
  NumericMatrix R  = RcppDeepState_NumericMatrix();
  R_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/R");
  R_stream << R;
  std::cout << "R values: "<< R << std::endl;
  R_stream.close();
  NumericVector R_row_index  = RcppDeepState_NumericVector();
  R_row_index_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/R_row_index");
  R_row_index_stream << R_row_index;
  std::cout << "R_row_index values: "<< R_row_index << std::endl;
  R_row_index_stream.close();
  NumericMatrix coeff_index  = RcppDeepState_NumericMatrix();
  coeff_index_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/coeff_index");
  coeff_index_stream << coeff_index;
  std::cout << "coeff_index values: "<< coeff_index << std::endl;
  coeff_index_stream.close();
  int NP0  = RcppDeepState_int();
  NP0_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/NP0");
  NP0_stream << NP0;
  std::cout << "NP0 values: "<< NP0 << std::endl;
  NP0_stream.close();
  NumericVector unreliability  = RcppDeepState_NumericVector();
  unreliability_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_pathmodel/inputs/unreliability");
  unreliability_stream << unreliability;
  std::cout << "unreliability values: "<< unreliability << std::endl;
  unreliability_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifiesurvey_rcpp_pathmodel(datalist,wgt1,wgtrep,vars_index,fayfac,NI,group_index1,group_values,L,L_row_index,NL,E,R,R_row_index,coeff_index,NP0,unreliability);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
