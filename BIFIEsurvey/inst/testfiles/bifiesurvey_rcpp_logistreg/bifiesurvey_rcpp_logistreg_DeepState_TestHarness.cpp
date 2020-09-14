#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifiesurvey_rcpp_logistreg(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector dep_index, Rcpp::NumericVector pre_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values, double eps, int maxiter);

TEST(BIFIEsurvey_deepstate_test,bifiesurvey_rcpp_logistreg_test){
  std::ofstream datalist_stream;
  std::ofstream wgt1_stream;
  std::ofstream wgtrep_stream;
  std::ofstream dep_index_stream;
  std::ofstream pre_index_stream;
  std::ofstream fayfac_stream;
  std::ofstream NI_stream;
  std::ofstream group_index1_stream;
  std::ofstream group_values_stream;
  std::ofstream eps_stream;
  std::ofstream maxiter_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix datalist  = RcppDeepState_NumericMatrix();
  datalist_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_logistreg/inputs/datalist");
  datalist_stream << datalist;
  std::cout << "datalist values: "<< datalist << std::endl;
  datalist_stream.close();
  NumericMatrix wgt1  = RcppDeepState_NumericMatrix();
  wgt1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_logistreg/inputs/wgt1");
  wgt1_stream << wgt1;
  std::cout << "wgt1 values: "<< wgt1 << std::endl;
  wgt1_stream.close();
  NumericMatrix wgtrep  = RcppDeepState_NumericMatrix();
  wgtrep_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_logistreg/inputs/wgtrep");
  wgtrep_stream << wgtrep;
  std::cout << "wgtrep values: "<< wgtrep << std::endl;
  wgtrep_stream.close();
  NumericVector dep_index  = RcppDeepState_NumericVector();
  dep_index_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_logistreg/inputs/dep_index");
  dep_index_stream << dep_index;
  std::cout << "dep_index values: "<< dep_index << std::endl;
  dep_index_stream.close();
  NumericVector pre_index  = RcppDeepState_NumericVector();
  pre_index_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_logistreg/inputs/pre_index");
  pre_index_stream << pre_index;
  std::cout << "pre_index values: "<< pre_index << std::endl;
  pre_index_stream.close();
  NumericVector fayfac  = RcppDeepState_NumericVector();
  fayfac_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_logistreg/inputs/fayfac");
  fayfac_stream << fayfac;
  std::cout << "fayfac values: "<< fayfac << std::endl;
  fayfac_stream.close();
  NumericVector NI  = RcppDeepState_NumericVector();
  NI_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_logistreg/inputs/NI");
  NI_stream << NI;
  std::cout << "NI values: "<< NI << std::endl;
  NI_stream.close();
  NumericVector group_index1  = RcppDeepState_NumericVector();
  group_index1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_logistreg/inputs/group_index1");
  group_index1_stream << group_index1;
  std::cout << "group_index1 values: "<< group_index1 << std::endl;
  group_index1_stream.close();
  NumericVector group_values  = RcppDeepState_NumericVector();
  group_values_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_logistreg/inputs/group_values");
  group_values_stream << group_values;
  std::cout << "group_values values: "<< group_values << std::endl;
  group_values_stream.close();
  double eps  = RcppDeepState_double();
  eps_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_logistreg/inputs/eps");
  eps_stream << eps;
  std::cout << "eps values: "<< eps << std::endl;
  eps_stream.close();
  int maxiter  = RcppDeepState_int();
  maxiter_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_logistreg/inputs/maxiter");
  maxiter_stream << maxiter;
  std::cout << "maxiter values: "<< maxiter << std::endl;
  maxiter_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifiesurvey_rcpp_logistreg(datalist,wgt1,wgtrep,dep_index,pre_index,fayfac,NI,group_index1,group_values,eps,maxiter);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
