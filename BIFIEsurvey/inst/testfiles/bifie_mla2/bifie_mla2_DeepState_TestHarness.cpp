#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifie_mla2(Rcpp::NumericMatrix X_list, Rcpp::NumericMatrix Z_list, Rcpp::NumericVector y_list, Rcpp::NumericVector wgttot, Rcpp::NumericVector wgtlev2, Rcpp::NumericVector wgtlev1, double globconv, int maxiter, Rcpp::NumericVector group, Rcpp::NumericVector group_values, Rcpp::NumericVector cluster, Rcpp::NumericMatrix wgtrep, int Nimp, Rcpp::NumericVector fayfac, Rcpp::NumericMatrix recov_constraint, int is_rcov_constraint);

TEST(BIFIEsurvey_deepstate_test,bifie_mla2_test){
  std::ofstream X_list_stream;
  std::ofstream Z_list_stream;
  std::ofstream y_list_stream;
  std::ofstream wgttot_stream;
  std::ofstream wgtlev2_stream;
  std::ofstream wgtlev1_stream;
  std::ofstream globconv_stream;
  std::ofstream maxiter_stream;
  std::ofstream group_stream;
  std::ofstream group_values_stream;
  std::ofstream cluster_stream;
  std::ofstream wgtrep_stream;
  std::ofstream Nimp_stream;
  std::ofstream fayfac_stream;
  std::ofstream recov_constraint_stream;
  std::ofstream is_rcov_constraint_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix X_list  = RcppDeepState_NumericMatrix();
  X_list_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/X_list");
  X_list_stream << X_list;
  std::cout << "X_list values: "<< X_list << std::endl;
  X_list_stream.close();
  NumericMatrix Z_list  = RcppDeepState_NumericMatrix();
  Z_list_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/Z_list");
  Z_list_stream << Z_list;
  std::cout << "Z_list values: "<< Z_list << std::endl;
  Z_list_stream.close();
  NumericVector y_list  = RcppDeepState_NumericVector();
  y_list_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/y_list");
  y_list_stream << y_list;
  std::cout << "y_list values: "<< y_list << std::endl;
  y_list_stream.close();
  NumericVector wgttot  = RcppDeepState_NumericVector();
  wgttot_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/wgttot");
  wgttot_stream << wgttot;
  std::cout << "wgttot values: "<< wgttot << std::endl;
  wgttot_stream.close();
  NumericVector wgtlev2  = RcppDeepState_NumericVector();
  wgtlev2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/wgtlev2");
  wgtlev2_stream << wgtlev2;
  std::cout << "wgtlev2 values: "<< wgtlev2 << std::endl;
  wgtlev2_stream.close();
  NumericVector wgtlev1  = RcppDeepState_NumericVector();
  wgtlev1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/wgtlev1");
  wgtlev1_stream << wgtlev1;
  std::cout << "wgtlev1 values: "<< wgtlev1 << std::endl;
  wgtlev1_stream.close();
  double globconv  = RcppDeepState_double();
  globconv_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/globconv");
  globconv_stream << globconv;
  std::cout << "globconv values: "<< globconv << std::endl;
  globconv_stream.close();
  int maxiter  = RcppDeepState_int();
  maxiter_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/maxiter");
  maxiter_stream << maxiter;
  std::cout << "maxiter values: "<< maxiter << std::endl;
  maxiter_stream.close();
  NumericVector group  = RcppDeepState_NumericVector();
  group_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/group");
  group_stream << group;
  std::cout << "group values: "<< group << std::endl;
  group_stream.close();
  NumericVector group_values  = RcppDeepState_NumericVector();
  group_values_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/group_values");
  group_values_stream << group_values;
  std::cout << "group_values values: "<< group_values << std::endl;
  group_values_stream.close();
  NumericVector cluster  = RcppDeepState_NumericVector();
  cluster_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/cluster");
  cluster_stream << cluster;
  std::cout << "cluster values: "<< cluster << std::endl;
  cluster_stream.close();
  NumericMatrix wgtrep  = RcppDeepState_NumericMatrix();
  wgtrep_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/wgtrep");
  wgtrep_stream << wgtrep;
  std::cout << "wgtrep values: "<< wgtrep << std::endl;
  wgtrep_stream.close();
  int Nimp  = RcppDeepState_int();
  Nimp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/Nimp");
  Nimp_stream << Nimp;
  std::cout << "Nimp values: "<< Nimp << std::endl;
  Nimp_stream.close();
  NumericVector fayfac  = RcppDeepState_NumericVector();
  fayfac_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/fayfac");
  fayfac_stream << fayfac;
  std::cout << "fayfac values: "<< fayfac << std::endl;
  fayfac_stream.close();
  NumericMatrix recov_raint  = RcppDeepState_NumericMatrix();
  recov_constraint_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/recov_constraint");
  recov_constraint_stream << recov_constraint;
  std::cout << "recov_constraint values: "<< recov_constraint << std::endl;
  recov_constraint_stream.close();
  int is_rcov_raint  = RcppDeepState_int();
  is_rcov_constraint_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_mla2/inputs/is_rcov_constraint");
  is_rcov_constraint_stream << is_rcov_constraint;
  std::cout << "is_rcov_constraint values: "<< is_rcov_constraint << std::endl;
  is_rcov_constraint_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifie_mla2(X_list,Z_list,y_list,wgttot,wgtlev2,wgtlev1,globconv,maxiter,group,group_values,cluster,wgtrep,Nimp,fayfac,recov_constraint,is_rcov_constraint);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
