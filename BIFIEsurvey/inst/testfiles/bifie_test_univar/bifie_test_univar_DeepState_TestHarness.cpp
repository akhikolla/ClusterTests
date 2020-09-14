#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::List bifie_test_univar(Rcpp::NumericMatrix mean1M, Rcpp::NumericMatrix sd1M, Rcpp::NumericMatrix sumweightM, int GG, Rcpp::NumericVector group_values, Rcpp::NumericMatrix mean1repM, Rcpp::NumericMatrix sd1repM, Rcpp::NumericMatrix sumweightrepM, Rcpp::NumericVector fayfac);

TEST(BIFIEsurvey_deepstate_test,bifie_test_univar_test){
  std::ofstream mean1M_stream;
  std::ofstream sd1M_stream;
  std::ofstream sumweightM_stream;
  std::ofstream GG_stream;
  std::ofstream group_values_stream;
  std::ofstream mean1repM_stream;
  std::ofstream sd1repM_stream;
  std::ofstream sumweightrepM_stream;
  std::ofstream fayfac_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix mean1M  = RcppDeepState_NumericMatrix();
  mean1M_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_test_univar/inputs/mean1M");
  mean1M_stream << mean1M;
  std::cout << "mean1M values: "<< mean1M << std::endl;
  mean1M_stream.close();
  NumericMatrix sd1M  = RcppDeepState_NumericMatrix();
  sd1M_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_test_univar/inputs/sd1M");
  sd1M_stream << sd1M;
  std::cout << "sd1M values: "<< sd1M << std::endl;
  sd1M_stream.close();
  NumericMatrix sumweightM  = RcppDeepState_NumericMatrix();
  sumweightM_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_test_univar/inputs/sumweightM");
  sumweightM_stream << sumweightM;
  std::cout << "sumweightM values: "<< sumweightM << std::endl;
  sumweightM_stream.close();
  int GG  = RcppDeepState_int();
  GG_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_test_univar/inputs/GG");
  GG_stream << GG;
  std::cout << "GG values: "<< GG << std::endl;
  GG_stream.close();
  NumericVector group_values  = RcppDeepState_NumericVector();
  group_values_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_test_univar/inputs/group_values");
  group_values_stream << group_values;
  std::cout << "group_values values: "<< group_values << std::endl;
  group_values_stream.close();
  NumericMatrix mean1repM  = RcppDeepState_NumericMatrix();
  mean1repM_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_test_univar/inputs/mean1repM");
  mean1repM_stream << mean1repM;
  std::cout << "mean1repM values: "<< mean1repM << std::endl;
  mean1repM_stream.close();
  NumericMatrix sd1repM  = RcppDeepState_NumericMatrix();
  sd1repM_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_test_univar/inputs/sd1repM");
  sd1repM_stream << sd1repM;
  std::cout << "sd1repM values: "<< sd1repM << std::endl;
  sd1repM_stream.close();
  NumericMatrix sumweightrepM  = RcppDeepState_NumericMatrix();
  sumweightrepM_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_test_univar/inputs/sumweightrepM");
  sumweightrepM_stream << sumweightrepM;
  std::cout << "sumweightrepM values: "<< sumweightrepM << std::endl;
  sumweightrepM_stream.close();
  NumericVector fayfac  = RcppDeepState_NumericVector();
  fayfac_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifie_test_univar/inputs/fayfac");
  fayfac_stream << fayfac;
  std::cout << "fayfac values: "<< fayfac << std::endl;
  fayfac_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifie_test_univar(mean1M,sd1M,sumweightM,GG,group_values,mean1repM,sd1repM,sumweightrepM,fayfac);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
