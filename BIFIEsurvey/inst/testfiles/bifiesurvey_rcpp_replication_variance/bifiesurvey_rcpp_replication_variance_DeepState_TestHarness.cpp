#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

Rcpp::NumericVector bifiesurvey_rcpp_replication_variance(Rcpp::NumericVector pars, Rcpp::NumericMatrix pars_repl, Rcpp::NumericVector fay_factor);

TEST(BIFIEsurvey_deepstate_test,bifiesurvey_rcpp_replication_variance_test){
  std::ofstream pars_stream;
  std::ofstream pars_repl_stream;
  std::ofstream fay_factor_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector pars  = RcppDeepState_NumericVector();
  pars_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_replication_variance/inputs/pars");
  pars_stream << pars;
  std::cout << "pars values: "<< pars << std::endl;
  pars_stream.close();
  NumericMatrix pars_repl  = RcppDeepState_NumericMatrix();
  pars_repl_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_replication_variance/inputs/pars_repl");
  pars_repl_stream << pars_repl;
  std::cout << "pars_repl values: "<< pars_repl << std::endl;
  pars_repl_stream.close();
  NumericVector fay_factor  = RcppDeepState_NumericVector();
  fay_factor_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BIFIEsurvey/inst/testfiles/bifiesurvey_rcpp_replication_variance/inputs/fay_factor");
  fay_factor_stream << fay_factor;
  std::cout << "fay_factor values: "<< fay_factor << std::endl;
  fay_factor_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    bifiesurvey_rcpp_replication_variance(pars,pars_repl,fay_factor);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
