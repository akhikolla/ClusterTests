#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::uvec myfind(IntegerVector evec, int e);

TEST(BCSub_deepstate_test,myfind_test){
  std::ofstream evec_stream;
  std::ofstream e_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  IntegerVector evec  = RcppDeepState_IntegerVector();
  evec_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BCSub/inst/testfiles/myfind/inputs/evec");
  evec_stream << evec;
  std::cout << "evec values: "<< evec << std::endl;
  evec_stream.close();
  int e  = RcppDeepState_int();
  e_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BCSub/inst/testfiles/myfind/inputs/e");
  e_stream << e;
  std::cout << "e values: "<< e << std::endl;
  e_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    myfind(evec,e);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
