#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::mat acvBlockMatrix(arma::mat acv);

TEST(beyondWhittle_deepstate_test,acvBlockMatrix_test){
  std::ofstream acv_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat acv  = RcppDeepState_mat();
  acv_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/acvBlockMatrix/inputs/acv");
  acv_stream << acv;
  std::cout << "acv values: "<< acv << std::endl;
  acv_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    acvBlockMatrix(acv);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
