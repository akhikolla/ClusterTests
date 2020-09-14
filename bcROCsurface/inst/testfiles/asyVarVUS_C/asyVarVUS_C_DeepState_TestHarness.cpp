#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector asyVarVUS_C(NumericVector tt, NumericMatrix D_hat, double mu_hat, NumericMatrix EstFunc, NumericMatrix Hess_inv, NumericMatrix Der_D1_hat, NumericMatrix Der_D2_hat, NumericMatrix Der_D3_hat);

TEST(bcROCsurface_deepstate_test,asyVarVUS_C_test){
  std::ofstream tt_stream;
  std::ofstream D_hat_stream;
  std::ofstream mu_hat_stream;
  std::ofstream EstFunc_stream;
  std::ofstream Hess_inv_stream;
  std::ofstream Der_D1_hat_stream;
  std::ofstream Der_D2_hat_stream;
  std::ofstream Der_D3_hat_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector tt  = RcppDeepState_NumericVector();
  tt_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/asyVarVUS_C/inputs/tt");
  tt_stream << tt;
  std::cout << "tt values: "<< tt << std::endl;
  tt_stream.close();
  NumericMatrix D_hat  = RcppDeepState_NumericMatrix();
  D_hat_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/asyVarVUS_C/inputs/D_hat");
  D_hat_stream << D_hat;
  std::cout << "D_hat values: "<< D_hat << std::endl;
  D_hat_stream.close();
  double mu_hat  = RcppDeepState_double();
  mu_hat_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/asyVarVUS_C/inputs/mu_hat");
  mu_hat_stream << mu_hat;
  std::cout << "mu_hat values: "<< mu_hat << std::endl;
  mu_hat_stream.close();
  NumericMatrix EstFunc  = RcppDeepState_NumericMatrix();
  EstFunc_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/asyVarVUS_C/inputs/EstFunc");
  EstFunc_stream << EstFunc;
  std::cout << "EstFunc values: "<< EstFunc << std::endl;
  EstFunc_stream.close();
  NumericMatrix Hess_inv  = RcppDeepState_NumericMatrix();
  Hess_inv_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/asyVarVUS_C/inputs/Hess_inv");
  Hess_inv_stream << Hess_inv;
  std::cout << "Hess_inv values: "<< Hess_inv << std::endl;
  Hess_inv_stream.close();
  NumericMatrix Der_D1_hat  = RcppDeepState_NumericMatrix();
  Der_D1_hat_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/asyVarVUS_C/inputs/Der_D1_hat");
  Der_D1_hat_stream << Der_D1_hat;
  std::cout << "Der_D1_hat values: "<< Der_D1_hat << std::endl;
  Der_D1_hat_stream.close();
  NumericMatrix Der_D2_hat  = RcppDeepState_NumericMatrix();
  Der_D2_hat_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/asyVarVUS_C/inputs/Der_D2_hat");
  Der_D2_hat_stream << Der_D2_hat;
  std::cout << "Der_D2_hat values: "<< Der_D2_hat << std::endl;
  Der_D2_hat_stream.close();
  NumericMatrix Der_D3_hat  = RcppDeepState_NumericMatrix();
  Der_D3_hat_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/asyVarVUS_C/inputs/Der_D3_hat");
  Der_D3_hat_stream << Der_D3_hat;
  std::cout << "Der_D3_hat values: "<< Der_D3_hat << std::endl;
  Der_D3_hat_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    asyVarVUS_C(tt,D_hat,mu_hat,EstFunc,Hess_inv,Der_D1_hat,Der_D2_hat,Der_D3_hat);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
