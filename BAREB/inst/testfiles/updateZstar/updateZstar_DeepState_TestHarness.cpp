#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::mat updateZstar(arma::mat mu_star, arma::mat delta, int n, int T0);

TEST(BAREB_deepstate_test,updateZstar_test){
  std::ofstream mu_star_stream;
  std::ofstream delta_stream;
  std::ofstream n_stream;
  std::ofstream T0_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat mu_star  = RcppDeepState_mat();
  mu_star_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updateZstar/inputs/mu_star");
  mu_star_stream << mu_star;
  std::cout << "mu_star values: "<< mu_star << std::endl;
  mu_star_stream.close();
  arma::mat delta  = RcppDeepState_mat();
  delta_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updateZstar/inputs/delta");
  delta_stream << delta;
  std::cout << "delta values: "<< delta << std::endl;
  delta_stream.close();
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updateZstar/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  int T0  = RcppDeepState_int();
  T0_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updateZstar/inputs/T0");
  T0_stream << T0;
  std::cout << "T0 values: "<< T0 << std::endl;
  T0_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    updateZstar(mu_star,delta,n,T0);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
