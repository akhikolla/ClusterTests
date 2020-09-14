#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::rowvec updatec(arma::mat Zstar, arma::mat mu, arma::mat D, double sigmac, double sigma_square, int n, int T0);

TEST(BAREB_deepstate_test,updatec_test){
  std::ofstream Zstar_stream;
  std::ofstream mu_stream;
  std::ofstream D_stream;
  std::ofstream sigmac_stream;
  std::ofstream sigma_square_stream;
  std::ofstream n_stream;
  std::ofstream T0_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat Zstar  = RcppDeepState_mat();
  Zstar_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updatec/inputs/Zstar");
  Zstar_stream << Zstar;
  std::cout << "Zstar values: "<< Zstar << std::endl;
  Zstar_stream.close();
  arma::mat mu  = RcppDeepState_mat();
  mu_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updatec/inputs/mu");
  mu_stream << mu;
  std::cout << "mu values: "<< mu << std::endl;
  mu_stream.close();
  arma::mat D  = RcppDeepState_mat();
  D_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updatec/inputs/D");
  D_stream << D;
  std::cout << "D values: "<< D << std::endl;
  D_stream.close();
  double sigmac  = RcppDeepState_double();
  sigmac_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updatec/inputs/sigmac");
  sigmac_stream << sigmac;
  std::cout << "sigmac values: "<< sigmac << std::endl;
  sigmac_stream.close();
  double sigma_square  = RcppDeepState_double();
  sigma_square_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updatec/inputs/sigma_square");
  sigma_square_stream << sigma_square;
  std::cout << "sigma_square values: "<< sigma_square << std::endl;
  sigma_square_stream.close();
  int n  = RcppDeepState_int();
  n_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updatec/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  int T0  = RcppDeepState_int();
  T0_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BAREB/inst/testfiles/updatec/inputs/T0");
  T0_stream << T0;
  std::cout << "T0 values: "<< T0 << std::endl;
  T0_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    updatec(Zstar,mu,D,sigmac,sigma_square,n,T0);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
