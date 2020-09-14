#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::mat Sigma_i_not_i(arma::mat x, int index);

TEST(BGGM_deepstate_test,Sigma_i_not_i_test){
  std::ofstream x_stream;
  std::ofstream index_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat x  = RcppDeepState_mat();
  x_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BGGM/inst/testfiles/Sigma_i_not_i/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  int index  = RcppDeepState_int();
  index_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BGGM/inst/testfiles/Sigma_i_not_i/inputs/index");
  index_stream << index;
  std::cout << "index values: "<< index << std::endl;
  index_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    Sigma_i_not_i(x,index);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
