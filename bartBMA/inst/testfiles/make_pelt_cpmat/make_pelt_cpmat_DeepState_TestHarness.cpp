#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

List make_pelt_cpmat(NumericMatrix data, NumericVector resp, double pen, int num_cp);

TEST(bartBMA_deepstate_test,make_pelt_cpmat_test){
  std::ofstream data_stream;
  std::ofstream resp_stream;
  std::ofstream pen_stream;
  std::ofstream num_cp_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix data  = RcppDeepState_NumericMatrix();
  data_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/make_pelt_cpmat/inputs/data");
  data_stream << data;
  std::cout << "data values: "<< data << std::endl;
  data_stream.close();
  NumericVector resp  = RcppDeepState_NumericVector();
  resp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/make_pelt_cpmat/inputs/resp");
  resp_stream << resp;
  std::cout << "resp values: "<< resp << std::endl;
  resp_stream.close();
  double pen  = RcppDeepState_double();
  pen_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/make_pelt_cpmat/inputs/pen");
  pen_stream << pen;
  std::cout << "pen values: "<< pen << std::endl;
  pen_stream.close();
  int num_cp  = RcppDeepState_int();
  num_cp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/make_pelt_cpmat/inputs/num_cp");
  num_cp_stream << num_cp;
  std::cout << "num_cp values: "<< num_cp << std::endl;
  num_cp_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    make_pelt_cpmat(data,resp,pen,num_cp);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
