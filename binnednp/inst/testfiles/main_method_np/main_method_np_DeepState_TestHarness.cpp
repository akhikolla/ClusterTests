#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

void main_method_np(int hn, int B, NumericVector hseq, NumericMatrix xbm, NumericVector y, NumericVector t, NumericVector MSE_J1, NumericVector MSE_J2, double J1_np, double J2_np, NumericVector J1b, NumericVector J2b);

TEST(binnednp_deepstate_test,main_method_np_test){
  std::ofstream hn_stream;
  std::ofstream B_stream;
  std::ofstream hseq_stream;
  std::ofstream xbm_stream;
  std::ofstream y_stream;
  std::ofstream t_stream;
  std::ofstream MSE_J1_stream;
  std::ofstream MSE_J2_stream;
  std::ofstream J1_np_stream;
  std::ofstream J2_np_stream;
  std::ofstream J1b_stream;
  std::ofstream J2b_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int hn  = RcppDeepState_int();
  hn_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/main_method_np/inputs/hn");
  hn_stream << hn;
  std::cout << "hn values: "<< hn << std::endl;
  hn_stream.close();
  int B  = RcppDeepState_int();
  B_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/main_method_np/inputs/B");
  B_stream << B;
  std::cout << "B values: "<< B << std::endl;
  B_stream.close();
  NumericVector hseq  = RcppDeepState_NumericVector();
  hseq_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/main_method_np/inputs/hseq");
  hseq_stream << hseq;
  std::cout << "hseq values: "<< hseq << std::endl;
  hseq_stream.close();
  NumericMatrix xbm  = RcppDeepState_NumericMatrix();
  xbm_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/main_method_np/inputs/xbm");
  xbm_stream << xbm;
  std::cout << "xbm values: "<< xbm << std::endl;
  xbm_stream.close();
  NumericVector y  = RcppDeepState_NumericVector();
  y_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/main_method_np/inputs/y");
  y_stream << y;
  std::cout << "y values: "<< y << std::endl;
  y_stream.close();
  NumericVector t  = RcppDeepState_NumericVector();
  t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/main_method_np/inputs/t");
  t_stream << t;
  std::cout << "t values: "<< t << std::endl;
  t_stream.close();
  NumericVector MSE_J1  = RcppDeepState_NumericVector();
  MSE_J1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/main_method_np/inputs/MSE_J1");
  MSE_J1_stream << MSE_J1;
  std::cout << "MSE_J1 values: "<< MSE_J1 << std::endl;
  MSE_J1_stream.close();
  NumericVector MSE_J2  = RcppDeepState_NumericVector();
  MSE_J2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/main_method_np/inputs/MSE_J2");
  MSE_J2_stream << MSE_J2;
  std::cout << "MSE_J2 values: "<< MSE_J2 << std::endl;
  MSE_J2_stream.close();
  double J1_np  = RcppDeepState_double();
  J1_np_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/main_method_np/inputs/J1_np");
  J1_np_stream << J1_np;
  std::cout << "J1_np values: "<< J1_np << std::endl;
  J1_np_stream.close();
  double J2_np  = RcppDeepState_double();
  J2_np_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/main_method_np/inputs/J2_np");
  J2_np_stream << J2_np;
  std::cout << "J2_np values: "<< J2_np << std::endl;
  J2_np_stream.close();
  NumericVector J1b  = RcppDeepState_NumericVector();
  J1b_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/main_method_np/inputs/J1b");
  J1b_stream << J1b;
  std::cout << "J1b values: "<< J1b << std::endl;
  J1b_stream.close();
  NumericVector J2b  = RcppDeepState_NumericVector();
  J2b_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/main_method_np/inputs/J2b");
  J2b_stream << J2b;
  std::cout << "J2b values: "<< J2b << std::endl;
  J2b_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    main_method_np(hn,B,hseq,xbm,y,t,MSE_J1,MSE_J2,J1_np,J2_np,J1b,J2b);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
