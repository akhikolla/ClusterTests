#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

arma::mat r_maternCov(arma::mat dist, double scale, double range, double smoothness, double nugget);

TEST(bisque_deepstate_test,r_maternCov_test){
  std::ofstream dist_stream;
  std::ofstream scale_stream;
  std::ofstream range_stream;
  std::ofstream smoothness_stream;
  std::ofstream nugget_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  arma::mat dist  = RcppDeepState_mat();
  dist_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bisque/inst/testfiles/r_maternCov/inputs/dist");
  dist_stream << dist;
  std::cout << "dist values: "<< dist << std::endl;
  dist_stream.close();
  double scale  = RcppDeepState_double();
  scale_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bisque/inst/testfiles/r_maternCov/inputs/scale");
  scale_stream << scale;
  std::cout << "scale values: "<< scale << std::endl;
  scale_stream.close();
  double range  = RcppDeepState_double();
  range_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bisque/inst/testfiles/r_maternCov/inputs/range");
  range_stream << range;
  std::cout << "range values: "<< range << std::endl;
  range_stream.close();
  double smoothness  = RcppDeepState_double();
  smoothness_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bisque/inst/testfiles/r_maternCov/inputs/smoothness");
  smoothness_stream << smoothness;
  std::cout << "smoothness values: "<< smoothness << std::endl;
  smoothness_stream.close();
  double nugget  = RcppDeepState_double();
  nugget_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bisque/inst/testfiles/r_maternCov/inputs/nugget");
  nugget_stream << nugget;
  std::cout << "nugget values: "<< nugget << std::endl;
  nugget_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    r_maternCov(dist,scale,range,smoothness,nugget);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
