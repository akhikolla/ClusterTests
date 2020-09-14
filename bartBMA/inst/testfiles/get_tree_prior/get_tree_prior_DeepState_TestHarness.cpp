#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double get_tree_prior(double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson, NumericMatrix tree_table, NumericMatrix tree_matrix, double alpha, double beta);

TEST(bartBMA_deepstate_test,get_tree_prior_test){
  std::ofstream spike_tree_stream;
  std::ofstream s_t_hyperprior_stream;
  std::ofstream p_s_t_stream;
  std::ofstream a_s_t_stream;
  std::ofstream b_s_t_stream;
  std::ofstream num_obs_stream;
  std::ofstream num_vars_stream;
  std::ofstream lambda_poisson_stream;
  std::ofstream tree_table_stream;
  std::ofstream tree_matrix_stream;
  std::ofstream alpha_stream;
  std::ofstream beta_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double spike_tree  = RcppDeepState_double();
  spike_tree_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_tree_prior/inputs/spike_tree");
  spike_tree_stream << spike_tree;
  std::cout << "spike_tree values: "<< spike_tree << std::endl;
  spike_tree_stream.close();
  int s_t_hyperprior  = RcppDeepState_int();
  s_t_hyperprior_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_tree_prior/inputs/s_t_hyperprior");
  s_t_hyperprior_stream << s_t_hyperprior;
  std::cout << "s_t_hyperprior values: "<< s_t_hyperprior << std::endl;
  s_t_hyperprior_stream.close();
  double p_s_t  = RcppDeepState_double();
  p_s_t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_tree_prior/inputs/p_s_t");
  p_s_t_stream << p_s_t;
  std::cout << "p_s_t values: "<< p_s_t << std::endl;
  p_s_t_stream.close();
  double a_s_t  = RcppDeepState_double();
  a_s_t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_tree_prior/inputs/a_s_t");
  a_s_t_stream << a_s_t;
  std::cout << "a_s_t values: "<< a_s_t << std::endl;
  a_s_t_stream.close();
  double b_s_t  = RcppDeepState_double();
  b_s_t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_tree_prior/inputs/b_s_t");
  b_s_t_stream << b_s_t;
  std::cout << "b_s_t values: "<< b_s_t << std::endl;
  b_s_t_stream.close();
  double num_obs  = RcppDeepState_double();
  num_obs_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_tree_prior/inputs/num_obs");
  num_obs_stream << num_obs;
  std::cout << "num_obs values: "<< num_obs << std::endl;
  num_obs_stream.close();
  double num_vars  = RcppDeepState_double();
  num_vars_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_tree_prior/inputs/num_vars");
  num_vars_stream << num_vars;
  std::cout << "num_vars values: "<< num_vars << std::endl;
  num_vars_stream.close();
  double lambda_poisson  = RcppDeepState_double();
  lambda_poisson_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_tree_prior/inputs/lambda_poisson");
  lambda_poisson_stream << lambda_poisson;
  std::cout << "lambda_poisson values: "<< lambda_poisson << std::endl;
  lambda_poisson_stream.close();
  NumericMatrix tree_table  = RcppDeepState_NumericMatrix();
  tree_table_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_tree_prior/inputs/tree_table");
  tree_table_stream << tree_table;
  std::cout << "tree_table values: "<< tree_table << std::endl;
  tree_table_stream.close();
  NumericMatrix tree_matrix  = RcppDeepState_NumericMatrix();
  tree_matrix_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_tree_prior/inputs/tree_matrix");
  tree_matrix_stream << tree_matrix;
  std::cout << "tree_matrix values: "<< tree_matrix << std::endl;
  tree_matrix_stream.close();
  double alpha  = RcppDeepState_double();
  alpha_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_tree_prior/inputs/alpha");
  alpha_stream << alpha;
  std::cout << "alpha values: "<< alpha << std::endl;
  alpha_stream.close();
  double beta  = RcppDeepState_double();
  beta_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/get_tree_prior/inputs/beta");
  beta_stream << beta;
  std::cout << "beta values: "<< beta << std::endl;
  beta_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    get_tree_prior(spike_tree,s_t_hyperprior,p_s_t,a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table,tree_matrix,alpha,beta);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
