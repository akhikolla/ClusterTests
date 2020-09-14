#ifndef UPDATE_H
#define UPDATE_H

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <vector>
# include <cstdlib>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
using namespace std;

using namespace Rcpp;
const double log_e = log(exp(1));
const double Inf = std::numeric_limits<double>::infinity();



// [[Rcpp::depends("RcppArmadillo")]]
//****************************************************************************80
int i4_uniform_ab ( int a, int b, int &seed );

double normal_01_cdf ( double x );
double normal_01_cdf_inv ( double cdf );
void normal_01_cdf_values ( int &n_data, double &x, double &fx );
double normal_01_mean ( );
double normal_01_moment ( int order );
double normal_01_pdf ( double x );
double normal_01_sample ( int &seed );
double normal_01_variance ( );

double normal_ms_cdf ( double x, double mu, double sigma );
double normal_ms_cdf_inv ( double cdf, double mu, double sigma );
double normal_ms_mean ( double mu, double sigma );
double normal_ms_moment ( int order, double mu, double sigma );
double normal_ms_moment_central ( int order, double mu, double sigma );
double normal_ms_moment_central_values ( int order, double mu, double sigma );
double normal_ms_moment_values ( int order, double mu, double sigma );
double normal_ms_pdf ( double x, double mu, double sigma );
double normal_ms_sample ( double mu, double sigma, int &seed );
double normal_ms_variance ( double mu, double sigma );

double r8_abs ( double x );
double r8_choose ( int n, int k );
double r8_factorial2 ( int n );
void r8_factorial2_values ( int &n_data, int &n, double &f );
double r8_huge ( );
double r8_log_2 ( double x );
double r8_mop ( int i );
double r8_uniform_01 ( int &seed );

double r8poly_value_horner ( int n, double a[], double x );

double *r8vec_linspace_new ( int n, double a_first, double a_last );
double r8vec_max ( int n, double x[] );
double r8vec_mean ( int n, double x[] );
double r8vec_min ( int n, double x[] );
double r8vec_variance ( int n, double x[] );


double truncated_normal_ab_cdf ( double x, double mu, double sigma, double a,
                                 double b );
void truncated_normal_ab_cdf_values ( int &n_data, double &mu, double &sigma,
                                      double &a, double &b, double &x, double &fx );
double truncated_normal_ab_cdf_inv ( double cdf, double mu, double sigma, double a,
                                     double b );
double truncated_normal_ab_mean ( double mu, double sigma, double a, double b );
double truncated_normal_ab_moment ( int order, double mu, double sigma, double a,
                                    double b );
double truncated_normal_ab_pdf ( double x, double mu, double sigma, double a,
                                 double b );
void truncated_normal_ab_pdf_values ( int &n_data, double &mu, double &sigma,
                                      double &a, double &b, double &x, double &fx );
double truncated_normal_ab_sample ( double mu, double sigma, double a, double b,
                                    int &seed );
double truncated_normal_ab_variance ( double mu, double sigma, double a, double b );

double truncated_normal_a_cdf ( double x, double mu, double sigma, double a );
void truncated_normal_a_cdf_values ( int &n_data, double &mu, double &sigma,
                                     double &a, double &x, double &fx );
double truncated_normal_a_cdf_inv ( double cdf, double mu, double sigma, double a );
double truncated_normal_a_mean ( double mu, double sigma, double a );
double truncated_normal_a_moment ( int order, double mu, double sigma, double a );
double truncated_normal_a_pdf ( double x, double mu, double sigma, double a );
void truncated_normal_a_pdf_values ( int &n_data, double &mu, double &sigma,
                                     double &a, double &x, double &fx );
double truncated_normal_a_sample ( double mu, double sigma, double a, int &seed );
double truncated_normal_a_variance ( double mu, double sigma, double a );

double truncated_normal_b_cdf ( double x, double mu, double sigma, double b );
void truncated_normal_b_cdf_values ( int &n_data, double &mu, double &sigma,
                                     double &b, double &x, double &fx );
double truncated_normal_b_cdf_inv ( double cdf, double mu, double sigma, double b );
double truncated_normal_b_mean ( double mu, double sigma, double b );
double truncated_normal_b_moment ( int order, double mu, double sigma, double b );
double truncated_normal_b_pdf ( double x, double mu, double sigma, double b );
void truncated_normal_b_pdf_values ( int &n_data, double &mu, double &sigma,
                                     double &b, double &x, double &fx );
double truncated_normal_b_sample ( double mu, double sigma, double b, int &seed );
double truncated_normal_b_variance ( double mu, double sigma, double b );

//****************************************************************************80
arma::rowvec nonnan_ind(arma::rowvec x);

arma::mat rowsome(arma::mat x,arma::rowvec ind);

arma::rowvec ton(int n);

arma::mat colsome(arma::mat x, arma::rowvec ind);

arma::mat subsome(arma::mat x, arma::rowvec ind);

double vectornorm(arma::rowvec x);


double dmvnrm_arma0(arma::rowvec x,
                    arma::rowvec mean,
                    double sigma_square,
                    bool logd);

double dmvnrm_arma(arma::rowvec x,
                   arma::rowvec mean,
                   double sigma_square,
                   bool logd);

int rmunoim(arma::rowvec probs);

int sampleint(int n);


// [[Rcpp::export]]
double kernelC(arma::rowvec x, arma::rowvec y, double theta, double tau);

arma::mat eye(int n);

// [[Rcpp::export]]
arma::mat updateC(arma::mat Z,double theta, double tau);

arma::rowvec getind(arma::rowvec x, double c);

arma::mat matrix(arma::rowvec x, int n);
arma::rowvec getgamma(arma::cube x, int i, int j);
arma::mat minusvalue(arma::mat x, int j, arma::rowvec ind, arma::rowvec y);

double findT(arma::rowvec x);

arma::mat removeii(arma::mat A, int a, int b);

arma::rowvec removei(arma::rowvec A, int a);

arma::vec cDmu(double x, arma::vec y);

arma::rowvec ind(arma::rowvec x, arma::rowvec index);

arma::mat choleskyDecompose(arma::mat A);
arma::mat inverse_lower(arma::mat A);
arma::mat inverse_via_chol(arma::mat A);

// [[Rcpp::export]]
List updateBeta(arma::mat X,
                arma::mat Y,
                arma::mat Z,
                arma::mat delta,
                arma::mat Beta,
                arma::cube Gamma,
                arma::rowvec E,
                arma::mat R,
                double S,
                arma::rowvec Ds,
                arma::mat mustar,
                arma::mat mu,
                double sigma,
                arma::rowvec c,
                arma::mat C,
                arma::mat step,
                arma::mat runif,
                int n, int m, int T0, int p, int q, arma::mat D,
                double theta, double tau);


List updateBetanomiss(arma::mat X,
                      arma::mat Y,
                      arma::mat Z,
                      arma::mat delta,
                      arma::mat Beta,
                      arma::cube Gamma,
                      arma::rowvec E,
                      arma::mat R,
                      double S,
                      arma::rowvec Ds,
                      arma::mat mu,
                      double sigma,
                      arma::mat C,
                      arma::mat step,
                      arma::mat runif,
                      int n, int m, int T0, int p, int q, arma::mat D,
                      double theta, double tau);

int Nuni(arma::rowvec x);

arma::mat getGammai(arma::cube Gamma, int j, int ds, int q);

arma::mat putvalue(arma::mat x, int k, arma::rowvec indi, arma::rowvec value);


arma::rowvec computemu(arma::mat mu, int k, arma::rowvec indi, arma::mat Zi, arma::rowvec gamma, arma::rowvec gammanew);

arma::rowvec computemustar(arma::mat mustar, int k, arma::rowvec indT, arma::rowvec c,
                           arma::mat D, arma::rowvec muk, arma::rowvec muknew, int T0);

arma::cube putGamma(arma::cube Gamma, arma::mat Gammai, int i, int q);
arma::rowvec ind1(arma::rowvec x, arma::rowvec indi, arma::rowvec index);
arma::mat removeiii(arma::mat x, int a);


// [[Rcpp::export]]
List updateGamma(arma::mat X,arma::mat Y,arma::mat Z,
                 arma::mat delta,arma::mat Beta,arma::cube Gamma,arma::rowvec E,arma::mat R,
                 double S,arma::rowvec Ds,arma::mat mu,arma::mat mustar,
                 double sigma,arma::rowvec c,
                 arma::cube step,arma::cube runif,
                 int n, int m, int T0, int p, int q, arma::mat D,
                 double theta, double tau);

List updateGammanomiss(arma::mat X,arma::mat Y,arma::mat Z,
                       arma::mat delta,arma::mat Beta,arma::cube Gamma,
                       arma::rowvec E,arma::mat R,
                       double S,arma::rowvec Ds,arma::mat mu, double sigma,
                       arma::cube step,arma::cube runif,
                       int n, int m, int T0, int p, int q, arma::mat D,
                       double theta, double tau);

// [[Rcpp::export]]
List updateE(arma::mat Beta,arma::cube Gamma,arma::rowvec w,
             arma::mat X, arma::mat Y,arma::mat Z, arma::mat delta,
             arma::rowvec E,arma::mat R,double S,arma::rowvec Ds,
             arma::mat mu,arma::mat mustar,double sigma,arma::rowvec c,
             int n, int m, int T0, int p, int q, arma::mat D);


List updateEnomiss(arma::mat Beta, arma::cube Gamma, arma::rowvec w,
                   arma::mat X, arma::mat Y,arma::mat Z,
                   arma::mat delta,arma::rowvec E,arma::mat R,
                   double S,arma::rowvec Ds,arma::mat mu,double sigma,
                   int n, int m, int T0, int p, int q, arma::mat D);

// [[Rcpp::export]]
arma::mat updatemu(arma::mat R,arma::mat Z,arma::mat X,
                   arma::cube Gamma,arma::rowvec K,arma::mat Beta,
                   arma::rowvec E,
                   int m, int n,int p, int q);


arma::rowvec putvalue_vec(arma::rowvec x, arma::rowvec indi, arma::colvec value);
arma::mat rbind(arma::mat x, arma::rowvec y, arma::rowvec z);
arma::mat rbind0(arma::mat x, arma::rowvec y);
arma::rowvec connect(arma::rowvec x, double y, double z);
arma::rowvec connect0(arma::rowvec x, double y);

double det_via_chol(arma::mat A);

List Split(arma::rowvec w,int K,arma::mat Gamma,arma::rowvec Beta,
           arma::mat X,arma::mat Y, arma::mat Z,arma::rowvec R,
           arma::mat delta,arma::rowvec mu, arma::mat mu_star,
           double c,double sigma_square, arma::mat C,
           double theta, double tau,
           int m, int n, int q, int T0, double hyper_delta);


List Splitnomiss(arma::rowvec w,int K,
                 arma::mat Gamma,arma::rowvec Beta,
                 arma::mat X,arma::mat Y,arma::mat Z,arma::rowvec R,
                 arma::mat delta,arma::rowvec mu,double sigma_square,arma::mat C,
                 double theta, double tau,
                 int m, int n, int q, int T0, double hyper_delta);




arma::rowvec sampleint2(int K);


List Merge(arma::rowvec w,
           int K,
           arma::mat Gamma,
           arma::rowvec Beta,
           arma::mat X,
           arma::mat Y,
           arma::mat Z,
           arma::rowvec R,
           arma::mat delta,
           arma::rowvec mu,
           arma::mat mu_star,
           double c,
           double sigma_square,
           arma::mat C,
           double theta, double tau,
           int m, int n, int q, int T0, double hyper_delta);

List Mergenomiss(arma::rowvec w,
                 int K,
                 arma::mat Gamma,
                 arma::rowvec Beta,
                 arma::mat X,
                 arma::mat Y,
                 arma::mat Z,
                 arma::rowvec R,
                 arma::mat delta,
                 arma::rowvec mu,
                 double sigma_square,
                 arma::mat C,
                 double theta, double tau,
                 int m, int n, int q, int T0, double hyper_delta);

// [[Rcpp::export]]
List RJi(arma::rowvec w,
         int K,
         arma::mat Gamma,
         arma::rowvec Beta,
         arma::mat X,
         arma::mat Y,
         arma::mat Z,
         arma::rowvec R,
         arma::mat delta,
         arma::rowvec mu,
         arma::mat mu_star,
         double c,
         double sigma_square,
         arma::mat C,
         double theta, double tau,
         int m, int n, int q, int T0, double hyper_delta );


List RJinomiss(arma::rowvec w,
               int K,
               arma::mat Gamma,
               arma::rowvec Beta,
               arma::mat X,
               arma::mat Y,
               arma::mat Z,
               arma::rowvec R,
               arma::mat delta,
               arma::rowvec mu,
               double sigma_square,
               arma::mat C,
               double theta, double tau,
               int m, int n, int q, int T0, double hyper_delta);


arma::rowvec updateRi(arma::rowvec w,
                      arma::mat Gamma,
                      arma::rowvec Beta,
                      arma::mat Y,
                      arma::mat Z,
                      arma::mat delta,
                      arma::rowvec mu,
                      arma::mat mu_star,
                      double c,
                      double sigma_square,
                      int K, arma::mat X,
                      int m, int n, int q);


arma::rowvec updateRinomiss(arma::rowvec w,
                            arma::mat Gamma,
                            arma::rowvec Beta,
                            arma::mat Y,
                            arma::mat Z,
                            arma::mat delta,
                            arma::rowvec mu,
                            double sigma_square,
                            int K, arma::mat X,
                            int m, int n, int q);

arma::rowvec sample_int(int size, arma::rowvec pro);


// [[Rcpp::export]]
arma::mat updateR(arma::mat w,
                  arma::cube Gamma,
                  arma::mat Beta,
                  arma::mat Y,
                  arma::mat Z,
                  arma::mat delta,
                  arma::mat mu,
                  arma::mat mu_star,
                  double c, int S,
                  double sigma_square,
                  arma::rowvec K,
                  arma::rowvec E,
                  arma::mat X,
                  int m, int n, int q, int p, int T0);

arma::mat updateRnomiss(arma::mat w,
                        arma::cube Gamma,
                        arma::mat Beta,
                        arma::mat Y,
                        arma::mat Z,
                        arma::mat delta,
                        arma::mat mu,
                        int S,
                        double sigma_square,
                        arma::rowvec K,
                        arma::rowvec E,
                        arma::mat X,
                        int m, int n, int q, int p, int T0);
// [[Rcpp::export]]
arma::mat updatemustar(arma::mat mu,
                       arma::rowvec c,
                       int n, int T0, arma::mat D);

double rtrun_norm( double mu, double sigma, double a, double b);

// [[Rcpp::export]]
arma::mat updateZstar(arma::mat mu_star, arma::mat delta,
                      int n, int T0);

arma::rowvec rmvnorm_row(arma::rowvec mean, arma::mat var);


// [[Rcpp::export]]
arma::rowvec updatec(arma::mat Zstar, arma::mat mu, arma::mat D,
                     double sigmac, double sigma_square, int n, int T0);


List updateBetanoDPP(arma::mat X,
                     arma::mat Y,
                     arma::mat Z,
                     arma::mat delta,
                     arma::mat Beta,
                     arma::cube Gamma,
                     arma::rowvec E,
                     arma::mat R,
                     double S,
                     arma::rowvec Ds,
                     arma::mat mustar,
                     arma::mat mu,
                     double sigma,
                     arma::rowvec c,
                     arma::mat step,
                     arma::mat runif,
                     int n, int m, int T0, int p, int q, arma::mat D);



List updateGammanoDPP(arma::mat X,
                      arma::mat Y,
                      arma::mat Z,
                      arma::mat delta,
                      arma::mat Beta,
                      arma::cube Gamma,
                      arma::rowvec E,
                      arma::mat R,
                      double S,
                      arma::rowvec Ds,
                      arma::mat mu,
                      arma::mat mustar,
                      double sigma,
                      arma::rowvec c,
                      arma::cube step,
                      arma::cube runif,
                      int n, int m, int T0, int p, int q, arma::mat D);


arma::cube spcopyGamma2(arma::cube Gamma, arma::rowvec Ds, int S, int q);

arma::cube spcopyGamma(arma::cube Gamma, int k, arma::mat Gammai, arma::rowvec Ds, int S, int q);





List updateRnoDPP(arma::mat Beta,
                  arma::cube Gamma,
                  arma::mat X,
                  arma::mat Y,
                  arma::mat Z,
                  arma::mat delta,
                  arma::rowvec E,
                  arma::mat R,
                  int S,
                  arma::rowvec Ds,
                  arma::mat mu,
                  arma::mat mustar,
                  double sigma,
                  double c,
                  arma::rowvec Ms,
                  int n, int m, int T0, int p, int q, arma::mat D, int Dmax);


List Split_empty(arma::rowvec w,
                 int K,
                 arma::mat Gamma,
                 arma::rowvec R,
                 arma::mat C,
                 double theta, double tau,
                 int m, int q, double hyper_delta);


List Merge_empty(arma::rowvec w,
                 int K,
                 arma::mat Gamma,
                 arma::rowvec R,
                 arma::mat C,
                 double theta, double tau,
                 int m,  int q, double hyper_delta);

// [[Rcpp::export]]
List RJi_empty(arma::rowvec w,
               int K,
               arma::mat Gamma,
               arma::rowvec R,
               arma::mat C,
               double theta, double tau,
               int m,  int q,  double hyper_delta);

#endif
