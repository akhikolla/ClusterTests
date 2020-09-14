
#ifndef _nor_generator_H
#define _nor_generator_H


#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>


arma::mat mvrnorm(arma::vec mu, arma::mat sigma);

arma::mat mvrnorm_Ind(arma::vec mu, double sigma2);


#endif
