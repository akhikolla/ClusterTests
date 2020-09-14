// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fastMean
NumericVector fastMean(NumericVector xx);
RcppExport SEXP _BayesESS_fastMean(SEXP xxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    rcpp_result_gen = Rcpp::wrap(fastMean(xx));
    return rcpp_result_gen;
END_RCPP
}
// fastProd
arma::mat fastProd(const arma::vec& X, const arma::vec& Y);
RcppExport SEXP _BayesESS_fastProd(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(fastProd(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// fastSum
NumericVector fastSum(NumericVector xx);
RcppExport SEXP _BayesESS_fastSum(SEXP xxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    rcpp_result_gen = Rcpp::wrap(fastSum(xx));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayesESS_fastMean", (DL_FUNC) &_BayesESS_fastMean, 1},
    {"_BayesESS_fastProd", (DL_FUNC) &_BayesESS_fastProd, 2},
    {"_BayesESS_fastSum", (DL_FUNC) &_BayesESS_fastSum, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesESS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
