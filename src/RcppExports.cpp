// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cumsumBorder
NumericVector cumsumBorder(NumericVector x, double border);
RcppExport SEXP RnaSeqSampleSize_cumsumBorder(SEXP xSEXP, SEXP borderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type border(borderSEXP);
    rcpp_result_gen = Rcpp::wrap(cumsumBorder(x, border));
    return rcpp_result_gen;
END_RCPP
}
