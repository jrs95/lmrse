// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// robustseEigen
NumericMatrix robustseEigen(NumericMatrix es, NumericMatrix xs, NumericMatrix cs, NumericMatrix xxs, NumericMatrix covs);
RcppExport SEXP _lmrse_robustseEigen(SEXP esSEXP, SEXP xsSEXP, SEXP csSEXP, SEXP xxsSEXP, SEXP covsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type es(esSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cs(csSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xxs(xxsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covs(covsSEXP);
    rcpp_result_gen = Rcpp::wrap(robustseEigen(es, xs, cs, xxs, covs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lmrse_robustseEigen", (DL_FUNC) &_lmrse_robustseEigen, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_lmrse(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
