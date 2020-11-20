// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bootstrapCpp
arma::mat bootstrapCpp(arma::mat X, int B, double dt);
RcppExport SEXP _cods_bootstrapCpp(SEXP XSEXP, SEXP BSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrapCpp(X, B, dt));
    return rcpp_result_gen;
END_RCPP
}
// bootstrapCppOld
arma::mat bootstrapCppOld(arma::mat X, int B, double dt);
RcppExport SEXP _cods_bootstrapCppOld(SEXP XSEXP, SEXP BSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrapCppOld(X, B, dt));
    return rcpp_result_gen;
END_RCPP
}
// exactSimulation
arma::mat exactSimulation(int N, double dt, arma::vec z0, arma::vec freq, arma::vec lvl, double s_phi, double s_gam, const char* model);
RcppExport SEXP _cods_exactSimulation(SEXP NSEXP, SEXP dtSEXP, SEXP z0SEXP, SEXP freqSEXP, SEXP lvlSEXP, SEXP s_phiSEXP, SEXP s_gamSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z0(z0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lvl(lvlSEXP);
    Rcpp::traits::input_parameter< double >::type s_phi(s_phiSEXP);
    Rcpp::traits::input_parameter< double >::type s_gam(s_gamSEXP);
    Rcpp::traits::input_parameter< const char* >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(exactSimulation(N, dt, z0, freq, lvl, s_phi, s_gam, model));
    return rcpp_result_gen;
END_RCPP
}
// factorialCpp
unsigned long long factorialCpp(unsigned long long x, unsigned long long result);
RcppExport SEXP _cods_factorialCpp(SEXP xSEXP, SEXP resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned long long >::type x(xSEXP);
    Rcpp::traits::input_parameter< unsigned long long >::type result(resultSEXP);
    rcpp_result_gen = Rcpp::wrap(factorialCpp(x, result));
    return rcpp_result_gen;
END_RCPP
}
// kronSum
arma::mat kronSum(arma::mat A);
RcppExport SEXP _cods_kronSum(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(kronSum(A));
    return rcpp_result_gen;
END_RCPP
}
// matPow
arma::mat matPow(arma::mat A, int k);
RcppExport SEXP _cods_matPow(SEXP ASEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(matPow(A, k));
    return rcpp_result_gen;
END_RCPP
}
// integralApprox
arma::mat integralApprox(arma::mat P, double t, bool verbose);
RcppExport SEXP _cods_integralApprox(SEXP PSEXP, SEXP tSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(integralApprox(P, t, verbose));
    return rcpp_result_gen;
END_RCPP
}
// Omega
arma::mat Omega(arma::mat P, double t, double s);
RcppExport SEXP _cods_Omega(SEXP PSEXP, SEXP tSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(Omega(P, t, s));
    return rcpp_result_gen;
END_RCPP
}
// johansenCpp
arma::mat johansenCpp(arma::mat X, int r, arma::mat A, arma::mat B, double dt);
RcppExport SEXP _cods_johansenCpp(SEXP XSEXP, SEXP rSEXP, SEXP ASEXP, SEXP BSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(johansenCpp(X, r, A, B, dt));
    return rcpp_result_gen;
END_RCPP
}
// oscillator
arma::mat oscillator(int N, float dt, arma::vec z0, arma::mat alpha, arma::mat beta, arma::vec omega, arma::vec freq, arma::vec lvl, arma::mat S_phi, arma::mat S_gam, const char* model);
RcppExport SEXP _cods_oscillator(SEXP NSEXP, SEXP dtSEXP, SEXP z0SEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP omegaSEXP, SEXP freqSEXP, SEXP lvlSEXP, SEXP S_phiSEXP, SEXP S_gamSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< float >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z0(z0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lvl(lvlSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S_phi(S_phiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S_gam(S_gamSEXP);
    Rcpp::traits::input_parameter< const char* >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(oscillator(N, dt, z0, alpha, beta, omega, freq, lvl, S_phi, S_gam, model));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _cods_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _cods_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _cods_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _cods_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}
// testCpp
arma::mat testCpp(arma::mat A);
RcppExport SEXP _cods_testCpp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(testCpp(A));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cods_bootstrapCpp", (DL_FUNC) &_cods_bootstrapCpp, 3},
    {"_cods_bootstrapCppOld", (DL_FUNC) &_cods_bootstrapCppOld, 3},
    {"_cods_exactSimulation", (DL_FUNC) &_cods_exactSimulation, 8},
    {"_cods_factorialCpp", (DL_FUNC) &_cods_factorialCpp, 2},
    {"_cods_kronSum", (DL_FUNC) &_cods_kronSum, 1},
    {"_cods_matPow", (DL_FUNC) &_cods_matPow, 2},
    {"_cods_integralApprox", (DL_FUNC) &_cods_integralApprox, 3},
    {"_cods_Omega", (DL_FUNC) &_cods_Omega, 3},
    {"_cods_johansenCpp", (DL_FUNC) &_cods_johansenCpp, 5},
    {"_cods_oscillator", (DL_FUNC) &_cods_oscillator, 11},
    {"_cods_rcpparma_hello_world", (DL_FUNC) &_cods_rcpparma_hello_world, 0},
    {"_cods_rcpparma_outerproduct", (DL_FUNC) &_cods_rcpparma_outerproduct, 1},
    {"_cods_rcpparma_innerproduct", (DL_FUNC) &_cods_rcpparma_innerproduct, 1},
    {"_cods_rcpparma_bothproducts", (DL_FUNC) &_cods_rcpparma_bothproducts, 1},
    {"_cods_testCpp", (DL_FUNC) &_cods_testCpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_cods(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
