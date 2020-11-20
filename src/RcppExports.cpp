// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "causalOT_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// cost_calculation_
Rcpp::NumericMatrix cost_calculation_(const Rcpp::NumericMatrix& A_, const Rcpp::NumericMatrix& B_, const double p);
RcppExport SEXP _causalOT_cost_calculation_(SEXP A_SEXP, SEXP B_SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type A_(A_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type B_(B_SEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(cost_calculation_(A_, B_, p));
    return rcpp_result_gen;
END_RCPP
}
// cost_mahal_
Rcpp::NumericMatrix cost_mahal_(const Rcpp::NumericMatrix& A_, const Rcpp::NumericMatrix& B_, const double p);
RcppExport SEXP _causalOT_cost_mahal_(SEXP A_SEXP, SEXP B_SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type A_(A_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type B_(B_SEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(cost_mahal_(A_, B_, p));
    return rcpp_result_gen;
END_RCPP
}
// kernel_calc_dose_
Rcpp::List kernel_calc_dose_(const Rcpp::NumericMatrix& X_, const Rcpp::NumericMatrix& z_, const double p, const Rcpp::NumericVector& theta_, const Rcpp::NumericVector& gamma_, const std::string& kernel_, const bool calc_covariance);
RcppExport SEXP _causalOT_kernel_calc_dose_(SEXP X_SEXP, SEXP z_SEXP, SEXP pSEXP, SEXP theta_SEXP, SEXP gamma_SEXP, SEXP kernel_SEXP, SEXP calc_covarianceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type z_(z_SEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma_(gamma_SEXP);
    Rcpp::traits::input_parameter< const std::string& >::type kernel_(kernel_SEXP);
    Rcpp::traits::input_parameter< const bool >::type calc_covariance(calc_covarianceSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_calc_dose_(X_, z_, p, theta_, gamma_, kernel_, calc_covariance));
    return rcpp_result_gen;
END_RCPP
}
// similarity_calc_dose_
Rcpp::List similarity_calc_dose_(const Rcpp::NumericMatrix& X_, const Rcpp::NumericMatrix& z_, const bool calc_covariance);
RcppExport SEXP _causalOT_similarity_calc_dose_(SEXP X_SEXP, SEXP z_SEXP, SEXP calc_covarianceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type z_(z_SEXP);
    Rcpp::traits::input_parameter< const bool >::type calc_covariance(calc_covarianceSEXP);
    rcpp_result_gen = Rcpp::wrap(similarity_calc_dose_(X_, z_, calc_covariance));
    return rcpp_result_gen;
END_RCPP
}
// kernel_calc_
Rcpp::List kernel_calc_(const Rcpp::NumericMatrix& X_, const Rcpp::IntegerVector& z, const double p, const Rcpp::NumericVector& theta_, const Rcpp::NumericVector& gamma_, const Rcpp::NumericVector& sigma_2_, const std::string& kernel_, const bool calc_covariance, const std::string& estimand);
RcppExport SEXP _causalOT_kernel_calc_(SEXP X_SEXP, SEXP zSEXP, SEXP pSEXP, SEXP theta_SEXP, SEXP gamma_SEXP, SEXP sigma_2_SEXP, SEXP kernel_SEXP, SEXP calc_covarianceSEXP, SEXP estimandSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma_(gamma_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type sigma_2_(sigma_2_SEXP);
    Rcpp::traits::input_parameter< const std::string& >::type kernel_(kernel_SEXP);
    Rcpp::traits::input_parameter< const bool >::type calc_covariance(calc_covarianceSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type estimand(estimandSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_calc_(X_, z, p, theta_, gamma_, sigma_2_, kernel_, calc_covariance, estimand));
    return rcpp_result_gen;
END_RCPP
}
// kernel_calc_pred_
Rcpp::List kernel_calc_pred_(const Rcpp::NumericMatrix& X_, const Rcpp::NumericMatrix& X_test_, const Rcpp::IntegerVector& z, const double p, const Rcpp::NumericVector& theta_, const Rcpp::NumericVector& gamma_, const Rcpp::NumericVector& sigma_2_, const std::string& kernel_, const bool calc_covariance, const std::string& estimand);
RcppExport SEXP _causalOT_kernel_calc_pred_(SEXP X_SEXP, SEXP X_test_SEXP, SEXP zSEXP, SEXP pSEXP, SEXP theta_SEXP, SEXP gamma_SEXP, SEXP sigma_2_SEXP, SEXP kernel_SEXP, SEXP calc_covarianceSEXP, SEXP estimandSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X_test_(X_test_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma_(gamma_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type sigma_2_(sigma_2_SEXP);
    Rcpp::traits::input_parameter< const std::string& >::type kernel_(kernel_SEXP);
    Rcpp::traits::input_parameter< const bool >::type calc_covariance(calc_covarianceSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type estimand(estimandSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_calc_pred_(X_, X_test_, z, p, theta_, gamma_, sigma_2_, kernel_, calc_covariance, estimand));
    return rcpp_result_gen;
END_RCPP
}
// kernel_update_
Rcpp::NumericMatrix kernel_update_(const Rcpp::NumericMatrix& sim_, const Rcpp::IntegerVector& z_, const double p, const Rcpp::NumericVector& theta_, const Rcpp::NumericVector& gamma_, const Rcpp::NumericVector& sigma_2_, const std::string& kernel_);
RcppExport SEXP _causalOT_kernel_update_(SEXP sim_SEXP, SEXP z_SEXP, SEXP pSEXP, SEXP theta_SEXP, SEXP gamma_SEXP, SEXP sigma_2_SEXP, SEXP kernel_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type sim_(sim_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type z_(z_SEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma_(gamma_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type sigma_2_(sigma_2_SEXP);
    Rcpp::traits::input_parameter< const std::string& >::type kernel_(kernel_SEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_update_(sim_, z_, p, theta_, gamma_, sigma_2_, kernel_));
    return rcpp_result_gen;
END_RCPP
}
// similarity_calc_
matrix similarity_calc_(const Rcpp::NumericMatrix& X_, const Rcpp::IntegerVector& z, const bool calc_covariance, const std::string& estimand);
RcppExport SEXP _causalOT_similarity_calc_(SEXP X_SEXP, SEXP zSEXP, SEXP calc_covarianceSEXP, SEXP estimandSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const bool >::type calc_covariance(calc_covarianceSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type estimand(estimandSEXP);
    rcpp_result_gen = Rcpp::wrap(similarity_calc_(X_, z, calc_covariance, estimand));
    return rcpp_result_gen;
END_RCPP
}
// marginal_lik_gp_
double marginal_lik_gp_(const Rcpp::NumericVector& y_, const Rcpp::NumericMatrix& K_);
RcppExport SEXP _causalOT_marginal_lik_gp_(SEXP y_SEXP, SEXP K_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y_(y_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type K_(K_SEXP);
    rcpp_result_gen = Rcpp::wrap(marginal_lik_gp_(y_, K_));
    return rcpp_result_gen;
END_RCPP
}
// kernel_calc_ot_
Rcpp::List kernel_calc_ot_(const Rcpp::NumericMatrix& X_, const Rcpp::IntegerVector& z, const double p, const Rcpp::NumericVector& theta_, const Rcpp::NumericVector& gamma_, const std::string& kernel_, const bool calc_covariance, const std::string& estimand);
RcppExport SEXP _causalOT_kernel_calc_ot_(SEXP X_SEXP, SEXP zSEXP, SEXP pSEXP, SEXP theta_SEXP, SEXP gamma_SEXP, SEXP kernel_SEXP, SEXP calc_covarianceSEXP, SEXP estimandSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma_(gamma_SEXP);
    Rcpp::traits::input_parameter< const std::string& >::type kernel_(kernel_SEXP);
    Rcpp::traits::input_parameter< const bool >::type calc_covariance(calc_covarianceSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type estimand(estimandSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_calc_ot_(X_, z, p, theta_, gamma_, kernel_, calc_covariance, estimand));
    return rcpp_result_gen;
END_RCPP
}
// entry
void entry(SEXP& xx, Rcpp::NumericMatrix& y, int colX_);
RcppExport SEXP _causalOT_entry(SEXP xxSEXP, SEXP ySEXP, SEXP colX_SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type colX_(colX_SEXP);
    entry(xx, y, colX_);
    return R_NilValue;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4barycenter_projection__mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4barycenter_projection_mahalanobis__mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4gp_hyper_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_causalOT_cost_calculation_", (DL_FUNC) &_causalOT_cost_calculation_, 3},
    {"_causalOT_cost_mahal_", (DL_FUNC) &_causalOT_cost_mahal_, 3},
    {"_causalOT_kernel_calc_dose_", (DL_FUNC) &_causalOT_kernel_calc_dose_, 7},
    {"_causalOT_similarity_calc_dose_", (DL_FUNC) &_causalOT_similarity_calc_dose_, 3},
    {"_causalOT_kernel_calc_", (DL_FUNC) &_causalOT_kernel_calc_, 9},
    {"_causalOT_kernel_calc_pred_", (DL_FUNC) &_causalOT_kernel_calc_pred_, 10},
    {"_causalOT_kernel_update_", (DL_FUNC) &_causalOT_kernel_update_, 7},
    {"_causalOT_similarity_calc_", (DL_FUNC) &_causalOT_similarity_calc_, 4},
    {"_causalOT_marginal_lik_gp_", (DL_FUNC) &_causalOT_marginal_lik_gp_, 2},
    {"_causalOT_kernel_calc_ot_", (DL_FUNC) &_causalOT_kernel_calc_ot_, 8},
    {"_causalOT_entry", (DL_FUNC) &_causalOT_entry, 3},
    {"_rcpp_module_boot_stan_fit4barycenter_projection__mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4barycenter_projection__mod, 0},
    {"_rcpp_module_boot_stan_fit4barycenter_projection_mahalanobis__mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4barycenter_projection_mahalanobis__mod, 0},
    {"_rcpp_module_boot_stan_fit4gp_hyper_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gp_hyper_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_causalOT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
