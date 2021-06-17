#include "causalOT_types.h"

// [[Rcpp::export]]
double cotDualL2_2_obj_(const SEXP & vars_, const SEXP & QQ,
                        const SEXP& cost_, double pf_,
                        double lambda) {
  // map R types to Eigen types
  const vecMap vars(Rcpp::as<vecMap>(vars_));
  const Eigen::MappedSparseMatrix<double> Q(Rcpp::as< Eigen::MappedSparseMatrix<double> >(QQ));
  const vecMap cost(Rcpp::as<vecMap>(cost_));
  double normalization = 1. /double(Q.rows());

  // calculate main objective and check for positivity
  vector eta = Q * vars - cost;
  vector pos = (eta.array() > 0.0).cast<double>();
  SpVec possp = pos.sparseView();
  double diff_sq = (eta.cwiseProduct(possp) ).squaredNorm();
  
  // add penalty factor to gradient
  double obj = pf_ * normalization;
  
  // add main objective
  obj += -0.5 / lambda * diff_sq * normalization;
  
  return(-obj); // negative b/c lbfgs minimizes by default and we need to maximize
}

// [[Rcpp::export]]
Rcpp::NumericVector cotDualL2_2_grad_(const SEXP & vars_, 
                                      const SEXP & QQ,
                        const SEXP& cost_, const SEXP & pf_,
                        double lambda) {
  // map R types to Eigen types
  const vecMap vars(Rcpp::as<vecMap>(vars_));
  const Eigen::MappedSparseMatrix<double> Q(Rcpp::as< Eigen::MappedSparseMatrix<double> >(QQ));
  const vecMap cost(Rcpp::as<vecMap>(cost_));
  const vecMap pf(Rcpp::as<vecMap>(pf_));
  
  // calculate main objective and check for positivity
  double normalization = 1./double(Q.rows());
  vector eta = Q * vars - cost;
  vector pos = (eta.array() > 0).cast<double>();
  SpVec possp = pos.sparseView();
  SpVec diff = eta.cwiseProduct(possp);

  // add penalty factors to the gradient
  vector grad = pf * normalization;

  // add main objective derivative to gradient
  grad += (Q.transpose() * -diff) / lambda * normalization ;

  // convert to R vector
  Rcpp::NumericVector output = Rcpp::wrap(-grad); 
  // negative b/c lbfgs minimizes by default and we need to maximize

  return(output);
}
