#include "causalOT_types.h"
#include "utils.h"


//Entropy Balancing Weights

//[[Rcpp::export]]
double entBW_obj_(const SEXP & vars_, 
                  const SEXP & A_,
                  double delta) {
  
  // map R types to Eigen types
  const vecMap beta_const(Rcpp::as<vecMap>(vars_));
  const matMap A(Rcpp::as<matMap>(A_));
  
  int K = A.cols();
  // int N = A.rows();
  vector beta = beta_const.head(K) - beta_const.tail(K);
  vector lin_pred = A * beta;
  // lin_pred.array() -= std::log(double(N));
  
  return(
    logSumExp(lin_pred) + (beta.array().abs() * delta).sum()
  );
}

//[[Rcpp::export]]
vector entBW_grad_(const SEXP & vars_, 
                   const SEXP & A_,
                   double delta) {
  // map R types to Eigen types
  const vecMap beta_const(Rcpp::as<vecMap>(vars_));
  const matMap A(Rcpp::as<matMap>(A_));
  
  int K = A.cols();
  int N = A.rows();
  vector beta = beta_const.head(K) - beta_const.tail(K);
  
  vector lin_pred = A * beta;
  vector grad_vec = vector::Zero(2 * K);
  vector cross_term = vector::Zero(K);
  double denom = logSumExp(lin_pred);
  
  for (int i = 0; i < 2*K; ++i) {
    if(beta_const(i) > 0) grad_vec(i) = delta;
  }
  cross_term = A.transpose() * (lin_pred.array() - denom).exp().matrix();
  grad_vec.head(K) += cross_term;
  grad_vec.tail(K) -= cross_term;
  
  return( grad_vec );
  
}
