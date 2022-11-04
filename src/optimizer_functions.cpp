#include "causalOT_types.h"
#include "utils.h"

// COT with entropy penalty, not debiased
// [[Rcpp::export]]
double cotEntropy_obj_(const SEXP & vars_, const SEXP & source_,
                       const SEXP & target_,
                       const SEXP& cost_, 
                       const SEXP& b_,
                       double delta,
                       double lambda) 
{
  
  // map R types to Eigen types
  const vecMap vars(Rcpp::as<vecMap>(vars_));
  const matMap source(Rcpp::as< matMap >(source_));
  const vecMap target(Rcpp::as< vecMap >(target_));
  const matMap cost(Rcpp::as<matMap>(cost_));
  const vecMap b(Rcpp::as<vecMap>(b_));
  
  // Dim of balance functions
  int K = source.cols();
  int N = cost.rows();
  int M = cost.cols();
  
  // Pull out variables and get linear predictor
  vector g = vars.head(M);
  double objective = g.dot(b);
  matrix eta = (-cost).rowwise() + g.transpose();
  
  if ( K > 0 ) {
    vector beta_const = vars.tail(2 * K);
    vector beta = beta_const.head(K) - beta_const.tail(K);
    eta.noalias() = eta.colwise() - source * beta;
    objective -= delta * beta_const.sum() + target.dot(beta);
  }
  
  eta.array() /= lambda;
  objective -= lambda * logSumExp(eta);
  
  return(-objective); // negative b/c lbfgs minimizes by default and we need to maximize
}

// [[Rcpp::export]]
Rcpp::NumericVector cotEntropy_grad_(const SEXP & vars_, const SEXP & source_,
                                     const SEXP & target_,
                                     const SEXP& cost_, 
                                     const SEXP& b_,
                                     double delta,
                                     double lambda) 
{
  // map R types to Eigen types
  const vecMap vars(Rcpp::as<vecMap>(vars_));
  const matMap source(Rcpp::as< matMap >(source_));
  const vecMap target(Rcpp::as< vecMap >(target_));
  const matMap cost(Rcpp::as<matMap>(cost_));
  const vecMap b(Rcpp::as<vecMap>(b_));
  
  // Dim of balance functions
  int K = source.cols();
  int N = cost.rows();
  int M = cost.cols();
  
  // gradient vector
  vector grad = vector::Zero(M + 2 * K);
  
  // Pull out variables and start forming dual form of transport matrix
  vector g = vars.head(M); // dual for margin constraint
  matrix eta = (-cost).rowwise() + g.transpose();
  
  if ( K > 0 ) { // if balance functions passed
    vector beta_const = vars.tail(2 * K); // get all of the 0 constrained coefficients
    vector beta = beta_const.head(K) - beta_const.tail(K); // get unconstrained coef
    eta.colwise() -= source * beta; // add linear predictor to cost
    // grad.tail(2*K).array() -= delta;
    for (int i = 0; i < (2 * K); ++i) {
      if (beta_const(i) > 0.) grad(M + i) -= delta; // derivative of penalty term
      // grad(M + i) -= delta; // derivative of penalty term
    }
  }
  
  eta.array() /= lambda;
  double e_sum = logSumExp(eta);
  eta.array() -= e_sum;
  matrix pi = eta.array().exp();
  
  grad.head(M) = b; // derivative of g, margin dual, first term
  grad.head(M) -= pi.colwise().sum(); // derivative of g, pi term
  
  if (K > 0) {
    vector pi_a = pi.rowwise().sum();
    vector cross = (source.transpose() * pi_a);
    grad.tail(K) += target;
    grad.block(M,0, K,1) -= target;
    
    grad.tail(K) -= cross; // negative vars derivative wrt transport matrix
    grad.block(M,0, K,1) += cross; // positive vars derivate wrt transport mat
  }
  // convert to R vector
  Rcpp::NumericVector output = Rcpp::wrap(-grad); // negative b/c lbfgs minimizes by default and we need to maximize
  
  return(output); 
}



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
