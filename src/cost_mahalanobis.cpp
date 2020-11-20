#include "causalOT_types.h"

matrix covariance(const refMatConst & samples) {
  int S = samples.cols();
  int d = samples.rows();
  // matrix c_samples(d, S);
  vector mean = samples.rowwise().mean();
  if(d != mean.rows()) Rcpp::stop("Dimension of mean vector not match dimension of samples vector!");
  
  matrix c_samples = samples.colwise() - mean;
  
  // for(int i = 0 ; i < S; i++){
  //   c_samples.col(i) = samples.col(i) - mean;
  // }
  return matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(c_samples, 1.0/double(S-1));
}

void cost_mahal_Lp(const refMatConst & A, const refMatConst & B, 
                        const matrix & L, matrix & cost_matrix, double p) {
  double p_inv = 1.0/p;
  
  for (int j = 0; j < B.cols(); j++) { 
    vector bvec = B.col(j);
    for (int i = 0; i < A.cols(); i++) {
      double cost_p = (L.triangularView<Eigen::Lower>()*(A.col(i)-bvec)).array().abs().pow(p).sum();
      cost_matrix(i,j) = std::pow(cost_p, p_inv);
    }
  }
}

void cost_mahal_L2(const refMatConst & A, const refMatConst & B, 
                         const matrix & L, matrix & cost_matrix) {
  for (int j = 0; j < B.cols(); j++) { 
    vector bvec = B.col(j);
    for (int i = 0; i < A.cols(); i++) {
      cost_matrix(i,j) = (L.triangularView<Eigen::Lower>()*(A.col(i)-bvec)).norm();
    }
  }
}

void cost_calculation_L2sq(const refMatConst & A, const refMatConst & B, 
                           const matrix & L, matrix & cost_matrix) {
  for (int j = 0; j < B.cols(); j++) { 
    vector bvec = B.col(j);
    for (int i = 0; i < A.cols(); i++) {
      cost_matrix(i,j) = (L.triangularView<Eigen::Lower>().solve(A.col(i)-bvec)).squaredNorm();
    }
  }
}

void cost_mahal_L1(const refMatConst & A, const refMatConst & B, 
                   const matrix & L, matrix & cost_matrix) {
  for (int j = 0; j < B.cols(); j++) { 
    vector bvec = B.col(j);
    for (int i = 0; i < A.cols(); i++) {
      cost_matrix(i,j) = (L.triangularView<Eigen::Lower>() * (A.col(i)-bvec)).cwiseAbs().sum();
    }
  }
}

//[[Rcpp::export]]
Rcpp::NumericMatrix cost_mahal_(const Rcpp::NumericMatrix & A_, 
                                const Rcpp::NumericMatrix & B_, 
                                const double p) {
  int N = A_.cols();
  int M = B_.cols();
  
  const matMap A(Rcpp::as<matMap >(A_));
  const matMap B(Rcpp::as<matMap >(B_));
  
  const matrix covA = covariance(A);
  const matrix covB = covariance(B);
  const matrix cov = 0.5 * covA + 0.5 * covB;
  const matrix L = cov.selfadjointView<Eigen::Lower>().llt().matrixL();
  const matrix L_inv = L.inverse();
  // Rcpp::Rcout << covA(0,0)<<std::endl;
  // Rcpp::Rcout << cov(0,0);
  
  matrix cost_matrix(N,M);
  
  if(p == 2.0) {
    cost_mahal_L2(A, B, L_inv, cost_matrix);
  } else if (p == 1.0){
    cost_mahal_L1(A, B, L_inv, cost_matrix);
  } else {
    cost_mahal_Lp(A, B, L_inv, cost_matrix, p);
  }
  
  return Rcpp::wrap(cost_matrix);
}
