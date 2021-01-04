#include "cost_lp.h"

void cost_calculation_Lp(const refMatConst & A, const refMatConst & B, matrix & cost_matrix, double p) {
  double p_inv = 1.0/p;
  
  for (int j = 0; j < B.cols(); j++) { 
    vector bvec = B.col(j);
    for (int i = 0; i < A.cols(); i++) {
      double cost_p = (A.col(i)-bvec).array().abs().pow(p).sum();
      cost_matrix(i,j) = std::pow(cost_p, p_inv);
    }
  }
}

void cost_calculation_L2(const refMatConst & A, const refMatConst & B, matrix & cost_matrix) {
  for (int j = 0; j < B.cols(); j++) { 
    vector bvec = B.col(j);
    for (int i = 0; i < A.cols(); i++) {
      cost_matrix(i,j) = (A.col(i)-bvec).norm();
    }
  }
}

void cost_calculation_L2sq(const refMatConst & A, const refMatConst & B, matrix & cost_matrix) {
  for (int j = 0; j < B.cols(); j++) { 
    vector bvec = B.col(j);
    for (int i = 0; i < A.cols(); i++) {
      cost_matrix(i,j) = (A.col(i)-bvec).squaredNorm();
    }
  }
}

void cost_calculation_L1(const refMatConst & A, const refMatConst & B, matrix & cost_matrix) {
  for (int j = 0; j < B.cols(); j++) { 
    vector bvec = B.col(j);
    for (int i = 0; i < A.cols(); i++) {
      cost_matrix(i,j) = (A.col(i)-bvec).cwiseAbs().sum();
    }
  }
}

//[[Rcpp::export]]
Rcpp::NumericMatrix cost_calculation_(const Rcpp::NumericMatrix & A_, const Rcpp::NumericMatrix & B_, const double p) {
  int N = A_.cols();
  int M = B_.cols();
  
  const matMap A(Rcpp::as<matMap >(A_));
  const matMap B(Rcpp::as<matMap >(B_));
  
  matrix cost_matrix(N,M);
  
  if(p == 2.0) {
    cost_calculation_L2(A, B, cost_matrix);
  } else if (p == 1.0){
    cost_calculation_L1(A, B, cost_matrix);
  } else {
    cost_calculation_Lp(A, B, cost_matrix, p);
  }
  
  return Rcpp::wrap(cost_matrix);
}
