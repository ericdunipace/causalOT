#include "causalOT_types.h"
#include "cost_lp.h"

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

matrix invsqrt(const refMatConst & A) {

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  if ((es.eigenvalues().array() <= 0).any()) Rcpp::warning("Covariance matrix is not positive definite. Using approximate inverse square-root based on SVD");
  
  Eigen::VectorXd D_half_inv = 1.0 / es.eigenvalues().array().abs().sqrt();
  // Rcpp::Rcout << D_half_inv << std::endl;
  return(es.eigenvectors() * D_half_inv.asDiagonal() * es.eigenvectors().transpose());
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
                                const double p,
                                const std::string estimand) {
  int N = A_.cols();
  int M = B_.cols();
  
  const matMap A(Rcpp::as<matMap >(A_));
  const matMap B(Rcpp::as<matMap >(B_));
  
  matrix cov(A.rows(), A.rows());
  
  // if (estimand == "ATE") {
  //   matrix C(A.rows(), A.cols() + B.cols());
  //   C << A, B;
  //   cov  = covariance(C);
  // } else if (estimand == "ATT" ) {
  //   cov = covariance(B);
  // } else if (estimand == "ATC" ) {
  //   cov = covariance(A);
  // } else {
  //   double frac_A = double(N)/double(N + M);
  //   double frac_B = double(M)/double(N + M);
  //   
  //   const matrix covA = covariance(A);
  //   const matrix covB = covariance(B);
  //   
  //   cov = frac_A * covA + frac_B * covB;
  // }
  
  if (estimand == "ATE") {
    cov = covariance(B);
  } else {
    matrix C(A.rows(), A.cols() + B.cols());
    C << A, B;
    cov  = covariance(C);
  }
  
  

  const matrix L_inv = invsqrt(cov); //cov inv sqrt
  
  matrix cost_matrix(N,M);
  
  // more memory efficient:
  // if(p == 2.0) {
  //   cost_mahal_L2(A, B, L_inv, cost_matrix);
  // } else if (p == 1.0){
  //   cost_mahal_L1(A, B, L_inv, cost_matrix);
  // } else {
  //   cost_mahal_Lp(A, B, L_inv, cost_matrix, p);
  // }
  matrix AL = L_inv * A;
  matrix BL = L_inv * B;
  // Rcpp::Rcout << std::endl << AL(0,0) << std::endl;
  if(p == 2.0) {
    cost_calculation_L2(AL, BL, cost_matrix);
  } else if (p == 1.0){
    cost_calculation_L1(AL, BL, cost_matrix);
  } else {
    cost_calculation_Lp(AL, BL, cost_matrix, p);
  }
  
  return Rcpp::wrap(cost_matrix);
}
