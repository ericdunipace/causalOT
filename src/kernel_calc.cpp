#include "causalOT_types.h"

void kernel_calculation(const matrix & A, matrix & kernel_matrix, double d,
                           double theta, double gamma) {
  //A must be demeaned and scaled by covariance prior to using
  int n = A.rows();
  kernel_matrix = matrix(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
  kernel_matrix.array() *= theta;
  kernel_matrix.array() += 1.0;
  kernel_matrix.array() = gamma * kernel_matrix.array().pow(d);
}

matrix covariance_kern(const refMatConst & samples) {
  int S = samples.rows();
  int d = samples.cols();
  rowVector mean = samples.colwise().mean();
  
  if(d != mean.cols()) Rcpp::stop("Dimension of mean vector does not match the dimension of samples vector!");
  matrix c_samples = samples.rowwise() - mean;
  return matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(c_samples.transpose(), 1.0/double(S-1));
}

//[[Rcpp::export]]
Rcpp::NumericMatrix kernel_calc_(const Rcpp::NumericMatrix & X_,  //confounders
                                 const Rcpp::NumericMatrix & z_,  //tx, a vector but easier if matrix
                                 const double d,
                                 const Rcpp::NumericVector  & theta_,
                                 const Rcpp::NumericVector & gamma_,
                                 const bool calc_covariance ) {
  
  int N = X_.rows();
  if(N != z_.rows() ) Rcpp::stop("Observations for X and z must be equal!");
  
  const matMap X(Rcpp::as<matMap >(X_));
  const matMap z(Rcpp::as<matMap >(z_));
  const double theta_x = theta_(0);
  const double theta_z = theta_(1);
  const double gamma_x = gamma_(0);
  const double gamma_z = gamma_(1);
  
  const rowVector mean_x = X.colwise().mean();
  const rowVector mean_z = z.colwise().mean();
    
  matrix A = X.rowwise() - mean_x;
  matrix B = z.rowwise() - mean_z;
  
  if(calc_covariance) {
    const matrix covX = covariance_kern(X);
    const matrix covZ = covariance_kern(z);
    
    const matrix temp_A = covX.llt().matrixL().solve(A.transpose());
    A = temp_A.transpose(); 
    
    B = 1.0/std::sqrt(covZ(0,0)) * B.array();
  }
  
  matrix kernel_matrix_X(N,N);
  matrix kernel_matrix_Z(N,N);
  
  kernel_calculation(A, kernel_matrix_X, d, theta_x, gamma_x);
  kernel_calculation(B, kernel_matrix_Z, d, theta_z, gamma_z);
  
  matrix kernel_matrix = kernel_matrix_X.array() * kernel_matrix_Z.array();
  
  return Rcpp::wrap(kernel_matrix);
}
