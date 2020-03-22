#include "causalOT_types.h"

const double log_2pi = std::log(2.0) + std::log(3.141592653589793115997963468544185161590576171875);

void kernel_calculation(const matrix & A, matrix & kernel_matrix, double p,
                           double theta, double gamma) {
  //A must be demeaned and scaled by covariance prior to using
  int n = A.rows();
  kernel_matrix = matrix(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
  kernel_matrix.array() *= theta;
  kernel_matrix.array() += 1.0;
  kernel_matrix.array() = gamma * kernel_matrix.array().pow(p);
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
Rcpp::NumericMatrix kernel_calc_dose_(const Rcpp::NumericMatrix & X_,  //confounders
                                 const Rcpp::NumericMatrix & z_,  //tx, a vector but easier if matrix
                                 const double p,
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
  
  kernel_calculation(A, kernel_matrix_X, p, theta_x, gamma_x);
  kernel_calculation(B, kernel_matrix_Z, p, theta_z, gamma_z);
  
  matrix kernel_matrix = kernel_matrix_X.array() * kernel_matrix_Z.array();
  
  return Rcpp::wrap(kernel_matrix);
}

//[[Rcpp::export]]
Rcpp::List similarity_calc_dose_(const Rcpp::NumericMatrix & X_,  //confounders
                 const Rcpp::NumericMatrix & z_,  //tx, a vector but easier if matrix,
                 const bool calc_covariance) { 
  
  int N = X_.rows();
  if(N != z_.rows() ) Rcpp::stop("Observations for X and z must be equal!");
  
  const matMap X(Rcpp::as<matMap >(X_));
  const matMap z(Rcpp::as<matMap >(z_));
  
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
  
  matrix kernel_matrix_X = matrix(N, N).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
  matrix kernel_matrix_Z = matrix(N, N).setZero().selfadjointView<Eigen::Lower>().rankUpdate(B);

  return (Rcpp::List::create(Rcpp::Named("X") = kernel_matrix_X, Rcpp::Named("Z") = kernel_matrix_Z));
}
  
  
  //[[Rcpp::export]]
  Rcpp::NumericMatrix kernel_calc_(const Rcpp::NumericMatrix & X_,  //confounders
                                        const Rcpp::IntegerVector & z,  //tx, a vector but easier if matrix
                                        const double p,
                                        const Rcpp::NumericVector  & theta_,
                                        const Rcpp::NumericVector & gamma_,
                                        const Rcpp::NumericVector & sigma_2_,
                                        const bool calc_covariance ) {
    
    int N = X_.rows();
    if(N != z.size() ) Rcpp::stop("Observations for X and treatement indicator must be equal!");
    if(N != sigma_2_.size() ) Rcpp::stop("Observations for X and variances must be equal!");
    
    const matMap X(Rcpp::as<matMap >(X_));
    const vecMap sigma_2(Rcpp::as<vecMap>(sigma_2_));
    // const matMap z(Rcpp::as<matMap >(z_));
    const double theta_0 = theta_(0);
    const double theta_1 = theta_(1);
    const double gamma_0 = gamma_(0);
    const double gamma_1 = gamma_(1);
    matrix theta = matrix::Zero(N,N);
    matrix gamma = matrix::Zero(N,N);
    
    const rowVector mean_x = X.colwise().mean();

    matrix A = X.rowwise() - mean_x;

    if(calc_covariance) {
      const matrix covX = covariance_kern(X);

      const matrix temp_A = covX.llt().matrixL().solve(A.transpose());
      A = temp_A.transpose(); 
      
    }
    
    for(int i = 0; i < N; i ++) {
      for(int j = 0; j < N; j ++) {
        if( z(i) ==1 & z(j) == 1 ) {
          theta(i,j) = theta_1;
          gamma(i,j) = gamma_1;
        }
        if( z(i) == 0 & z(j) == 0 ) {
          theta(i,j) = theta_0;
          gamma(i,j) = gamma_0;
        }
      }
    }
    
    
    matrix kernel_matrix_X = matrix(N, N).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
   
    matrix kernel_matrix = sigma_2.asDiagonal();
    kernel_matrix.array() += (kernel_matrix_X.array() *theta.array() + 1.0).pow(p) * gamma.array();
    
    return Rcpp::wrap(kernel_matrix);
  }
//[[Rcpp::export]]
Rcpp::NumericMatrix kernel_update_(const Rcpp::NumericMatrix & sim_,  //similarity matrix
                                 const Rcpp::IntegerVector & z_,  //tx
                                 const double p,
                                 const Rcpp::NumericVector  & theta_,
                                 const Rcpp::NumericVector & gamma_,
                                 const Rcpp::NumericVector & sigma_2_) 
{
  
  int N = sim_.rows();
  if(N != z_.size() ) Rcpp::stop("Observations for X and treatement indicator must be equal!");
  if(N != sigma_2_.size() ) Rcpp::stop("Observations for X and variances must be equal!");
  
  const matMap sim(Rcpp::as<matMap >(sim_));
  const vecMap sigma_2(Rcpp::as<vecMap>(sigma_2_));
  // const matMap z(Rcpp::as<matMap >(z_));
  const double theta_0 = theta_(0);
  const double theta_1 = theta_(1);
  const double gamma_0 = gamma_(0);
  const double gamma_1 = gamma_(1);
  matrix theta = matrix::Zero(N,N);
  matrix gamma = matrix::Zero(N,N);
  
  for(int i = 0; i < N; i ++) {
    for(int j = 0; j < N; j ++) {
      if( z_(i) ==1 & z_(j) == 1 ) {
        theta(i,j) = theta_1;
        gamma(i,j) = gamma_1;
      }
      if( z_(i) == 0 & z_(j) == 0 ) {
        theta(i,j) = theta_0;
        gamma(i,j) = gamma_0;
      }
    }
  }
  
  matrix kernel_matrix = sigma_2.asDiagonal();

  kernel_matrix.array() += (sim.array() * theta.array() + 1.0).pow(p) * gamma.array();
  
  return Rcpp::wrap(kernel_matrix);
}

//[[Rcpp::export]]
matrix similarity_calc_(const Rcpp::NumericMatrix & X_,  //confounders
                                 const bool calc_covariance) { 
  
  int N = X_.rows();

  const matMap X(Rcpp::as<matMap >(X_));

  const rowVector mean_x = X.colwise().mean();

  matrix A = X.rowwise() - mean_x;

  if(calc_covariance) {
    const matrix covX = covariance_kern(X);

    const matrix temp_A = covX.llt().matrixL().solve(A.transpose());
    A = temp_A.transpose(); 
    
  }

  
  return (matrix(N, N).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A));
}


//[[Rcpp::export]]
double marginal_lik_gp_(const Rcpp::NumericVector & y_,
                        const Rcpp::NumericMatrix & K_
                        ) {
  
  const matMap K(Rcpp::as<matMap>(K_));
  const vecMap y(Rcpp::as<vecMap>(y_));
  double N = double(y_.size());
  // Eigen::JacobiSVD svd = K.jacobiSvd();
  llt chol = K.llt();
  matrix L = chol.matrixL();
  
  double log_prob = -y.dot(chol.solve(y));
  log_prob += - L.diagonal().array().log().sum() - 0.5 * N * log_2pi;
  //note, should be -0.5 * log(det(K)) but simplifies to-0.5 * 2* log(det(L)) = log(det(L)), where L is lower cholesky
  
  return(log_prob);
}


// double marginal_lik_gp_grad_(const Rcpp::NumericVector & y_,
//                         const Rcpp::NumericMatrix & K_
// ) {
//   
//   const matMap K(Rcpp::as<matMap>(K_));
//   const vecMap y(Rcpp::as<vecMap>(y_));
//   double N = double(y_.size());
//   llt L = K.llt();
//   vector alpha = L.solve(y);
//   matrix grad_mat =  alpha * alpha.transpose() - L.matrixLLT().inverse();
//   
//   
//   double log_prob = y.transpose() * L.solve(y);
//   log_prob += -0.5 * std::log(L.matrixL().determinant()) * 2.0 - 0.5 * N * std::log(2.0 * M_PI);
//   
//   return(log_prob);
// }
  