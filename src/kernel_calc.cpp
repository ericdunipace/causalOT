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

matrix mean_kern(const refMatConst & X, const Rcpp::IntegerVector & z,
                 const std::string & estimand) {
  
  int d = X.cols();
  int N = X.rows();
  
  rowVector mean_x = rowVector::Zero(d);
  
  if(estimand == "ATE") {
    mean_x = X.colwise().mean();
  } else if (estimand == "ATT") {
    int nt = Rcpp::sum(z);
    
    for(int i = 0; i < N; i ++) {
      if(z(i) == 1) {
        mean_x += X.row(i) / double(nt);
      }
    }
  } else if (estimand == "ATC") {
    int nc = N - Rcpp::sum(z);
    
    for(int i = 0; i < N; i ++) {
      if(z(i) == 0) {
        mean_x += X.row(i) / double(nc);
      }
    }
  }
  return(mean_x);
}

void covariance_kern(matrix & A, 
                     const Rcpp::IntegerVector & z,
                     const std::string & estimand) {
  int N = A.rows();
  int d = A.cols();

  matrix covX(d,d);
  
  if(estimand == "ATE") {
    covX = matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.transpose(), 1.0/double(N-1));
    
  } else if (estimand == "ATT") {
    int nt = Rcpp::sum(z);
    matrix At(nt,d);
    int idx = 0;
    
    for (int i = 0; i < N; i ++) {
      if(z(i) == 1) {
        At.row(idx) = A.row(i);
        idx++;
      }
    }
    
    covX = matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(At.transpose(), 1.0/double(nt-1));
    
  } else if (estimand == "ATC") {
    int nc = N - Rcpp::sum(z);
    matrix Ac(nc,d);
    int idx = 0;
    
    for(int i = 0; i < N; i ++) {
      if(z(i) == 0) {
        Ac.row(idx) = A.row(i);
        idx++;
      }
    }
    covX = matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(Ac.transpose(), 1.0/double(nc-1));
  }
  
  matrix temp_A = covX.llt().matrixL().solve(A.transpose());
  A = temp_A.transpose(); 
  // Rcpp::Rcout << A(0,0) << ", ";
  
}


matrix covariance_centered(matrix & A) {
  int N = A.rows();
  int d = A.cols();
  // rowVector mean = samples.colwise().mean();
  
  matrix covX = matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.transpose(), 1.0/double(N-1));
  
  return(covX);
  
}

//[[Rcpp::export]]
Rcpp::List kernel_calc_dose_(const Rcpp::NumericMatrix & X_,  //confounders
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
    const matrix covX = covariance_centered(A);
    const matrix covZ = covariance_centered(B);
    
    const matrix temp_A = covX.llt().matrixL().solve(A.transpose());
    A = temp_A.transpose(); 
    
    B = 1.0/std::sqrt(covZ(0,0)) * B.array();
  }
  
  matrix kernel_matrix_X(N,N);
  matrix kernel_matrix_Z(N,N);
  
  kernel_calculation(A, kernel_matrix_X, p, theta_x, gamma_x);
  kernel_calculation(B, kernel_matrix_Z, p, theta_z, gamma_z);
  
  matrix kernel_matrix = kernel_matrix_X.array() * kernel_matrix_Z.array();
  
  return Rcpp::List::create(Rcpp::Named("cov_kernel") = kernel_matrix, 
                            Rcpp::Named("mean_kernel") = kernel_matrix.colwise().mean());
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
    const matrix covX = covariance_centered(A);
    const matrix covZ = covariance_centered(B);
    
    const matrix temp_A = covX.llt().matrixL().solve(A.transpose());
    A = temp_A.transpose(); 
    
    B = 1.0/std::sqrt(covZ(0,0)) * B.array();
  }
  
  matrix kernel_matrix_X = matrix(N, N).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
  matrix kernel_matrix_Z = matrix(N, N).setZero().selfadjointView<Eigen::Lower>().rankUpdate(B);

  return (Rcpp::List::create(Rcpp::Named("X") = kernel_matrix_X, Rcpp::Named("Z") = kernel_matrix_Z));
}
  
  
//[[Rcpp::export]]
Rcpp::List kernel_calc_(const Rcpp::NumericMatrix & X_,  //confounders
                                      const Rcpp::IntegerVector & z,  //tx, a vector but easier if matrix
                                      const double p,
                                      const Rcpp::NumericVector  & theta_,
                                      const Rcpp::NumericVector & gamma_,
                                      const Rcpp::NumericVector & sigma_2_,
                                      const bool calc_covariance,
                                      const std::string & estimand) {
  
  int N = X_.rows();
  int d = X_.cols();
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
  
  rowVector mean_x = mean_kern(X, z, estimand);

  matrix A = X.rowwise() - mean_x;
  

  if(calc_covariance) {
    covariance_kern(A, z, estimand);
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
  // Rcpp::Rcout << kernel_matrix_X(0,0) << ", " << sigma_2(0) << "\n";
  // Rcpp::Rcout << "theta" << theta(0,0) << ", gamma " << gamma(0,0) << "\n";
  matrix kernel_matrix = sigma_2.asDiagonal();
  kernel_matrix.array() += (kernel_matrix_X.array() *theta.array() + 1.0).pow(p) * gamma.array();
  // Rcpp::Rcout << (kernel_matrix_X.array() *theta.array() + 1.0)(0,0) << " \n";
  // Rcpp::Rcout << (kernel_matrix_X.array() *theta.array() + 1.0).pow(p)(0,0) << " \n";
  // Rcpp::Rcout << kernel_matrix(0,0) << "\n";
  vector mean_kernel(N);
  for (int n = 0; n < N; n ++) {
    if (z(n) == 0) {
      mean_kernel(n) = ((theta_0 * kernel_matrix_X.col(n).array() + 1.0).pow(p) * gamma_0).mean();
    } else if ( z(n) == 1) {
      mean_kernel(n) = ((theta_1 * kernel_matrix_X.col(n).array() + 1.0).pow(p) * gamma_1).mean();
    } else {
      Rcpp::stop("Invalid value in z!");
    }
  }
  
  return (
      Rcpp::List::create(Rcpp::Named("cov_kernel") = kernel_matrix, 
                            Rcpp::Named("mean_kernel") = mean_kernel)
            );
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
                        const Rcpp::IntegerVector & z,
                                 const bool calc_covariance,
                                 const std::string & estimand) { 
  
  int N = X_.rows();

  const matMap X(Rcpp::as<matMap >(X_));

  const rowVector mean_x = mean_kern(X, z, estimand);

  matrix A = X.rowwise() - mean_x;

  if(calc_covariance) {
    covariance_kern(A, z, estimand);
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


//[[Rcpp::export]]
Rcpp::List kernel_calc_ot_(const Rcpp::NumericMatrix & X_,  //confounders
                                 const Rcpp::IntegerVector & z,  //tx, a vector but easier if matrix
                                 const double p,
                                 const Rcpp::NumericVector  & theta_,
                                 const Rcpp::NumericVector & gamma_,
                                 const bool calc_covariance,
                                 const std::string & estimand) {
  
  int N = X_.rows();
  int d = X_.cols();

  if(N != z.size() ) Rcpp::stop("Observations for X and treatement indicator must be equal!");
  // if(N != sigma_2_.size() ) Rcpp::stop("Observations for X and variances must be equal!");
  int N_t = Rcpp::sum(z);
  int N_c = N - N_t;
  
  const matMap X(Rcpp::as<matMap >(X_));
  // const vecMap sigma_2(Rcpp::as<vecMap>(sigma_2_));
  // const matMap z(Rcpp::as<matMap >(z_));
  const double theta_0 = theta_(0);
  const double theta_1 = theta_(1);
  const double gamma_0 = gamma_(0);
  const double gamma_1 = gamma_(1);
  // matrix theta = matrix::Zero(N,N);
  // matrix gamma = matrix::Zero(N,N);
  Rcpp::List output(2);
  
  rowVector mean_x = mean_kern(X, z, estimand);
  
  matrix A = X.rowwise() - mean_x;
  
  
  if(calc_covariance) {
    covariance_kern(A, z, estimand);
  }
  
  
  matrix kernel_matrix_X = matrix(N, N).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
  
  // Block of size (p,q), starting at (i,j)	
  // matrix.block(i,j,p,q);
  // matrix.block<p,q>(i,j);
  
  if(estimand == "ATT") {
    matrix temp_sim_t = kernel_matrix_X.block(0,N_c,N_c,N_t);
    matrix kernel_matrix_t = (temp_sim_t.array() * theta_0 + 1.0).pow(p) * gamma_0;
    
    output[0] = Rcpp::wrap(kernel_matrix_t);
    output[1] = R_NilValue;
  } else if (estimand == "ATC") {
    matrix temp_sim_c = kernel_matrix_X.block(0,N_c,N_c,N_t);
    matrix kernel_matrix_c = (temp_sim_c.array() * theta_1 + 1.0).pow(p) * gamma_1;
    
    output[0] = Rcpp::wrap(kernel_matrix_c);
    output[1] = R_NilValue;
  } else if (estimand == "ATE") {
    matrix temp_sim_c = kernel_matrix_X.block(0,0,N_c,N);
    matrix temp_sim_t = kernel_matrix_X.block(N_c,0,N_t,N);
    // Rcpp::Rcout << temp_sim_c(0,0) << ", ";
    matrix kernel_matrix_c = (temp_sim_c.array() * theta_0 + 1.0).pow(p) * gamma_0;
    matrix kernel_matrix_t = (temp_sim_t.array() * theta_1 + 1.0).pow(p) * gamma_1;
    // matrix kernel_matrix_c = temp_sim_c;
    // matrix kernel_matrix_t = temp_sim_t;
    // Rcpp::Rcout << kernel_matrix_c(0,0);
    output[0] = Rcpp::wrap(kernel_matrix_c);
    output[1] = Rcpp::wrap(kernel_matrix_t);
    
  }
  
  
  
  
  return Rcpp::wrap(output);
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
  