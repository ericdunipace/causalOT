#include "causalOT_types.h"

const double log_2pi = std::log(2.0) + std::log(3.141592653589793115997963468544185161590576171875);

void kernel_calculation(const matrix & A, matrix & kernel_matrix, double p,
                           double theta, double gamma, const std::string & kernel_) {
  //A can be demeaned and scaled by covariance prior to using
  int n = A.rows();
  if(kernel_ == "polynomial") {
    kernel_matrix = matrix(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
    kernel_matrix.array() *= theta;
    kernel_matrix.array() += 1.0;
    kernel_matrix.array() = gamma * kernel_matrix.array().pow(p);
  } else if (kernel_ == "RBF") {
    for(int i = 0; i < n; i++) kernel_matrix.col(i) = (A.rowwise() - A.row(i)).matrix().rowwise().squaredNorm();
    kernel_matrix.array() *= theta * (-0.5);
    kernel_matrix.array() = gamma * kernel_matrix.array().exp();
  } else if (kernel_ == "linear") {
    kernel_matrix = matrix(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
  }
  
}

rowVector mean_kern(const refMatConst & X, const Rcpp::IntegerVector & z,
                 const std::string & estimand,
                 const bool calc_covariance) {
  
  int d = X.cols();
  int N = X.rows();
  rowVector mean_x = rowVector::Zero(d);
  if(!calc_covariance) return(mean_x);
  // if(kernel_type == "exponential") {
  //   
  // } else if (kernel_type == "polynomial") {
  mean_x = X.colwise().mean();
    // if(estimand == "ATE") {
    //   mean_x = X.colwise().mean();
    // } else if (estimand == "ATT") {
    //   int nt = Rcpp::sum(z);
    //   
    //   for(int i = 0; i < N; i ++) {
    //     if(z(i) == 1) {
    //       mean_x += X.row(i) / double(nt);
    //     }
    //   }
    // } else if (estimand == "ATC") {
    //   int nc = N - Rcpp::sum(z);
    //   
    //   for(int i = 0; i < N; i ++) {
    //     if(z(i) == 0) {
    //       mean_x += X.row(i) / double(nc);
    //     }
    //   }
    // }
  // }
  
  return(mean_x);
}

void covariance_kern(matrix & A, 
                     const Rcpp::IntegerVector & z,
                     const std::string & estimand) {
  int N = A.rows();
  int d = A.cols();

  matrix covX(d,d);
  
  covX = matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.transpose(), 1.0/double(N-1));
  // if(estimand == "ATE") {
  //   covX = matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.transpose(), 1.0/double(N-1));
  //   
  // } else if (estimand == "ATT") {
  //   int nt = Rcpp::sum(z);
  //   matrix At(nt,d);
  //   int idx = 0;
  //   
  //   for (int i = 0; i < N; i ++) {
  //     if(z(i) == 1) {
  //       At.row(idx) = A.row(i);
  //       idx++;
  //     }
  //   }
  //   
  //   covX = matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(At.transpose(), 1.0/double(nt-1));
  //   
  // } else if (estimand == "ATC") {
  //   int nc = N - Rcpp::sum(z);
  //   matrix Ac(nc,d);
  //   int idx = 0;
  //   
  //   for(int i = 0; i < N; i ++) {
  //     if(z(i) == 0) {
  //       Ac.row(idx) = A.row(i);
  //       idx++;
  //     }
  //   }
  //   covX = matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(Ac.transpose(), 1.0/double(nc-1));
  // }
  
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
                                 const std::string & kernel_,
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
  
  kernel_calculation(A, kernel_matrix_X, p, theta_x, gamma_x, kernel_);
  kernel_calculation(B, kernel_matrix_Z, p, theta_z, gamma_z, kernel_);
  
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
                                      const std::string & kernel_,
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
  
  rowVector mean_x = mean_kern(X, z, estimand, calc_covariance);
  
  matrix A = X.rowwise() - mean_x;
  
  
  if(calc_covariance) {
    covariance_kern(A, z, estimand);
  }
  
  matrix kernel_matrix_X(N,N);
  matrix kernel_matrix = sigma_2.asDiagonal();
  vector mean_kernel(N);
  
  if(kernel_ == "RBF" || kernel_ == "polynomial") {
    matrix theta = matrix::Zero(N,N);
    matrix gamma = matrix::Zero(N,N);
    
    for(int i = 0; i < N; i ++) {
      for(int j = 0; j < N; j ++) {
        if( (z(i) ==1) & (z(j) == 1) ) {
          theta(i,j) = theta_1;
          gamma(i,j) = gamma_1;
        }
        if( ( z(i) == 0 ) & ( z(j) == 0 ) ) {
          theta(i,j) = theta_0;
          gamma(i,j) = gamma_0;
        }
      }
    }
    if(kernel_ == "RBF") {
      for(int n = 0; n < N; n ++) kernel_matrix_X.col(n) = (A.rowwise() - A.row(n)).matrix().rowwise().squaredNorm();
      kernel_matrix.array() += (- 0.5 * kernel_matrix_X.array() *theta.array() ).exp() * gamma.array();
      for (int n = 0; n < N; n ++) {
        if (z(n) == 0) {
          mean_kernel(n) = ((-0.5 * theta_0 * kernel_matrix_X.col(n).array()).exp() * gamma_0).mean();
        } else if ( z(n) == 1) {
          mean_kernel(n) = ((-0.5 * theta_1 * kernel_matrix_X.col(n).array() ).exp() * gamma_1).mean();
        } else {
          Rcpp::stop("Invalid value in z!");
        }
      }
    } else if (kernel_ == "polynomial") {
      kernel_matrix_X = matrix(N, N).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
      kernel_matrix.array() += (kernel_matrix_X.array() *theta.array() + 1.0).pow(p) * gamma.array();
      for (int n = 0; n < N; n ++) {
        if (z(n) == 0) {
          mean_kernel(n) = ((theta_0 * kernel_matrix_X.col(n).array() + 1.0).pow(p) * gamma_0).mean();
        } else if ( z(n) == 1) {
          mean_kernel(n) = ((theta_1 * kernel_matrix_X.col(n).array() + 1.0).pow(p) * gamma_1).mean();
        } else {
          Rcpp::stop("Invalid value in z!");
        }
      }
  }
  } else if (kernel_ == "linear") {
    kernel_matrix = matrix(N, N).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
    mean_kernel = kernel_matrix.colwise().mean();
  }
  
  return (
      Rcpp::List::create(Rcpp::Named("cov_kernel") = kernel_matrix, 
                            Rcpp::Named("mean_kernel") = mean_kernel)
            );
}

//[[Rcpp::export]]
Rcpp::List kernel_calc_pred_(const Rcpp::NumericMatrix & X_,  //confounders
                             const Rcpp::NumericMatrix & X_test_, //test points
                              const Rcpp::IntegerVector & z,  //tx, a vector but easier if matrix
                              const double p,
                              const Rcpp::NumericVector  & theta_,
                              const Rcpp::NumericVector & gamma_,
                              const Rcpp::NumericVector & sigma_2_,
                              const std::string & kernel_,
                              const bool calc_covariance,
                              const std::string & estimand) {
  
  int N = X_.rows();
  int N1 = sum(z);
  int N0 = N - N1;
  int Ntest = X_test_.rows();
  int d = X_.cols();
  if(N != z.size() ) Rcpp::stop("Observations for X and treatement indicator must be equal!");
  // if(N != sigma_2_.size() ) Rcpp::stop("Observations for X and variances must be equal!");
  
  const matMap X(Rcpp::as<matMap >(X_));
  const matMap X_test(Rcpp::as<matMap >(X_test_));
  // const vecMap sigma_2(Rcpp::as<vecMap>(sigma_2_));
  // const matMap z(Rcpp::as<matMap >(z_));
  const double theta_0 = theta_(0);
  const double theta_1 = theta_(1);
  const double gamma_0 = gamma_(0);
  const double gamma_1 = gamma_(1);
  vector sigma_1(N1);
  vector sigma_0(N0);
  
  
  matrix A = X;
  matrix Atest = X_test.transpose();
  
  if(calc_covariance) {
    rowVector mean_x = mean_kern(X, z, estimand, calc_covariance);
    
    A = A.rowwise() - mean_x;
    Atest = Atest.colwise() - mean_x.transpose();
    
    llt chol = covariance_centered(A).llt();
    Atest.noalias() = chol.matrixL().solve(Atest);
    matrix A_temp = chol.matrixL().solve(A.transpose());
    A = A_temp.transpose();
  }
  sigma_0.fill(sigma_2_(0));
  sigma_1.fill(sigma_2_(1));
  matrix kernel_matrix_X0(N0,N0);
  matrix kernel_matrix_X1(N1,N1);
  matrix cross_matrix_X0(N0,Ntest);
  matrix cross_matrix_X1(N1,Ntest);
  matrix kernel_matrix_0 = sigma_0.asDiagonal();
  matrix kernel_matrix_1 = sigma_1.asDiagonal();
  vector mean_kernel(N);
  
  int count1 = 0;
  int count0 = 0;
  matrix A1(d,N1);
  matrix A0(d,N0);
  for(int n = 0; n < N; n ++) {
    if(z(n) == 1) {
      A1.col(count1) = A.row(n);
      count1++;
    } else if (z(n) == 0) {
      A0.col(count0) = A.row(n);
      count0++;
    } else {
      Rcpp::stop("Invalid value in z!");
    }
  }
  
  if(kernel_ == "RBF") {
    
    for(int n = 0; n < N1; n++){
      kernel_matrix_X1.col(n) = (A1.colwise() - A1.col(n)).matrix().colwise().squaredNorm();
    }
    for(int n = 0; n < N0; n++){
      kernel_matrix_X0.col(n) = (A0.colwise() - A0.col(n)).matrix().colwise().squaredNorm();
    }
    kernel_matrix_1.array() += (- 0.5 * kernel_matrix_X1.array() *theta_1 ).exp() * gamma_1;
    kernel_matrix_0.array() += (- 0.5 * kernel_matrix_X0.array() *theta_0 ).exp() * gamma_0;
    for(int n = 0; n < Ntest; n ++) {
      cross_matrix_X1.col(n) = (A1.colwise() - Atest.col(n)).matrix().colwise().squaredNorm();
      cross_matrix_X0.col(n) = (A0.colwise() - Atest.col(n)).matrix().colwise().squaredNorm();
    }
    cross_matrix_X0.array() = (- 0.5 * cross_matrix_X0.array() *theta_0 ).exp() * gamma_0;
    cross_matrix_X1.array() = (- 0.5 * cross_matrix_X1.array() *theta_1 ).exp() * gamma_1;
  } else if (kernel_ == "polynomial") {
    kernel_matrix_X0 = matrix(N0, N0).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A0.transpose());
    kernel_matrix_X1 = matrix(N1, N1).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A1.transpose());
    kernel_matrix_1.array() += (kernel_matrix_X1.array() *theta_1 + 1.0).pow(p) * gamma_1;
    kernel_matrix_0.array() += (kernel_matrix_X0.array() *theta_0 + 1.0).pow(p) * gamma_0;
    cross_matrix_X1 = ((A1.transpose() * Atest).array() *theta_1 + 1.0).pow(p) * gamma_1;
    cross_matrix_X0 = ((A0.transpose() * Atest).array() *theta_0 + 1.0).pow(p) * gamma_0;
  } else if (kernel_ == "linear") {
    kernel_matrix_X0 = matrix(N0, N0).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A0.transpose());
    kernel_matrix_X1 = matrix(N1, N1).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A1.transpose());
    cross_matrix_X1 = A1.transpose() * Atest;
    cross_matrix_X0 = A0.transpose() * Atest;
  } else {
    Rcpp::stop("Current kernel choice not supported");
  }
  
  return (
      Rcpp::List::create(Rcpp::Named("0") = 
        Rcpp::List::create(Rcpp::Named("cov") = kernel_matrix_0,
                           Rcpp::Named("cross") = cross_matrix_X0), 
       Rcpp::Named("1") = 
         Rcpp::List::create(Rcpp::Named("cov") = kernel_matrix_1,
                            Rcpp::Named("cross") = cross_matrix_X1)
       )
  );
}

//[[Rcpp::export]]
Rcpp::NumericMatrix kernel_update_(const Rcpp::NumericMatrix & sim_,  //similarity matrix
                                 const Rcpp::IntegerVector & z_,  //tx
                                 const double p,
                                 const Rcpp::NumericVector  & theta_,
                                 const Rcpp::NumericVector & gamma_,
                                 const Rcpp::NumericVector & sigma_2_,
                                 const std::string & kernel_) 
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
      if( (z_(i) == 1 ) & ( z_(j) == 1 ) ) {
        theta(i,j) = theta_1;
        gamma(i,j) = gamma_1;
      }
      if( (z_(i) == 0 ) & ( z_(j) == 0) ) {
        theta(i,j) = theta_0;
        gamma(i,j) = gamma_0;
      }
    }
  }
  
  matrix kernel_matrix = sigma_2.asDiagonal();

  if(kernel_== "RBF") {
    kernel_matrix.array() += (sim.array() * theta.array() * -0.5).exp() * gamma.array();
  } else if (kernel_== "polynomial") {
    kernel_matrix.array() += (sim.array() * theta.array() + 1.0).pow(p) * gamma.array();
  } else if (kernel_== "linear" ) {
    kernel_matrix += sim;
  }
  
  return Rcpp::wrap(kernel_matrix);
}

//[[Rcpp::export]]
matrix similarity_calc_(const Rcpp::NumericMatrix & X_,  //confounders
                        const Rcpp::IntegerVector & z,
                                 const bool calc_covariance,
                                 const std::string & estimand) { 
  
  int N = X_.rows();

  const matMap X(Rcpp::as<matMap >(X_));

  matrix A = X;

  if(calc_covariance) {
    const rowVector mean_x = mean_kern(X, z, estimand, calc_covariance);
    
    A = X.rowwise() - mean_x;
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
                                 const std::string & kernel_,
                                 const bool calc_covariance,
                                 const std::string & estimand                             ) {
  
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
  
  rowVector mean_x = mean_kern(X, z, estimand, calc_covariance);//, kernel_type);
  
  matrix A = X.rowwise() - mean_x;
  
  
  if(calc_covariance) {
    covariance_kern(A, z, estimand);
  }
  
  
  matrix kernel_matrix_X(N,N);
  
  if( (kernel_== "polynomial") | (kernel_=="linear") ) {
    kernel_matrix_X = matrix(N, N).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
  } else if (kernel_ == "RBF") {
    for(int n = 0; n < N; n ++) kernel_matrix_X.col(n) = (A.rowwise() - A.row(n)).matrix().rowwise().squaredNorm();
  }
  
  // Block of size (p,q), starting at (i,j)	
  // matrix.block(i,j,p,q);
  // matrix.block<p,q>(i,j);
  
  if(estimand == "ATT") {
    matrix temp_sim_t = kernel_matrix_X.block(0,N_c,N_c,N_t);
    matrix kernel_matrix_t(N_c, N_t);
    if(kernel_== "polynomial") {
      kernel_matrix_t = (temp_sim_t.array() * theta_0 + 1.0).pow(p) * gamma_0;
    } else if (kernel_== "RBF") {
      kernel_matrix_t = (temp_sim_t.array() * theta_0 * -0.5).exp() * gamma_0;
    } else if (kernel_ == "linear") {
      kernel_matrix_t = temp_sim_t;
    }
    
    output[0] = Rcpp::wrap(kernel_matrix_t);
    output[1] = R_NilValue;
  } else if (estimand == "ATC") {
    matrix temp_sim_c = kernel_matrix_X.block(0,N_c,N_c,N_t);
    matrix kernel_matrix_c(N_c, N_t);
    if(kernel_== "polynomial") {
      kernel_matrix_c = (temp_sim_c.array() * theta_1 + 1.0).pow(p) * gamma_1;
    } else if (kernel_== "RBF") {
      kernel_matrix_c = (temp_sim_c.array() * theta_1 * -0.5).exp() * gamma_1;
    } else if (kernel_ == "linear") {
      kernel_matrix_c = temp_sim_c;
    }
    output[0] = Rcpp::wrap(kernel_matrix_c);
    output[1] = R_NilValue;
  } else if (estimand == "ATE") {
    matrix temp_sim_c = kernel_matrix_X.block(0,0,N_c,N);
    matrix temp_sim_t = kernel_matrix_X.block(N_c,0,N_t,N);
    // Rcpp::Rcout << temp_sim_c(0,0) << ", ";
    matrix kernel_matrix_c(N_c, N);
    matrix kernel_matrix_t(N_t, N);
    
    if(kernel_== "polynomial") {
      kernel_matrix_c = (temp_sim_c.array() * theta_1 + 1.0).pow(p) * gamma_1;
      kernel_matrix_t = (temp_sim_t.array() * theta_0 + 1.0).pow(p) * gamma_0;
    } else if (kernel_== "RBF") {
      kernel_matrix_c = (temp_sim_c.array() * theta_1 * -0.5).exp() * gamma_1;
      kernel_matrix_t = (temp_sim_t.array() * theta_0 * -0.5).exp() * gamma_0;
    } else if (kernel_== "linear") {
      kernel_matrix_c = temp_sim_c;
      kernel_matrix_t = temp_sim_t;
    }
    // matrix kernel_matrix_c = temp_sim_c;
    // matrix kernel_matrix_t = temp_sim_t;
    // Rcpp::Rcout << kernel_matrix_c(0,0);
    output[0] = Rcpp::wrap(kernel_matrix_c);
    output[1] = Rcpp::wrap(kernel_matrix_t);
    
  }
  
  return Rcpp::wrap(output);
  
}
