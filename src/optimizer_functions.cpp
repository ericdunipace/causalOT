#include "causalOT_types.h"


// [[Rcpp::export]]
double cotDualL2_obj_(const SEXP & vars_, const SEXP & QQ,
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
Rcpp::NumericVector cotDualL2_grad_(const SEXP & vars_, 
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


vector cot_d_omega_entropy(vector &  eta, double lambda) {
  // return((eta.array()/lambda).exp());
  
  double max_val = eta.maxCoeff()/lambda;
  vector out = eta.array()/lambda - max_val;
  double l_total = std::log(out.array().exp().sum()) + max_val;
  
  return((out.array() + max_val - l_total).exp());
}

vector cot_d_omega_L2(vector &  eta, double lambda) {
  vector out(eta.rows());
  
  for(int i = 0; i < eta.rows(); i ++) {
    if(eta(i) > 0.) {
      out(i) = eta(i) / lambda;
    } else {
      out(i) = 0.0;
    }
  }
  
  return(out);
}

// vector cot_omega_entropy(vector &  eta, double lambda) {
//   return(lambda * (eta.array()/lambda).exp());
// }
double cot_omega_entropy(vector &  eta, double lambda) {
  
  double max_val = eta.maxCoeff()/lambda;
  vector out = eta.array()/lambda - max_val;
  double l_total = std::log(out.array().exp().sum()) + max_val;
  
  return(l_total * lambda);
}

// vector cot_omega_L2(vector &  eta, double lambda) {
//   
//   vector out(eta.rows());
//   
//   for(int i = 0; i < eta.rows(); i ++) {
//     if(eta(i) > 0.) {
//       double eta_d = eta(i);
//       out(i) = eta_d * eta_d * 0.5 / lambda;
//     } else {
//       out(i) = 0.0;
//     }
//   }
//   
//   return(out);
// }

double cot_omega_L2(vector &  eta, double lambda) {
  
  double out = 0.0;
  
  for(int i = 0; i < eta.rows(); i ++) {
    if ( eta(i) > 0.0 ) {
      double eta_d = eta(i);
      out += eta_d * eta_d * 0.5 / lambda;
    }
  }
  
  return(out);
}

// template <typename T> double sgn(T val) { // from https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
//   return (T(0) < val) - (val < T(0));
// }

Rcpp::XPtr<dblomegaPtr> getPtrOmega (std::string penalty) {
  if (penalty == "entropy") {
    return(Rcpp::XPtr<dblomegaPtr>(new dblomegaPtr(&cot_omega_entropy)) );
  } else if (penalty == "L2") {
    return(Rcpp::XPtr<dblomegaPtr>(new dblomegaPtr(&cot_omega_L2)) );
  } else {
    Rcpp::stop("Penalty function not found!");
  }  
}

Rcpp::XPtr<omegaPtr> getPtrOmega_d (std::string penalty) {
  if (penalty == "entropy") {
    return(Rcpp::XPtr<omegaPtr>(new omegaPtr(&cot_d_omega_entropy)) );
  } else if (penalty == "L2") {
    return ( Rcpp::XPtr<omegaPtr>(new omegaPtr(&cot_d_omega_L2)) );
  } else {
    Rcpp::stop("Penalty function not found!");
  }  
}


double sbw_omega_entropy(vector &  eta, double lambda) {
  int N = eta.rows();
  double l_a = -std::log(double(N));
  
  vector out = -eta.array() + l_a;
  double max_val = out.maxCoeff();
  out.array() -= max_val;
  double l_total = std::log(out.array().exp().sum()) + max_val;
  
  return(l_total);
}


double sbw_omega_L2(vector &  eta, double lambda) {
  int N = eta.rows();
  double ndiv = 1.0/double(N);
  
  double out = 0.0;
  
  for(int i = 0; i < eta.rows(); i ++) {
    if(eta(i) > 0.) {
      double eta_d = eta(i);
      out += eta_d * eta_d * 0.25 - eta_d * ndiv;
    } 
  }
  
  return(out);
}


vector sbw_d_omega_entropy(vector &  eta, double lambda) {
  int N = eta.rows();
  double l_a = -std::log(double(N));
  
  vector out = -eta.array() + l_a;
  double max_val = out.maxCoeff();
  out.array() -= max_val;
  double l_total = std::log(out.array().exp().sum()) + max_val;
  out.array() += max_val - l_total;
  
  return(out.array().exp());
}


vector sbw_d_omega_L2(vector &  eta, double lambda) {
  int N = eta.rows();
  double ndiv = 1.0/double(N);
  
  vector out(N);
  
  for(int i = 0; i < eta.rows(); i ++) {
    if(eta(i) > 0.) {
      out(i) = eta(i) * 0.5 - ndiv;
    } else {
      out(i) = 0.0;
    }
  }
  
  return(out);
}

// template <typename T> double sgn(T val) { // from https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
//   return (T(0) < val) - (val < T(0));
// }

Rcpp::XPtr<dblomegaPtr> getPtrSbwOmega (std::string penalty) {
  if (penalty == "entropy") {
    return(Rcpp::XPtr<dblomegaPtr>(new dblomegaPtr(&sbw_omega_entropy)) );
  } else if (penalty == "L2") {
    return(Rcpp::XPtr<dblomegaPtr>(new dblomegaPtr(&sbw_omega_L2)) );
  } else {
    Rcpp::stop("Penalty function not found!");
  }  
}

Rcpp::XPtr<omegaPtr> getPtrSbwOmega_d (std::string penalty) {
  if (penalty == "entropy") {
    return(Rcpp::XPtr<omegaPtr>(new omegaPtr(&sbw_d_omega_entropy)) );
  } else if (penalty == "L2") {
    return ( Rcpp::XPtr<omegaPtr>(new omegaPtr(&sbw_d_omega_L2)) );
  } else {
    Rcpp::stop("Penalty function not found!");
  }  
}


// [[Rcpp::export]]
double cotDual_obj_(const SEXP & vars_, const SEXP & QQ,
                      const SEXP& cost_, const SEXP & pf_,
                      double lambda, 
                      std::string penalty) {
  // map R types to Eigen types
  const vecMap vars(Rcpp::as<vecMap>(vars_));
  const Eigen::MappedSparseMatrix<double> Q(Rcpp::as< Eigen::MappedSparseMatrix<double> >(QQ));
  const vecMap cost(Rcpp::as<vecMap>(cost_));
  const vecMap pf(Rcpp::as<vecMap>(pf_));
  
  // get omega function
  Rcpp::XPtr<dblomegaPtr> xpfun = getPtrOmega(penalty);
  dblomegaPtr omega_fun = *xpfun;
  
  // calculate omega gradient
  vector eta = Q * vars - cost; 
  double obj =  pf.dot(vars) - omega_fun(eta, lambda);
  
  return(-obj); // negative b/c lbfgs minimizes by default and we need to maximize
}

// [[Rcpp::export]]
Rcpp::NumericVector cotDual_grad_(const SEXP & vars_, 
                                    const SEXP & QQ,
                                    const SEXP& cost_, const SEXP & pf_,
                                    double lambda,
                                    std::string penalty) {
  // map R types to Eigen types
  const vecMap vars(Rcpp::as<vecMap>(vars_));
  const Eigen::MappedSparseMatrix<double> Q(Rcpp::as< Eigen::MappedSparseMatrix<double> >(QQ));
  const vecMap cost(Rcpp::as<vecMap>(cost_));
  const vecMap pf(Rcpp::as<vecMap>(pf_));
  
  // get omega function
  Rcpp::XPtr<omegaPtr> xpfun_d = getPtrOmega_d(penalty);
  omegaPtr d_omega_fun = *xpfun_d;
  
  // calculate  gradient
  vector eta = Q * vars - cost; 
  vector grad = -Q.transpose() * d_omega_fun(eta, lambda);
  grad += pf;
  
  // convert to R vector
  Rcpp::NumericVector output = Rcpp::wrap(-grad); 
  
  return(output); // negative b/c lbfgs minimizes by default and we need to maximize
}


// [[Rcpp::export]]
double sbwDual_obj_(const SEXP & vars_, const SEXP & QQ,
                    const SEXP & pf_,
                    std::string penalty) {
  // map R types to Eigen types
  const vecMap vars(Rcpp::as<vecMap>(vars_));
  const matMap Q(Rcpp::as< matMap >(QQ));
  const vecMap pf(Rcpp::as<vecMap>(pf_));
  
  // get omega function
  Rcpp::XPtr<dblomegaPtr> xpfun = getPtrSbwOmega(penalty);
  dblomegaPtr omega_fun = *xpfun;
  
  // calculate omega gradient
  vector eta = Q * vars; 
  double obj =  pf.dot(vars) - omega_fun(eta, 1.0);
  
  return(-obj); // negative b/c lbfgs minimizes by default and we need to maximize
}

// [[Rcpp::export]]
Rcpp::NumericVector sbwDual_grad_(const SEXP & vars_, 
                                  const SEXP & QQ,
                                  const SEXP & pf_,
                                  std::string penalty) {
  // map R types to Eigen types
  const vecMap vars(Rcpp::as<vecMap>(vars_));
  const matMap Q(Rcpp::as< matMap >(QQ));
  const vecMap pf(Rcpp::as<vecMap>(pf_));
  
  // get omega function
  Rcpp::XPtr<omegaPtr> xpfun_d = getPtrOmega_d(penalty);
  omegaPtr d_omega_fun = *xpfun_d;
  
  // calculate omega gradient
  vector eta = Q * vars;
  vector grad = -Q.transpose() * d_omega_fun(eta, 1.0) ;
  grad += pf;
  
  // convert to R vector
  Rcpp::NumericVector output = Rcpp::wrap(-grad); 

  
  return(output); // negative b/c lbfgs minimizes by default and we need to maximize
}


// [[Rcpp::export]]
double otDualL2_obj_(const SEXP & vars_, 
                     const SEXP & a_,
                     const SEXP & b_,
                        const SEXP& cost_,
                        double lambda) {
  // map R types to Eigen types
  const vecMap a(Rcpp::as<vecMap>(a_));
  const vecMap b(Rcpp::as<vecMap>(b_));
  const vecMap vars(Rcpp::as<vecMap>(vars_));
  const matMap cost(Rcpp::as<matMap>(cost_));
  
  // dimensions
  int N = a.rows();
  int M = b.rows();
  
  // pull out relevant portions of vars
  const double * f = &vars[0];
  const double * g = &vars[N];
  
  // for (int i = 0; i < N; i ++) f[i] = &vars[i];
  // for (int j = 0; j < M; j ++) g[j] = &vars[j+N];
  
  double obj = 0.0;
  
  for(int i = 0; i < N; i ++) obj += *(f+i) * a(i);
  for(int j = 0; j < M; j ++) obj += *(g+j) * b(j);
  
  for (int j = 0; j < M; j ++) {
    double cur_g = *(g +j);
    for (int i = 0; i < N; i ++) {
      double eta = *(f + i) + cur_g - cost(i,j);
      if (eta > 0.0) obj -= eta * eta * 0.5 / lambda;
    }
  }
  
  // for (int i = 0; i < N; i ++) f(i) = vars(i);
  // for (int j = 0; j < M; j ++) g(j) = vars(j+N);
  // 
  
  // // main part of objective
  // double obj = f.dot(a) + g.dot(b);
  // 
  // // calculate portion from dual penalty
  // for(int j = 0; j < M; j ++) {
  //   double cur_g = g(j);
  //   for(int i = 0; i < N; i ++) {
  //     double eta = f(i) + cur_g - cost(i,j);
  //     if (eta > 0.0) obj -= eta * eta * 0.5 / lambda;
  //     }
  //   }
  // 
  delete[] f;
  delete[] g;
  return(-obj); // negative b/c lbfgs minimizes by default and we need to maximize
}

// [[Rcpp::export]]
Rcpp::NumericVector otDualL2_grad_(const SEXP & vars_, 
                                     const SEXP & a_,
                                     const SEXP & b_,
                                     const SEXP& cost_,
                                     double lambda) {
  // map R types to Eigen types
  const vecMap a(Rcpp::as<vecMap>(a_));
  const vecMap b(Rcpp::as<vecMap>(b_));
  const vecMap vars(Rcpp::as<vecMap>(vars_));
  const matMap cost(Rcpp::as<matMap>(cost_));
  
  // dimensions
  int N = a.rows();
  int M = b.rows();
  
  // pull out relevant portions of vars
  // const double *f[N];
  // const double *g[M];
  const double * f = &vars[0];
  const double * g = &vars[N];
  
  // for (int i = 0; i < N; i ++) *f[i] = &vars[i];
  // for (int j = 0; j < M; j ++) *g[j] = &vars[j+N];
  
  // save first part of gradient
  vector grad(N+M);
  grad.block(0,0,N,1) = a;
  grad.block(N,0,M,1) = b;
  
  
  // add main objective derivative to gradient
  for(int j = 0; j < M; j ++) {
    double cur_g = *(g+j);
    for(int i = 0; i < N; i ++) {
      double eta = *(f + i) + cur_g - cost(i,j);
      if (eta > 0.0) {
        double delta = eta / lambda;
        grad(i) -= delta;
        grad(j+N) -= delta;
      }
    }
  }
  delete f;
  delete g;
  
  // // pull out relevant portions of vars
  // vector f(N);
  // vector g(M);
  // 
  // for (int i = 0; i < N; i ++) (&f)[i] = (&vars)[i];
  // for (int j = 0; j < M; j ++) (&g)[j] = (&vars)[j+N];
  // 
  // // save first part of gradient
  // vector f_grad = a;
  // vector g_grad = b;
  // 
  // // add main objective derivative to gradient
  // for(int j = 0; j < M; j ++) {
  //   double cur_g = g(j);
  //   for(int i = 0; i < N; i ++) {
  //     double eta = f(i) + cur_g - cost(i,j);
  //     if (eta > 0.0) {
  //       double delta = eta / lambda;
  //       f_grad(i) -= delta;
  //       g_grad(j) -= delta;
  //     }
  //   }
  // }
  // 
  // // concatenate the gradients
  // vector grad(N+M);
  // for (int i = 0; i < N; i ++) grad(i) = f_grad(i);
  // for (int j = 0; j < M; j ++) grad(j+N) = g_grad(j);
  
  return(Rcpp::wrap(-grad)); // negative b/c lbfgs minimizes by default and we need to maximize
}


// [[Rcpp::export]]
double otDualL2_obj_self_(const SEXP & vars_, 
                       const SEXP & a_,
                       const SEXP& cost_,
                       double lambda) {
  // map R types to Eigen types
  const vecMap a(Rcpp::as<vecMap>(a_));
  const vecMap f(Rcpp::as<vecMap>(vars_));
  const matMap cost(Rcpp::as<matMap>(cost_));
  
  // dimensions
  int N = a.rows();

  
  // main part of objective
  double obj = 2 * f.dot(a);
  
  // calculate portion from dual penalty
  for(int j = 0; j < N; j ++) {
    double cur_g = f(j);
    for(int i = 0; i < N; i ++) {
      double eta = f(i) + cur_g - cost(i,j);
      if (eta > 0.0) obj -= eta * eta * 0.5 / lambda;
    }
  }
  
  return(-obj); // negative b/c lbfgs minimizes by default and we need to maximize
}


// [[Rcpp::export]]
Rcpp::NumericVector otDualL2_grad_self_ (const SEXP & vars_, 
                                     const SEXP & a_,
                                     const SEXP& cost_,
                                     double lambda) {
  // map R types to Eigen types
  const vecMap a(Rcpp::as<vecMap>(a_));
  const vecMap f(Rcpp::as<vecMap>(vars_));
  const matMap cost(Rcpp::as<matMap>(cost_));
  
  // dimensions
  int N = a.rows();
  
  
  // save first part of gradient
  vector f_grad = a.array() * 2;

  // add main objective derivative to gradient
  for(int j = 0; j < N; j ++) {
    double cur_g = f(j);
    for(int i = 0; i < N; i ++) {
      double eta = f(i) + cur_g - cost(i,j);
      if (eta > 0.0) {
        f_grad(i) -= 2 * eta  / lambda;
      }
    }
  }
  
  // concatenate the gradients
  
  return(Rcpp::wrap(-f_grad)); // negative b/c lbfgs minimizes by default and we need to maximize
}


// pointers



//[[Rcpp::export]]
Rcpp::XPtr<gradCrossPtr> otL2_grad_(){
  return(Rcpp::XPtr<gradCrossPtr>(new gradCrossPtr(&otDualL2_grad_)));
}



//[[Rcpp::export]]
Rcpp::XPtr<objCrossPtr> otL2_obj_(){
  return(Rcpp::XPtr<objCrossPtr>(new objCrossPtr(&otDualL2_obj_)));
}


//[[Rcpp::export]]
Rcpp::XPtr<gradSelfPtr> otL2_grad_self_(){
  return(Rcpp::XPtr<gradSelfPtr>(new gradSelfPtr(&otDualL2_grad_self_)));
}


//[[Rcpp::export]]
Rcpp::XPtr<objSelfPtr> otL2_obj_self_(){
  return(Rcpp::XPtr<objSelfPtr>(new objSelfPtr(&otDualL2_obj_self_)));
}
