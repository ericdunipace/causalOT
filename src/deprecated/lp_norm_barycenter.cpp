#include "causalOT_types.h"

matrix weighted_mean_output(const matrix & Y,
                            const matrix & gamma) {
  
  return( Y.transpose() * gamma );
  
}

vector lp_norm_p(const ColXpr & Z, 
                 const matrix & Y, 
                 double p) {
  
  return((Y.rowwise() - Z.transpose()).array().abs().pow(p).matrix().rowwise().sum());
  
}


void Lp_weight(const matrix & Y,
                 matrix & Z,
                 const matrix & gamma,
                 matrix & weight,
                 double p,
                 int N,
                 int M,
                 int D) {
  double pm2 = p - 2.0;
  weight = gamma;
  Rcpp::Rcout <<  weight.col(0).sum()<< "\n\n";
  for(int m = 0; m < M; m++) {
    weight.col(m) *= lp_norm_p(Z.col(m), Y, pm2);
  }
  Rcpp::Rcout <<  weight.col(0).sum()<< "\n\n";
  Rcpp::Rcout <<  weight.rows() << ", " << weight.cols()<< "\n\n";
  Rcpp::Rcout << lp_norm_p(Z.col(0), Y, pm2)(0,0) << "\n";
  weight.array() = weight.array().rowwise()/(weight.colwise().sum().array());
  // Rcpp::Rcout << weight.col(0).sum() << "\n";
}

bool converged(matrix& Z_new, matrix & Z_old, double prec) {
  Rcpp::Rcout<< Z_new.isApprox(Z_old, prec) << ", ";
  return(Z_new.isApprox(Z_old, prec));
}

//[[Rcpp::export]]
Rcpp::List Lp_norm_min_(const Rcpp::NumericMatrix & Y_, 
            const double p,
            const Rcpp::NumericMatrix & gamma_,
            int iterations,
            double precision
) {
  int N = gamma_.rows();
  int M = gamma_.cols();
  int D = Y_.cols();
  if(N != Y_.rows() ) Rcpp::stop("Observations for Y and gamma must be equal!");
  
  const matMap Y(Rcpp::as<matMap >(Y_));
  const matMap gamma(Rcpp::as<matMap >(gamma_));
  matrix weight = gamma;
  matrix Z = weighted_mean_output(Y, gamma);
  matrix Z_old = matrix::Zero(D, M);
  
  int i = 0;
  
  // if(Z.rows() != D) Rcpp::stop("rows: setup Z matrix wrong!");
  // if(Z.cols() != M) Rcpp::stop("cols: setup Z matrix wrong!");
  for(i = 0; i < iterations; i++) {
    Lp_weight(Y, Z, gamma, weight, p, N, M, D);
    Z = weighted_mean_output(Y, weight);
    if(converged(Z, Z_old, precision)) {
      break;
    } else {
      Z_old = Z;
    }
  }
  Z.transposeInPlace();
  return(Rcpp::List::create(Rcpp::Named("z") = Z,
                            Rcpp::Named("iter") = i));
}