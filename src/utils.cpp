#include "utils.h"

//[[Rcpp::export]]
double logSumExp(vector & x_) {
  double alpha = x_(0);
  double r = 0.0;
  int N = x_.rows();
  // for(auto x : x_) {
  for (int i = 0; i < N; ++i) {
    double x = x_(i);
    if (x <= alpha) {
      r += std::exp(x - alpha);
    } else {
      r *= std::exp(alpha - x);
      r += 1.0;
      alpha = x;
    }
  }
  return(std::log(r) + alpha);
}

double logSumExp(matrix & x_) {
  double alpha = x_(0,0);
  double r = 0.0;
  int N = x_.rows();
  int M = x_.cols();
  
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i) {
      double x = x_(i,j);
      if (x <= alpha) {
        r += std::exp(x - alpha);
      } else {
        r *= std::exp(alpha - x);
        r += 1.0;
        alpha = x;
      }
    }
  }
  return(std::log(r) + alpha);
}

//[[Rcpp::export]]
vector rowLogSumExp(matrix & x_) {
  int N = x_.rows();
  int M = x_.cols();
  
  vector alpha_ = x_.col(0);
  vector r_ = vector::Zero(N);
  
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i) {
      double x = x_(i,j);
      double alpha = alpha_(i);
      double r = r_(i); 
      
      if (x <= alpha) {
        r += std::exp(x - alpha);
      } else {
        r *= std::exp(alpha - x);
        r += 1.0;
        alpha_(i) = x;
      }
      r_(i) = r;
    }
  }
  vector out = r_.array().log() + alpha_.array();
  return(out);
}

//[[Rcpp::export]]
vector colLogSumExp(matrix & x_) {
  int N = x_.rows();
  int M = x_.cols();
  
  vector out = vector::Zero(M);
  
  for (int j = 0; j < M; ++j) {
    double alpha = x_(0,j);
    double r = 1.0; 
    for (int i = 1; i < N; ++i) {
      double x = x_(i,j);
      if (x <= alpha) {
        r += std::exp(x - alpha);
      } else {
        r *= std::exp(alpha - x);
        r += 1.0;
        alpha = x;
      }
    }
    out(j) = std::log(r) + alpha;
  }
  return(out);
}

