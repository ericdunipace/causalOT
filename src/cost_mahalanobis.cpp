#include <RcppEigen.h>

//[[Rcpp::depends(RcppEigen)]]
typedef Eigen::VectorXd vector;
typedef Eigen::Matrix<long double, Eigen::Dynamic,  1> vectorLD;
typedef Eigen::VectorXi vectorI;

typedef Eigen::MatrixXd matrix;
typedef Eigen::MatrixXi matrixI;
typedef Eigen::Matrix<long double, Eigen::Dynamic,  Eigen::Dynamic> matrixLD;

typedef Eigen::Ref<matrix> refMat;
typedef Eigen::Ref<matrixI> refMatI;

typedef Eigen::Ref<Eigen::ArrayXi> refArrayI;
typedef Eigen::Ref<Eigen::ArrayXd> refArray;

typedef Eigen::Ref<vector> refVec;
typedef Eigen::Ref<vectorI> refVecI;

typedef Eigen::Ref<const matrix> refMatConst;
typedef Eigen::Ref<const matrixI> refMatConstI;
typedef Eigen::Ref<const Eigen::ArrayXd> refArrayConst;
typedef Eigen::Ref<const Eigen::ArrayXi> refArrayConstI;

typedef Eigen::Ref<const vector> refVecConst;
typedef Eigen::Ref<const vectorI> refVecConstI;

typedef matrix::ColXpr ColXpr;
typedef matrixI::ColXpr ColXprI;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rowMat;
typedef Eigen::LLT<matrix> llt;
typedef Eigen::LDLT<matrix> ldlt;

typedef Eigen::Map<matrix> matMap;
typedef Eigen::Map<const matrix> matMapConst;
typedef Eigen::Map<Eigen::Matrix<long double, Eigen::Dynamic,  Eigen::Dynamic>> matMapLD;
typedef Eigen::Map<rowMat> rowMatMap;

typedef Eigen::Map<Eigen::VectorXd> vecMap;
typedef Eigen::Map<vectorLD> vecMapLD;
typedef Eigen::Map<const vector> vecMapConst;
typedef Eigen::Map<Eigen::VectorXi> vecMapI;
typedef Eigen::Map<const Eigen::VectorXi> vecMapConstI;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;

matrix covariance(const refMatConst & samples) {
  int S = samples.cols();
  int d = samples.rows();
  matrix c_samples(d, S);
  vector mean = samples.rowwise().mean();
  
  if(d != mean.rows()) Rcpp::stop("Dimension of mean vector not match dimension of samples vector!");
  for(int i = 0 ; i < S; i++){
    c_samples.col(i) = samples.col(i) - mean;
  }
  return matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(c_samples, 1.0/double(S-1));
}

void cost_mahal_Lp(const refMatConst & A, const refMatConst & B, 
                        const matrix & L, matrix & cost_matrix, double p) {
  double p_inv = 1.0/p;
  
  for (int j = 0; j < B.cols(); j++) { 
    vector bvec = B.col(j);
    for (int i = 0; i < A.cols(); i++) {
      double cost_p = (L.triangularView<Eigen::Lower>().solve(A.col(i)-bvec)).array().pow(p).sum();
      cost_matrix(i,j) = std::pow(cost_p, p_inv);
    }
  }
}

void cost_mahal_L2(const refMatConst & A, const refMatConst & B, 
                         const matrix & L, matrix & cost_matrix) {
  for (int j = 0; j < B.cols(); j++) { 
    vector bvec = B.col(j);
    for (int i = 0; i < A.cols(); i++) {
      cost_matrix(i,j) = (L.triangularView<Eigen::Lower>().solve(A.col(i)-bvec)).norm();
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
      cost_matrix(i,j) = (L.triangularView<Eigen::Lower>().solve(A.col(i)-bvec)).cwiseAbs().sum();
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
  const matrix L = cov.llt().matrixL();
  // Rcpp::Rcout << covA(0,0)<<std::endl;
  // Rcpp::Rcout << cov(0,0);
  
  matrix cost_matrix(N,M);
  
  if(p == 2.0) {
    cost_mahal_L2(A, B, L, cost_matrix);
  } else if (p == 1.0){
    cost_mahal_L1(A, B, L, cost_matrix);
  } else {
    cost_mahal_Lp(A, B, L, cost_matrix, p);
  }
  
  return Rcpp::wrap(cost_matrix);
}
