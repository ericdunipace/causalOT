#include <Rcpp.h>


// [[Rcpp::export]]
void entry(SEXP & xx, Rcpp::NumericMatrix & y, int colX_) {
  Rcpp::NumericMatrix x(xx);
  int N = y.rows();
  int M = y.cols();
  int colX = colX_ - 1;
  if(colX > x.cols()) Rcpp::stop("Column out of range");
  if((colX + M )> x.cols()) Rcpp::stop("Max column out of range");
  if(N > x.rows()) Rcpp::stop("Max rows out of range");
  for(int i = 0; i < M; i ++){
    for(int j = 0; j < N; j ++) x(j, colX + i) = y(j,i);
  }
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// n <- 100
// X <- matrix(0, n,n)
// Y <- matrix(rnorm(n*n),n,n)
// all(X == 0)
// all.equal(X, Y)
// entry(X, Y, 1)
// all(X == 0)
// all.equal(X, Y)
// */
