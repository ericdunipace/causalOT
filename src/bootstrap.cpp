#include "causalOT_types.h"

//[[Rcpp::export]]
Rcpp::NumericVector bootStrap_(Rcpp::List & w_list, 
                               int nboot,
                               Rcpp::Environment & object
                    ) {
  // data and lengths
  int n_w = w_list.length();
  int n = int(object["n"]);
  int m = int(object["m"]);
  Rcpp::NumericVector a = object["a"];
  Rcpp::NumericVector b = object["b"];
  
  // output vector
  Rcpp::NumericVector means(n_w);
  Rcpp::IntegerVector a_counts(n); // holds boot measure counts
  Rcpp::IntegerVector b_counts(m); // holds boot measure counts
  Rcpp::NumericVector w_tilde(n);
  Rcpp::NumericVector b_tilde(m);
  
  // evaluation function from R6 object
  Rcpp::Function evalFunction = object["evalBoot"];

    // output vector
  means.fill(0.0);
  
  // denominator for mean
  double n_bd = double(nboot);
  
  for(int i = 0; i < nboot; ++i) {
    GetRNGstate();
    rmultinom(n, a.begin(), n, a_counts.begin());
    rmultinom(m, b.begin(), m, b_counts.begin());
    PutRNGstate();
    b_tilde = Rcpp::as<Rcpp::NumericVector>(b_counts)/double(m);
    
    for (int j = 0; j < n_w; ++j) {
      w_tilde = Rcpp::as<Rcpp::NumericVector>(a_counts) *
        Rcpp::as<Rcpp::NumericVector>(w_list[j]);
      w_tilde = w_tilde /  Rcpp::sum( w_tilde );
      double evalout = Rcpp::as<double>(evalFunction(w_tilde, b_tilde));
      means(j) += evalout/n_bd;
    }
    
  }
  
  return(means);
}


//[[Rcpp::export]]
Rcpp::NumericVector sbw_oop_bs_(Rcpp::List & w_list, 
                               int nboot,
                               matrix & source,
                               vector & target,
                               Rcpp::NumericVector & a
) {
  // data and lengths
  int n_w = w_list.length();
  int n = source.rows();
  
  // output vector
  Rcpp::NumericVector means(n_w);
  Rcpp::IntegerVector a_counts(n); // holds boot measure counts
  matrix w_matrix(n,n_w);
  vector w_tilde(n);
  
  // holder for weights
  for(int j = 0; j < n_w; j++) w_matrix.col(j) = Rcpp::as<vector>(w_list[j]);
  
  // output vector
  means.fill(0.0);
  
  // denominator for mean
  double n_bd = double(nboot);
  
  for (int i = 0; i < nboot; ++i) {
    // draw bootstrap distribution
    rmultinom(n, a.begin(), n, a_counts.begin());
    vector a_count_dbl = Rcpp::as<vector>(a_counts);
    
    // eval each weight
    for (int j = 0; j < n_w; ++j) {
      w_tilde = a_count_dbl.array() * w_matrix.col(j).array();
      double evalout = ((source.transpose() * w_tilde) - target).array().abs().mean();
      means(j) += evalout/n_bd;
    }
    
  }
  
  return(means);
}
