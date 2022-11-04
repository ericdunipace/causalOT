#include <Rcpp.h>

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