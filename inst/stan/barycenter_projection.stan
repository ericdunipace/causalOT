functions {
  real  norm_lp_lpdf(matrix y,
               real p,
               matrix z,
               matrix gamma
               ) {
       real loss = 0.0;
       int N = rows(gamma);
       int M = cols(gamma);
       int D = cols(z);
       for(m in 1:M){
         for(n in 1:N) {
           row_vector[D] abs_diff = fabs(z[m] - y[n]);
           for(d in 1:D) loss += pow(abs_diff[d], p) * gamma[n,m];
         }
         
         
       }
       return (-loss);
  }
  
  matrix  norm_lp_m2(matrix y,
               real p,
               matrix z,
               matrix gamma
               ) {
       int N = rows(gamma);
       int M = cols(gamma);
       int D = cols(z);
       matrix[N,M] weights;
       real pm2 = p - 2.0;
       for(m in 1:M){
         for(n in 1:N) {
           row_vector[D] abs_diff = fabs(z[m] - y[n]);
           real temp_weight = 0.0;
           for(d in 1:D)  temp_weight += pow(abs_diff[d], pm2);
           weights[n,m] = gamma[n,m] * temp_weight;
         }
         
         
       }
       return weights;
  }
  
  real  wt_l2_lpdf(matrix y,
               matrix z,
               matrix gamma
               ) {
       real loss = 0.0;
       int N = rows(gamma);
       int M = cols(gamma);
       int D = cols(z);
       for(m in 1:M){
         for(n in 1:N) {
           loss += squared_distance(z[m], y[n]) * gamma[n,m];
         }
         
         
       }
       return loss;
  }
}

data {
  int<lower=0> N;
  int<lower=0> M;
  int<lower=0> D;
  matrix[N,D] y;
  matrix[N,M] gamma;
  real p;
}
parameters {
  matrix[M,D] z;
}
model {
  // matrix[N,M] weight = norm_lp_m2(y,p,z,gamma);
  
  y ~ norm_lp(p,z,gamma);
}
