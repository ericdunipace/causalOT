functions {
  
  matrix pol_kern( matrix cost,
                        real theta,
                        real gamma,
                        real sigma,
                        real p
                        ) { //kernel
  int N = rows(cost);
  
  matrix[N, N] K;
  
  for (j in 1:N) {
    for (i in j:N) {
      K[i,j] = gamma * gamma * pow(1.0 + theta * theta * cost[i,j], p);
      K[j,i] = K[i,j];
    }
  }
  for (i in 1:N) {
    K[i,i] += sigma * sigma;
  }
  return K;
}


matrix rbf_kern( matrix cost,
                        real theta,
                        real gamma,
                        real sigma
                        ) { //kernel
  int N = rows(cost);
  
  matrix[N, N] K;
  
  for (j in 1:(N-1)) {
    for (i in (j+1):N) {
      real calc = gamma * gamma * exp(- 0.5 * 1/(theta * theta) * cost[i,j]);
      K[i,j] = calc;
      K[j,i] = calc;
    }
  }
  for (i in 1:N) {
    K[i,i] = gamma * gamma + sigma * sigma;
  }
   
  return K;
}

matrix linear_kern( matrix cost,
                        // real theta_0,
                        // real theta_1,
                        // real gamma_0,
                        // real gamma_1,
                        real sigma
                        ) { //kernel
  int N = rows(cost);
  
  matrix[N, N] K = cost;
  
  for (i in 1:N) {
    K[i,i] += sigma * sigma;
  }
   
  return K;
}
}

data {
  int<lower=0> N;
  vector[N] y;
  matrix[N,N] discrep;
  real p;
  int<lower = 1, upper = 3> kernel; //1 = polynomial, 2 = expo, 3 = linear
  real a;
  real b;
}

transformed data {
  vector[N] mu = rep_vector(0.0, N);
  int num_param = 1;
  real delta_x = sqrt(max(discrep));
  int use_gamma = 0;
  if(kernel == 3) num_param = 0;
  
  if(a != 0.0) use_gamma = 1;
}

parameters {
  real<lower=0.0> sigma_half;
  real<lower=0.0> theta_half[num_param];
  real<lower=0.0> gamma_half[num_param];
}

transformed parameters {
  matrix[N, N] L;
  {
    matrix[N,N] Sigma;
  
    if(kernel == 1 ) {
      Sigma = pol_kern(discrep, 
                      theta_half[1],
                      gamma_half[1],
                      sigma_half,
                      p);
    } else if (kernel == 2) {
      Sigma = rbf_kern(discrep, 
                      theta_half[1],
                      gamma_half[1],
                      sigma_half);
    } else if (kernel == 3) {
      Sigma = linear_kern(discrep, sigma_half);
    }
    L = cholesky_decompose(Sigma);
  }
} 

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and covariance matrix 'Sigma'.
model {
  
  sigma_half ~ normal(0.0, 0.5);
  
  if(kernel == 2) {
      if(use_gamma == 0){
        theta_half[1] ~ normal(0.0, delta_x/3.0);
      } else if (use_gamma == 1) {
        theta_half[1] ~ inv_gamma(a, b);
      }
      gamma_half[1] ~ normal(0.0, 2.5);
  }
  
  // print("before control", target());
  y ~ multi_normal_cholesky(mu, L);
  // print("before treatment", target());

}
// marginal lik for the polynomial method so we can select power p
generated quantities {
  real marg_lik = multi_normal_cholesky_lpdf(y | mu, L);
  // real theta = square(theta_half);
  // real sigma = square(sigma_half);
}

