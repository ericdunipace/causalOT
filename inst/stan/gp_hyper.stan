functions {
  
matrix pol_kern( matrix cost,
               int[] z,
                        real theta_0,
                        real theta_1,
                        real gamma_0,
                        real gamma_1,
                        real sigma_0,
                        real sigma_1,
                        real p
                        ) { //kernel
  int N = rows(cost);
  
  matrix[N, N] K;
  
  for (i in 1:N) {
    for (j in 1:N) {
      if ( (z[i] == 1) && ( z[j] == 1) ) {
        K[i,j] = gamma_1 * pow(1.0 + theta_1 * cost[i,j], p);
      } else if ( (z[i] == 0) && (z[j] == 0) ) {
        K[i,j] = gamma_0 * pow(1.0 + theta_0 * cost[i,j], p);
      } else {
        K[i,j] = 0.0;
      }
    }
  }
  for (i in 1:N) {
    if(z[i] == 1) {
       K[i,i] += sigma_1;
    } else if (z[i] == 0) {
       K[i,i] += sigma_0;
    }
   
  }
  return K;
}

matrix rbf_kern( matrix cost,
               int[] z,
                        real theta_0,
                        real theta_1,
                        real gamma_0,
                        real gamma_1,
                        real sigma_0,
                        real sigma_1
                        ) { //kernel
  int N = rows(cost);
  
  matrix[N, N] K;
  
  for (i in 1:N) {
    for (j in 1:N) {
      if ( (z[i] == 1) && ( z[j] == 1) ) {
        K[i,j] = gamma_1 * exp(- 0.5 * theta_1 * cost[i,j]);
      } else if ( (z[i] == 0) && (z[j] == 0) ) {
        K[i,j] = gamma_0 * exp(- 0.5 * theta_0 * cost[i,j]);;
      } else {
        K[i,j] = 0.0;
      }
    }
  }
  for (i in 1:N) {
    if(z[i] == 1) {
       K[i,i] += sigma_1;
    } else if (z[i] == 0) {
       K[i,i] += sigma_0;
    }
   
  }
  return K;
}

  matrix pol_kern_dose(matrix cost_0, matrix cost_1,
                            real theta_0,
                            real theta_1,
                            real gamma_0,
                            real gamma_1,
                            real sigma,
                            real p) { //kernel
    int N = rows(cost_0);

    matrix[N, N] K;

    // real sigma_2 = sigma;

    for(i in 1:N) {
      for(j in 1:N) {
        K[i,j] = gamma_0 * pow(1.0 + theta_0 * cost_0[i,j], p) * gamma_1 * pow(1.0 + theta_1 * cost_1[i,j], p);
      }
    }

    return K;
  }
}

data {
  int<lower=0> N;
  vector[N] y;
  matrix[N,N] discrep;
  matrix[N,N] discrep_z;
  int z[N];
  real p;
  int is_dose;
  int kernel;
}

transformed data {
  vector[N] mu;
  
  for(n in 1:N) mu[n] = 0.0;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> sigma_0;
  real<lower=0> sigma_1[1-is_dose]; //only defined in binary tx eg
  real<lower=0> theta_0;
  real<lower=0> theta_1;
  real<lower=0> gamma_0;
  real<lower=0> gamma_1;
}

transformed parameters {
  matrix[N,N] Sigma;
  
  if(is_dose == 1) { //dose kernel with one variance
    Sigma = pol_kern_dose(discrep_z, discrep,
                            theta_0,
                            theta_1,
                            gamma_0,
                            gamma_1,
                            sigma_0,
                            p);
  } else { //non-dose with two variances
    if(kernel == 1 ) {
      Sigma = pol_kern(discrep, z,
                            theta_0,
                            theta_1,
                            gamma_0,
                            gamma_1,
                            sigma_0,
                            sigma_1[1],
                            p);
    } else if (kernel == 2) {
      Sigma = rbf_kern(discrep, z,
                            theta_0,
                            theta_1,
                            gamma_0,
                            gamma_1,
                            sigma_0,
                            sigma_1[1]);
    }
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and covariance matrix 'Sigma'.
model {
  y ~ multi_normal(mu, Sigma);
}

generated quantities {
  real marg_lik = multi_normal_lpdf(y | mu, Sigma);
}

