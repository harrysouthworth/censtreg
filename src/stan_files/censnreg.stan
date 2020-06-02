data {
  int<lower=0> N_obs;
  int<lower=0> N_cens;
  int<lower=0> K;
  matrix[N_obs, K] x_obs;
  matrix[N_cens, K] x_cens;
  real y_obs[N_obs];
  real L[N_cens];
  real sigma_params[2];
}
parameters {
  real y_cens[N_cens];
  for (i in 1:N_cens){
    real<upper=L[i]> tmp;
    tmp = y_cens[i];
    y_cens[i] = tmp;
  }
  //vector[N_cens]<upper=L> y_cens;
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ lognormal(sigma_params[1], sigma_params[2]);
  y_obs ~ normal(x_obs * beta, sigma);
  y_cens ~ normal(x_cens * beta, sigma);
}
