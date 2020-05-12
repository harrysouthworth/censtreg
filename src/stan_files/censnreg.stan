data {
  int<lower=0> N_obs;
  int<lower=0> N_cens;
  int<lower=0> K;
  matrix[N_obs, K] x_obs;
  matrix[N_cens, K] x_cens;
  real y_obs[N_obs];
  real L;

  //real lognu_params[3];
  real sigma_params[2];
}
parameters {
  real<upper=L> y_cens[N_cens];
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ lognormal(sigma_params[1], sigma_params[2]);
  y_obs ~ normal(x_obs * beta, sigma);
  y_cens ~ normal(x_cens * beta, sigma);
}
