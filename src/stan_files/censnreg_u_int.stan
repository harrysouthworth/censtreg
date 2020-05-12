data {
  int<lower=0> N_obs;
  int<lower=0> N_cens;
  int<lower=0> K;
  matrix[N_obs, K] x_obs;
  matrix[N_cens, K] x_cens;

  real y_obs[N_obs];
  real U;

  real sigma_params[2];
}
parameters {
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ lognormal(sigma_params[1], sigma_params[2]);
  y_obs ~ normal(x_obs * beta, sigma);
  target += normal_lccdf(U | x_cens * beta, sigma);
}
