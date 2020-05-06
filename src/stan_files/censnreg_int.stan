data {
  int<lower=0> N_obs;
  int<lower=0> N_cens;
  int<lower=0> K;
  matrix[N_obs, K] x_obs;

  real y_obs[N_obs];
  real L;

  real sigma_params[2];
}
parameters {
  real<upper=min(y_obs)> L;
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ lognormal(sigma_params[1], sigma_params[2]);
  L ~ normal(mu, sigma);
  y_obs ~ normal(x_obs * beta, sigma);
  target += N_cens * normal_lcdf(L | x_obs * beta, sigma);
}
