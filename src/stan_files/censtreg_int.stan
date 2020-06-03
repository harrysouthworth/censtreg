data {
  int<lower=0> N_obs;
  int<lower=0> N_cens;
  int<lower=0> K;
  matrix[N_obs, K] x_obs;
  matrix[N_cens, K] x_cens;
  real y_obs[N_obs];
  real L[N_cens];

  real lognu_params[3];
  real sigma_params[2];
}
parameters {
  vector[K] beta;
  real<lower=0> lognu;
  real<lower=0> sigma;
}
transformed parameters{
  real nu = exp(lognu);
}
model {
  sigma ~ lognormal(sigma_params[1], sigma_params[2]);
  lognu ~ student_t(lognu_params[1], lognu_params[2], lognu_params[3]);
  y_obs ~ student_t(exp(lognu), x_obs * beta, sigma);
  target += student_t_lcdf(L | exp(lognu), x_cens * beta, sigma);
}
