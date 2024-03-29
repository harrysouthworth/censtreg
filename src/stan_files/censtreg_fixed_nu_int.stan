data {
  int<lower=0> N_obs;
  int<lower=0> N_cens;
  int<lower=0> K;
  matrix[N_obs, K] x_obs;
  matrix[N_cens, K] x_cens;
  real<lower=0> nu;
  real y_obs[N_obs];
  vector[N_cens] L;

  real sigma_params[2];
}
parameters {
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ lognormal(sigma_params[1], sigma_params[2]);
  y_obs ~ student_t(nu, x_obs * beta, sigma);
  target += student_t_lcdf(L | nu, x_cens * beta, sigma);
}
