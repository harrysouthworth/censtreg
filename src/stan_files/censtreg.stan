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
transformed data{
  real maxL = max(L);
}
parameters {
  real<upper=max(L)> y_tmp[N_cens];
  vector[K] beta;
  /* log(1) = 0 so constrain tail weight at Cauchy. Work with log(nu) because
     experience shows the distribution can be extremely long-tailed to the right
     otherwise. */
  real<lower=0> lognu;
  real<lower=0> sigma;
}
transformed parameters{
  real nu = exp(lognu);
  vector[N_cens] y_cens = y_tmp + L - maxL;
}
model {
  sigma ~ lognormal(sigma_params[1], sigma_params[2]);
  lognu ~ student_t(lognu_params[1], lognu_params[2], lognu_params[3]);
  y_obs ~ student_t(exp(lognu), x_obs * beta, sigma);
  y_cens ~ student_t(exp(lognu), x_cens * beta, sigma);
}
