data {
  int<lower=0> N_obs;
  int<lower=0>N_cens;
  int<lower=0> K;
  real<lower=0> nu;
  matrix[N_obs, K] x_obs;
  matrix[N_cens, K] x_cens;
  real y_obs[N_obs];
  vector[N_cens] L;
  real sigma_params[2];
}
transformed data{
  real maxL = max(L);
}
parameters {
  vector<upper=max(L)>[N_cens] y_tmp;
  vector[K] beta;
  real<lower=0> sigma;
}
transformed parameters{
  vector[N_cens] y_cens = y_tmp + L - maxL;
}
model {
  sigma ~ lognormal(sigma_params[1], sigma_params[2]);
  y_obs ~ student_t(nu, x_obs * beta, sigma);
  y_cens ~ student_t(nu, x_cens * beta, sigma);
}
