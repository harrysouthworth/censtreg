## Models for data censored from below
stan_censtreg <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censtreg.stan"),
                                   model_name = "censtreg")
stan_censtreg_fixed_nu <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censtreg_fixed_nu.stan"),
                                            model_name = "censtreg_fixed_nu")
stan_censnreg <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censnreg.stan"),
                                   model_name = "censnreg")

# Models for data censored from above. Presumably these could be done by negating the response
# and then monkeying about with the returned model object. Given the small amount of effor to
# define these as Stan models, the monkeying about part, and then testing the monkeying about,
# seemed like it migth be more unreliable or more work.
stan_censtreg_u <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censtreg_u.stan"),
                                     model_name = "censtreg_u")
stan_censtreg_fixed_nu_u <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censtreg_fixed_nu_u.stan"),
                                              model_name = "censtreg_fixed_nu_u")
stan_censnreg_u <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censnreg_u.stan"),
                                     model_name = "censnreg_u")
