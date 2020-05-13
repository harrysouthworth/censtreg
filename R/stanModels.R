## Models for data censored from below, estimating the censored values
stan_censtreg <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censtreg.stan"),
                                   model_name = "censtreg")
stan_censtreg_fixed_nu <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censtreg_fixed_nu.stan"),
                                            model_name = "censtreg_fixed_nu")
stan_censnreg <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censnreg.stan"),
                                   model_name = "censnreg")

################################################################################
## Models for data censored from below, estimating the censored values
stan_censnreg_int <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censnreg_int.stan"),
                                       model_name = "censnreg_int")

stan_censtreg_int <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censtreg_int.stan"),
                                       model_name = "censtreg_int")

stan_censtreg_fixed_nu_int <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censtreg_fixed_nu_int.stan"),
                                                model_name = "censtreg_fixed_nu_int")
