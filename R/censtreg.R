#' Linear regression for data with left-censored response and t-distributed residuals
#' @param formula A formula.
#' @para data A data frame.
#' @param lower An integer dictating the left censoring threshold. There is no
#'   default and this is a required argument. Note that if the formula contains
#'   a transformation of the response, you must remember to transform the lower
#'   limit as well.
#' @param chains The number of Markov chains to run. Defaults to \code{chains = NULL}
#'   and the function guesses the number based on the number of available cores,
#'   specifically, the number of available cores less 1.
#' @param cores The number of cores to use. Defaults to \code{cores = NULL} and
#'   the function uses as many cores as chains, or the number of available cores
#'   less 1.
#' @iter The number of steps to run each Markov chain for. Defaults to
#'   \code{iter = 2000}.
#' @param lognu_params Prior parameters for the $t$-distribution used as the prior for
#'   $log nu$. Defaults to \code{lognu_params = c(nu = 6, mu = 6, sigma = 10)}.
#' @param sigma_params Prior parameters for the lognormal distribution used as the
#'   prior for $sigma$. Defaults to \code{sigma_params = c(mu = 1, sigma = 10)}.
#' @details The function uses the arguments to construct a call to \code{rstan::stan}.
#'   Not many of the available options for \code{rstan::stan} are manipulable
#'   through this function, but such things can be easily added if desired.
#'   The function assumes the residuals are t-distributed and the kurtosis parameter
#'   is treated as random. Computations are done on the scale of the log of the
#'   kurtosis parameter, and a t-distribution centered at 6 with kurtosis parameter
#'   6 is used as a prior because experience suggests taking logs and giving the
#'   sampler some clues about likely values aids chain convergence.
#'   The Stan code is largely copied from Section 4.1, "Censored Data" of the
#'   Stan User Guide, the differences being the ability to let the location
#'   depend on covariates, the use of the t-distribution instead of the Gaussian,
#'   and the prior on the kurtosis parameter.
#' @note It is in principle straightforward to generatlize the function to deal
#'   with right-censoring as well, and to pass priors in.#'
#' @export censtreg
stan_censtreg <- rstan::stan_model(file = file.path(here::here(), "src/stan_files/censtreg.stan"),
                                   model_name = "censtreg")

censtreg <- function(formula, data, lower, chains=NULL, cores=NULL,
                     iter = 2000, warmup = 1000,
                     lognu_params = c(nu = 6, mu = 6, sigma = 10),
                     sigma_params = c(mu = 1, sigma = 10)){
  thecall <- match.call()

  if (is.null(chains)){
    chains <- parallel::detectCores() - 1
    if (chains == 0){ # Only 1 core
      chains <- 1
    }
  }

  if (is.null(cores)){
    cores <- min(parallel::detectCores() - 1, chains)
  }

  # Try to check if y is transformed
  yname <- deparse(as.list(formula)[[2]])
  ybracket <- grepl("\\(", yname)
  ypower <- grepl("\\^", yname)

  lname <- deparse(thecall$lower)
  lbracket <- grepl("\\(", lname)
  lpower <- grepl("\\^", lname)

  if (ybracket & !lbracket | ypower & !lpower){
    warning("response appears to be transformed but not the lower limit")
  }

  y <- model.response(model.frame(formula, data))
  X <- model.matrix(formula, data)

  K <- ncol(X)

  i <- y > lower

  if (sum(i, na.rm=TRUE) == 0 | sum(i, na.rm=TRUE) == length(na.omit(y))){
    stop("either there are no observations beneath the threshold, or none above it")
  }

  sdata <- list(y_obs = y[i],
                x_obs = X[i, ], x_cens = X[!i, ],
                K = ncol(X), N_obs = sum(i), N_cens = sum(!i),
                L = lower,
                lognu_params = lognu_params, sigma_params = sigma_params)

  #o <- stan(file = file.path(here::here(), "stan/censtreg.stan"), data = sdata,
  #          iter = iter, cores = cores, chains = chains)
  o <- rstan::sampling(stan_censtreg, data = sdata,
                       cores = cores, chains = chains, iter = iter, warmup = warmup)

  o <- list(model = o, call = thecall, formula = formula, data = data, lower = lower, names = colnames(X))

  class(o) <- "censtreg"
  o
}

#' @method print censtreg
print.censtreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  print(x$call)
  cat("\n")

  o <- summary(x$model)$summary
  o <- o[substring(rownames(o), 1, 7) != "y_cens[", c(1, 3)]
  o <- o[!(rownames(o) %in% c("lognu", "lp__")), ]

  rownames(o)[substring(rownames(o), 1, 5) == "beta["] <- x$names

  print(t(o), digits = digits)
  invisible()
}


#' @method summary censtreg
#' @aliases print.summary.censtreg
summary.censtreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  o <- summary(x$model)$summary
  o <- o[substring(rownames(o), 1, 7) != "y_cens[", ]
  o <- o[!(rownames(o) %in% c("lognu", "lp__")), ]
  rownames(o)[substring(rownames(o), 1, 5) == "beta["] <- x$names

  class(o) <- "summary.censtreg"
  o
}

print.summary.censtreg <- function(object, digits = max(3L, getOption("digits") - 3L), ...){
  object <- unclass(object)
  print(t(object), digits = digits)
  invisible()
}

#' @method plot censtreg
plot.censtreg <- function(x, y, what = "trace"){
  m <- x$model
  s <- summary(m)$summary

  what <- c(rownames(s)[startsWith(rownames(s), "beta[")], "sigma", "lognu")

  traceplot(m, pars = what)
}

#' Simulate left-censored data to use with \code{censtreg}.
#' @param n The number of datapoints to simulate. Defaults to \code{n = 1000}.
#' @details A single predictor is simulated, the intercept, slope, scale and
#'   kurtosis parameters being randomly sampled. The returned object is a list
#'   containing the data frame and the parameters used to simulate it.
#' @export simCensData
simCensData <- function(n=1000){
  s <- sample(2:7, size=1)
  nu <- sample(2:10, size = 1)
  a <- sample(1:5, size = 1)
  b <- sample(seq(.5, 1.5, by=.1), size = 1)

  x <- rt(n, df=nu)

  y <- a + b * x + rt(n, df=nu) * s

  list(data.frame(x, y), pars = c(a=a, b=b, s=s, nu=nu))
}

checkCenstreg <- function(model, data){
  pars <- data$pars

  s <- t(summary(model))

  low <- s[4, ] > pars
  high <- s[8, ] < pars

  sum(c(low, high))
}

################################################################################
# I ran this and got 4 cases where the 95% posterior interval excluded a
# parameter. From manual inspection, in each case the true value was very
# close to one of the limits of the interval and the posterior means were
# fairly close to the truth
if (FALSE){
  errors <- vector(length=100)
  pars <- matrix(ncol=4, nrow=100)
  estimates <- vector(mode="list", length = 100)

  for (i in 1:20){
    dd <- simCensData(n = 2000)
    plot(dd[[1]])

    o <- censtreg(y ~ x, data = dd[[1]], chains = 4)
    print(summary(o))
    dd$pars
    plot(o)

    errors[i] <- checkCenstreg(o, dd)
    pars[i, ] <- dd$pars
    estimates[[i]] <- summary(o)
  }

  pairs(jitter(pars), col = ifelse(errors == 0, 1, 3))
}
