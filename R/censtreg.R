#' Linear regression for data with left-censored response and t-distributed residuals
#' @param formula A formula.
#' @param data A data frame. This should contain all the data, including the
#'   censored values. Supposing censoring is from below, setting all response
#'   variable values that are censored to be some value beneath the limit will
#'   result in the correct model. Missing values are not treated as being
#'   censored, but as missing. The function determines what is censored by the
#'   values of the repsonse and the censoring limit, not by missing values or
#'   censoring flags.
#' @param limit An integer dictating the censoring threshold. There is no
#'   default and this is a required argument. Note that if the formula contains
#'   a transformation of the response, you must remember to transform the
#'   limit as well. The arguement \code{upper} dictates if the censoring
#'   limit is a lower limit (i.e. left-censoring, censoring from below) or
#'   an upper limit (right-censoring or censoring from above).
#' @param upper Logical defaulting to \code{upper = FALSE} indicating that
#'   \code{limit} defines a lower limit, not upper limit, for censoring (i.e.
#'   left-censoring, not right-censoring).
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
#' @param nu A fixed value of kurtosis parameter nu, which might be helpful if there
#'   are convergence issues and fixing this parameter is acceptable. Putting
#'   \code{nu = Inf} results in the residuals being treated as Gaussian. Defaults
#'   to \code{nu = NULL} and nu is treated as a random parameter in the model.
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
#' @note If the kurtosis parameter is allowed to vary, it is not uncommon for
#'   some of the Markov chains to exhibit weird behaviour. One solution is to
#'   specify a tighter prior. Presumably, if you're letting the kurtosis vary,
#'   you've reason to believe it is somewhere between 1 and 10 at the widest,
#'   and more likely between 4 and 6, so priors that constrain the majority
#'   of the mass to those kinds of intervals should only be weakly informative.
#' @export censtreg
censtreg <- function(formula, data, limit, upper = FALSE, chains=NULL, cores=NULL,
                     iter = 2000, warmup = 1000,
                     lognu_params = c(nu = 6, mu = 6, sigma = 10),
                     sigma_params = c(mu = 1, sigma = 10),
                     nu = NULL){
  thecall <- match.call()

  stanmod <- getCensModel(nu, upper)

  if (is.null(chains)){
    chains <- parallel::detectCores() - 1
    if (chains == 0){ # Only 1 core
      chains <- 1
    }
  }

  if (is.null(cores)){
    cores <- min(parallel::detectCores() - 1, chains)
  }

  checkForTransformations(formula, thecall)

  y <- model.response(model.frame(formula, data))
  X <- model.matrix(formula, data)

  if (!upper){
    bl <- data.frame(index = (1:length(y))[y < limit],
                     observed = y[y < limit])
  } else {
    bl <- data.frame(index = (1:length(y))[y > limit],
                     observed = y[y > limit])
  }

  K <- ncol(X)

  if (upper){
    i <- y < limit
  } else {
    i <- y > limit
  }

  if (sum(i, na.rm=TRUE) == 0 | sum(i, na.rm=TRUE) == length(na.omit(y))){
    stop("either there are no observations beneath the threshold, or none above it")
  }

  sdata <- list(y_obs = y[i],
                x_obs = X[i, ], x_cens = X[!i, ],
                K = ncol(X), N_obs = sum(i), N_cens = sum(!i),
                L = limit,
                lognu_params = lognu_params, sigma_params = sigma_params,
                nu = nu)

  o <- rstan::sampling(stanmod, data = sdata,
                       cores = cores, chains = chains, iter = iter, warmup = warmup)

  o <- list(model = o, call = thecall, formula = formula, data = data,
            limit = limit, upper = upper, names = colnames(X), breaches = bl)

  class(o) <- "censtreg"
  o
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

    o <- censtreg(y ~ x, data = dd[[1]], chains = 4, lower = 0)
    print(summary(o))
    dd$pars
    plot(o)

    errors[i] <- checkCenstreg(o, dd)
    pars[i, ] <- dd$pars
    estimates[[i]] <- summary(o)
  }

  pairs(jitter(pars), col = ifelse(errors == 0, 1, 3))
}
