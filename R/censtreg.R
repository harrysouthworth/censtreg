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
#' @param method Either "estimate" the default, in which case the censored values
#'   are treated like parameters and are estimate; or "integrate", in which case
#'   the censored values are integrated out.
#' @param iter The number of steps to run each Markov chain for. Defaults to
#'   \code{iter = 2000}.
#' @param lognu_params Prior parameters for the $t$-distribution used as the prior for
#'   $log nu$. Defaults to \code{lognu_params = c(nu = 6, mu = 6, sigma = 10)}.
#' @param sigma_params Prior parameters for the lognormal distribution used as the
#'   prior for $sigma$. Defaults to \code{sigma_params = c(mu = 1, sigma = 100)}.
#' @param nu A fixed value of kurtosis parameter nu, which might be helpful if there
#'   are convergence issues and fixing this parameter is acceptable. Putting
#'   \code{nu = Inf} results in the residuals being treated as Gaussian. Defaults
#'   to \code{nu = NULL} and nu is treated as a random parameter in the model.
#' @param silent Whether to display Stan's update messages to screen. Defaults
#'   to \code{silent = FALSE}. You'd want to know if the chains were going wrong,
#'   but silencing them can be useful when code is embedded in knitr or
#'   rmarkdown code. If \code{silence = TRUE}, the function tells
#'   \code{rstan::sampling} not to be verbose, not to show warnings, not to
#'   print open progress, and not to refresh. Sheesh!
#' @param ... Additional arguments to \code{rstan::sampling}.
#' @details The function uses the arguments to construct a call to \code{rstan::stan}.
#'   Not many of the available options for \code{rstan::stan} are manipulable
#'   through this function, but such things can be easily added if desired.
#'   If the kurtosis parameter
#'   is treated as random. Computations are done on the scale of the log of the
#'   kurtosis parameter, and a t-distribution centered at log(6) with kurtosis parameter
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
#' @seealso fitted.censtreg
#' @export censtreg
censtreg <- function(formula, data, limit, upper = FALSE, chains=NULL, cores=NULL,
                     method = "estimate",
                     iter = 2000, warmup = 1000,
                     lognu_params = c(nu = 6, mu = log(6), sigma = 1),
                     sigma_params = c(mu = 1, sigma = 100),
                     nu = NULL, silent = FALSE, ...){
  thecall <- match.call()

  checkPriorParams(lognu_params, sigma_params)

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

  if (any(y == limit)){
    stop("values of y exactly at the limit: make them unambiguous by putting them beyond the limit")
  }

  if (sum(is.na(y)) > 0 | sum(is.na(X)) > 0){
    stop("censtreg does not support missing values")
  }

  if (colnames(X)[1] == "(Intercept)"){
    rX <- X[, -1, drop = FALSE]
  } else {
    rX <- X
  }

  u <- apply(rX, 2, function(z) length(unique(z)))
  if (min(u) == 1){
    stop("One or more covariates has zero variance")
  }

  K <- ncol(X)

  if (upper){
    y <- -y
    limit <- -limit
    bl <- data.frame(index = (1:length(y))[y < limit], observed = -y[y < limit])
  } else {
    bl <- data.frame(index = (1:length(y))[y < limit], observed = y[y < limit])
  }


  i <- y > limit

  if (sum(i) == 0){
    stop("there are no uncensored values")
  }

  if (length(limit) == 1){
    limit <- rep(limit, sum(!i))
  }

  # getCensModel reports messages for debugging & testing purposes. We don't want them here
  stanmod <- suppressMessages(getCensModel(nu, method))

  sdata <- list(y_obs = y[i],
                x_obs = X[i, , drop = FALSE], x_cens = X[!i, , drop = FALSE],
                K = ncol(X), N_obs = sum(i), N_cens = sum(!i),
                L = limit,
                lognu_params = lognu_params, sigma_params = sigma_params,
                nu = nu)

  if (silent){
    o <- rstan::sampling(stanmod, data = sdata,
                         cores = cores, chains = chains,
                         iter = iter, warmup = warmup,
                         refresh = 0, show_messages = FALSE, verbose = FALSE,
                         open_progress = FALSE, ...)
  } else {
    o <- rstan::sampling(stanmod, data = sdata,
                         cores = cores, chains = chains,
                         iter = iter, warmup = warmup, ...)
  }

  if (upper){
    limit <- -limit
    # o@samples is a list with one entry for each chain
    # Each chain is a list, not matrix, with an entry for each parameter
    o@sim$samples <- lapply(o@sim$samples, function(X){
      which <- (1:length(X))[substring(names(X), 1, 5) == "beta["]
      X[which] <- lapply(X[which], function(A) -A)

      X
    })
  }

  o <- list(model = o, call = thecall, formula = formula, data = data,
            limit = limit, upper = upper, names = colnames(X), breaches = bl)

  class(o) <- "censtreg"
  o
}
