#' Linear regression for data with left-censored response and t-distributed residuals
#' @param formula A formula.
#' @param data A data frame. This should contain all the data, including the
#'   censored values. Supposing censoring is from below, setting all response
#'   variable values that are censored to be some value beneath the limit will
#'   result in the correct model. Missing values are not allowed in either the
#'   response or the predictors. The function determines what is censored by the
#'   values of the response and the censoring limit.
#' @param limit A vector that gives the censoring
#'   point, or a single number giving the censoring point. The former allows
#'   differing censoring points over observations.
#'   Note that if the formula contains a transformation of the response, you
#'   must remember to transform the
#'   limit as well (which is why the argument requires a numeric vector, not the
#'   name of a column in the data). The argument \code{upper} dictates if the censoring
#'   limit is a lower limit (i.e. left-censoring, censoring from below) or
#'   an upper limit (right-censoring or censoring from above).
#' @param censored A string naming a column in \code{data} that indicates if
#'   a response value is censored (\code{censored == 1}) or not (\code{censored == 0}).
#' @param upper Logical defaulting to \code{upper = FALSE} indicating that
#'   \code{limit} defines a lower limit, not upper limit, for censoring (i.e.
#'   left-censoring, not right-censoring).
#' @param chains The number of Markov chains to run. Defaults to \code{chains = NULL}
#'   and the function guesses the number based on the number of available cores,
#'   specifically, the number of available cores less 1.
#' @param cores The number of cores to use. Defaults to \code{cores = NULL} and
#'   the function uses as many cores as chains, or the number of available cores
#'   less 1.
#' @param method Either "estimate", in which case the censored values
#'   are treated like parameters and are estimate; or "integrate", in which case
#'   the censored values are integrated out. Defaults to \code{method = "integrate"}.
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
censtreg <- function(formula, data, censored, limit, upper = FALSE, chains=NULL, cores=NULL,
                     method = "integrate",
                     iter = 2000, warmup = 1000,
                     lognu_params = c(nu = 6, mu = log(6), sigma = 1),
                     sigma_params = c(mu = 1, sigma = 100),
                     nu = NULL, silent = FALSE, ...){
  thecall <- match.call()

  if (missing(censored) | missing(limit)){
    stop("censored (string naming column) and limit (vector giving censoring thresholds) must be provided")
  }

  if (length(limit) == 1){
    limit <- rep(limit, length = nrow(data))
  } else if (length(limit) != nrow(data)){
    stop("limit should be a single number or a vector of length nrow(data)")
  }

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

  #checkForTransformations(formula, thecall)
  mf <- model.frame(formula, data)

  y <- model.response(model.frame(formula, data))
  X <- model.matrix(formula, data)

  if (any(is.na(y)) | any(is.na(X))){
    stop("missing values are not allowed")
  }

  xlevels <- lapply(mf, function(X){
    if (is.factor(X) | is.ordered(X)){
      levels(X)
    } else if (is.character(X)){
      unique(X)
    } else {
      NULL
    }
  })
  xlevels[sapply(xlevels, is.null)] <- NULL

  censored <- data[, censored]
  if (!all(sort(unique(censored)) == c(0, 1))){
    stop("censoring variable should have values 0 or 1 only")
  }

  if (colnames(X)[1] == "(Intercept)"){
    rX <- X[, -1, drop = FALSE]
  } else {
    rX <- X
  }

  if (ncol(X) > 1){
    u <- apply(rX, 2, function(z) length(unique(z)))

    if (min(u) == 1){
      stop("One or more covariates has zero variance")
    }
  }

  K <- ncol(X)

  if (upper){
    y <- -y
    limit <- -limit
  }

  i <- censored == 0
  if (any(y[!i] > limit[!i])){
    warning("values labeled as censored but not in breach of limit")
  }

  if (sum(is.na((y[!i]))) > 0 | sum(is.na(X[!i, ])) > 0){
    stop("missing values in uncensored response and associated predictors not allowed")
  }

  if (sum(i) == length(y)){
    stop("there are no uncensored values")
  } else if (sum(i) == 0){
    stop("there are no censored values")
  }

  # Vectors of length 1 need special treatment
  N_cens <- sum(!i)
  if (N_cens == 1){
    L <- array(limit[!i], dim = 1)
  } else {
    L <- limit[!i]
  }

  sdata <- list(y_obs = y[i],
                x_obs = X[i, , drop = FALSE], x_cens = X[!i, , drop = FALSE],
                K = ncol(X), N_obs = sum(i), N_cens = N_cens,
                L = L,
                lognu_params = lognu_params, sigma_params = sigma_params,
                nu = nu)

  # getCensModel reports messages for debugging & testing purposes. We don't want them here
  stanmod <- suppressMessages(getCensModel(nu, method))

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
      which <- (1:length(X))[substring(names(X), 1, 5) %in% c("beta[", "y_cen")]
      X[which] <- lapply(X[which], function(A) -A)
      X
    })
  }

  breaches <- data[censored == 1, ]

  o <- list(model = o, call = thecall, formula = formula, data = data,
            limit = limit, upper = upper, names = colnames(X), breaches = breaches,
            xlevels = xlevels)

  class(o) <- "censtreg"
  o
}
