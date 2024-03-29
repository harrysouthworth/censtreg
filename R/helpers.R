#' Get the compiled Stan model given nu
getCensModel_estimate <- function(nu){
  if (is.null(nu)){
    message("stan_censtreg")
    stanmod <- stan_censtreg
  } else if (is.infinite(nu)){
    message("stan_censnreg")
    stanmod <- stan_censnreg
  } else { # nu is a specified constant
    if (nu < 1){
      stop("nu < 1 not permitted")
    }
    message("stan_censtreg_fixed_nu")
    stanmod <- stan_censtreg_fixed_nu
  }
}

getCensModel_integrate <- function(nu){
  if (is.null(nu)){
    message("stat_censtreg_int")
    stanmod <- stan_censtreg_int
  } else if (is.infinite(nu)){
    message("stat_censnreg_int")
    stanmod <- stan_censnreg_int
  } else { # nu is a specified constant
    if (nu < 1){
      stop("nu < 1 not permitted")
    }
    message("stan_censtreg_fixed_nu_int")
    stanmod <- stan_censtreg_fixed_nu_int
  }
}

getCensModel <- function(nu, method){
  if (method == "estimate"){
    getCensModel_estimate(nu)
  } else if (method == "integrate"){
    getCensModel_integrate(nu)
  } else {
    stop("method argument in call to censtreg should be either 'estimate' or 'integrate'")
  }
}



checkForTransformations <- function(formula, thecall){
  # Try to check if y is transformed
  yname <- deparse(as.list(formula)[[2]])
  ybracket <- grepl("\\(", yname)
  ypower <- grepl("\\^", yname)

  lname <- deparse(thecall$limit)
  lbracket <- grepl("\\(", lname)
  lpower <- grepl("\\^", lname)

  if (ybracket & !lbracket | ypower & !lpower){
    warning("response appears to be transformed but not the lower limit")
  }

  invisible()
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

  list(data = data.frame(x, y), pars = c(a=a, b=b, s=s, nu=nu))
}

checkCenstreg <- function(model, data){
  pars <- data$pars

  s <- t(summary(model))

  low <- s[4, ] > pars
  high <- s[8, ] < pars

  sum(c(low, high))
}

checkPriorParams <- function(lognu_params, sigma_params){
  if (length(lognu_params) != 3 | !is.numeric(lognu_params)){
    stop("lognu_params should be a numeric vector of length 3")
  }

  if (length(sigma_params) != 2 | !is.numeric(sigma_params)){
    stop("sigma_params should be a numeric vector of length 2")
  }

  invisible()
}
