#' Get the compiled Stan model given nu and upper
getCensModel <- function(nu, upper){
  if (is.null(nu)){
    if (upper) stanmod <- stan_censtreg_u
    else stanmod <- stan_censtreg
  } else if (is.infinite(nu)){
    if (upper) stanmod <- stan_censnreg_u
    else stanmod <- stan_censnreg
  } else { # nu is a specified constant
    if (nu < 1){
      stop("nu < 1 not permitted")
    }
    if (upper0) stanmod <- stan_censtreg_fixed_nu_u
    else stanmod <- stan_censtreg_fixed_nu
  }
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
