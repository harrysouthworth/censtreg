library(censtreg)
library(testthat)

test_that("Left-censored regression: posterior means match generating mechanism", {

  set.seed(20200512)

  between <- function(x, lower, upper){
    x < upper & x > lower
  }

  m <- o <- list(); i <- 1

  for (i in 1:10){
    message(i)
    a <- sample(-10:10, size = 1)
    b <- runif(1, -.5, 1.5)
    s <- abs(rnorm(1, 1, 4))
    nu <- runif(1, 20, 1e12)

    x <- rt(1000, nu) * s
    y <- a +  b * x + rt(1000, nu) * s

    d <- data.frame(x = x, y = y)#ifelse (y >= th, y, -Inf))
    d$th <- quantile(y, prob = runif(1, .05, .4))
    d$censored <- as.numeric(d$y <= d$th)

    m[[i]] <- censtreg(y ~ x, data = d, censored = "censored", limit = d$th,
                       silent = TRUE, nu = Inf)

    o[[i]] <- censtreg(y ~ x, data = d, censored = "censored", limit = d$th,
                       silent = TRUE, method = "integrate", nu = Inf)

    sm <- coef(summary(m[[i]], probs = c(.01, .99)))
    lo <- sm[, 4]
    hi <- sm[, 5]

    truth <- c(a, b, s)

    testthat::expect_true(all(between(truth, lo, hi)),
                          label = paste0("Coefficients look ok: left censored, estimated. i = ", i))

    sm <- coef(summary(o[[i]], probs = c(.01, .99)))
    lo <- sm[, 4]
    hi <- sm[, 5]

    testthat::expect_true(all(between(truth, lo, hi)),
                          label = paste0("Coefficients look ok: left censored, integrated. i = ", i))

  }
})

test_that("Right-censored regression looks ok", {
  set.seed(20200513)

  between <- function(x, lower, upper){
    x < upper & x > lower
  }

  m <- o <- list(); i <- 1

  for (i in 1:10){
    message(i)
    a <- sample(-10:10, size = 1)
    b <- runif(1, -.5, 1.5)
    s <- abs(rnorm(1, 1, 4))
    nu <- runif(1, 20, 100)

    x <- rt(1000, nu) * s
    y <- a +  b * x + rt(1000, nu) * s

    th <- quantile(y, prob = runif(1, .6, .95))

    d <- data.frame(x = x, y = y, th = unname(th), censored = as.numeric(y > th))

    m[[i]] <- censtreg(y ~ x, data = d, censored = "censored", limit = d$th,
                       silent = TRUE, upper = TRUE, nu = Inf)
    o[[i]] <- censtreg(y ~ x, data = d, censored = "censored", limit = d$th,
                       method = "integrate", silent = TRUE, upper = TRUE, nu = Inf)

    sm <- coef(summary(m[[i]], probs = c(.01, .99)))
    truth <- c(a, b, s)
    lo <- sm[, 4]
    hi <- sm[, 5]

    testthat::expect_true(all(between(truth, lo, hi)),
                          label = "Coefficients look ok: right censored")

    sm <- coef(summary(o[[i]], probs = c(.01, .99)))
    lo <- sm[, 4]
    hi <- sm[, 5]

    testthat::expect_true(all(between(truth, lo, hi)),
                          label = "Coefficients look ok: right censored")
  }
})

test_that("Fail if no threshold breaches", {
  x <- rt(100, 4)
  y <- x + rt(100, 4)

  d <- data.frame(x = x, y = y)

  th <- min(y) * 2

  expect_error(censtreg(y ~ x, data = d, censored == "censored", limit = d$th),
               label = "No values beneath threshold")

  th <- max(y) * 2

  expect_error(censtreg(y ~ x, data = d, censored = "censored", limit = d$th,
                        upper = TRUE), label = "No values above threshold")
})

test_that("Simulated censored values look sane", {
  a <- sample(-10:10, size = 1)
  b <- runif(1, -.5, 1.5)
  s <- abs(rnorm(1, 1, 4))
  nu <- runif(1, 20, 100)

  x <- rt(1000, nu) * s
  y <- a +  b * x + rt(1000, nu) * s

  th <- quantile(y, prob = runif(1, .6, .95))

  d <- data.frame(x = x, y = y, th = unname(th), censored = as.numeric(y > th))

  m <- censtreg(y ~ x, data = d, censored = "censored", limit = d$th,
                silent = TRUE, upper = TRUE, nu = Inf)

  f <- fitted(m, method = "impute", m = 1000)
  expect_true(all(f > unique(th)), label = "Simulated censored values > limit: upper")
  expect_true(min(f - th) < .001, label = "Minimum simulated censored values close to limit: upper")

  # Need to add test of unimodality - I've read that there's an example in the
  # book by Efron and Tibshirani

  #hist(f[, 1])

  a <- sample(-10:10, size = 1)
  b <- runif(1, -.5, 1.5)
  s <- abs(rnorm(1, 1, 4))
  nu <- runif(1, 20, 100)

  x <- rt(1000, nu) * s
  y <- a +  b * x + rt(1000, nu) * s

  th <- quantile(y, prob = runif(1, .05, .4))

  d <- data.frame(x = x, y = y, th = unname(th), censored = as.numeric(y < th))

  m <- censtreg(y ~ x, data = d, censored = "censored", limit = d$th,
                silent = TRUE, upper = FALSE, nu = Inf)

  f <- fitted(m, method = "impute", m = 1000)
  expect_true(all(f < unique(th)), label = "Simulated censored values are < limit: !upper")
  expect_true(min(f - th) < .001, label = "Simulated censored values close to limit: !upper")
})
