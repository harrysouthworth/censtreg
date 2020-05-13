library(censtreg)
library(testthat)

test_that("Left-censored regression: posterior means match generating mechanism", {
  set.seed(20200330)

  between <- function(x, lower, upper){
    x < upper & x > lower
  }

  m <- o <- list(); i <- 1

  for (i in 1:10){
    a <- sample(-10:10, size = 1)
    b <- runif(1, -.5, 1.5)
    s <- abs(rnorm(1, 1, 4))
    nu <- runif(1, 3, 7)

    x <- rt(1000, nu) * s
    y <- a +  b * x + rt(1000, nu) * s

    th <- quantile(y, prob = runif(1, .05, .4))

    d <- data.frame(x = x, y = ifelse (y >= th, y, th))

    m[[i]] <- censtreg(y ~ x, data = d, limit = th, silent = TRUE,
                       lognu_params = c(nu = 6, mu = log(6), s = 1))

    o[[i]] <- censtreg(y ~ x, data = d, limit = th, silent = TRUE,
                       lognu_params = c(nu = 6, mu = log(6), s = 1),
                       method = "integrate")

    sm <- coef(summary(m[[i]], probs = c(.01, .99)))
    lo <- sm[, 4]
    hi <- sm[, 5]

    truth <- c(a, b, s, nu)

    testthat::expect_true(all(between(truth, lo, hi)),
                          label = "Coefficients look ok: left censored")

    sm <- coef(summary(o[[i]], probs = c(.01, .99)))
    lo <- sm[, 4]
    hi <- sm[, 5]

    testthat::expect_true(all(between(truth, lo, hi)),
                          label = "Coefficients look ok: left censored")
  }
})

test_that("Right-censored regression looks ok", {
  set.seed(20200330)

  between <- function(x, lower, upper){
    x < upper & x > lower
  }

  m <- o <- list(); i <- 1

  for (i in 1:10){
    message(paste("THIS IS RUN", i, " OF 10"))

    a <- sample(-10:10, size = 1)
    b <- runif(1, -.5, 1.5)
    s <- abs(rnorm(1, 1, 4))
    nu <- runif(1, 4, 7)

    x <- rt(1000, nu) * s
    y <- a +  b * x + rt(1000, nu) * s

    th <- quantile(y, prob = runif(1, .6, .95))

    d <- data.frame(x = x, y = ifelse (y < th, y, 2 * th))

    m[[i]] <- censtreg(y ~ x, data = d, limit = th, silent = TRUE, upper = TRUE,
                       lognu_params = c(nu = 6, mu = log(6), s = 1))
    o[[i]] <- censtreg(y ~ x, data = d, limit = th, silent = TRUE, upper = TRUE,
                       lognu_params = c(nu = 6, mu = log(6), s = 1),
                       method = "integrate")

    sm <- coef(summary(m[[i]], probs = c(.01, .99)))
    truth <- c(a, b, s, nu)

    testthat::expect_true(all(between(truth, sm[, 4], sm[, 5])),
                          label = "Coefficients look ok: right censored")

    sm <- coef(summary(o[[i]], probs = c(.01, .99)))

    testthat::expect_true(all(between(truth, sm[, 4], sm[, 5])),
                          label = "Coefficients look ok: right censored")
  }
})

test_that("Fail if no threshold breaches", {
  x <- rt(100, 4)
  y <- x + rt(100, 4)

  d <- data.frame(x = x, y = y)

  th <- min(y) * 2

  expect_error(censtreg(y ~ x, data = d, limit = th),
               label = "No values beneath threshold")

  th <- max(y) * 2

  expect_error(censtreg(y ~ x, data = d, limit = th, upper = TRUE),
               label = "No values above threshold")
})


