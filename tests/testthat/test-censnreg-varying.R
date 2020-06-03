library(censtreg)
library(testthat)

test_that("Left-censored regression, varying limits: posterior means match generating mechanism", {
  set.seed(20200514)

  between <- function(x, lower, upper){
    x < upper & x > lower
  }

  m <- o <- list(); i <- 1

  for (i in 1:10){
    a <- sample(-10:10, size = 1)
    b <- runif(1, -.5, 1.5)
    s <- abs(rnorm(1, 1, 4))
    nu <- runif(1, 20, 100)

    x <- rt(1000, nu) * s
    y <- a +  b * x + rt(1000, nu) * s

    th <- quantile(y, prob = runif(length(y), .05, .4))

    d <- data.frame(x = x, y = ifelse (y >= th, y, -Inf), th = th,
                    censored = as.numeric(y < th))


    m[[i]] <- censtreg(y ~ x, data = d, censored = "censored", limit = d$th,
                       silent = TRUE, nu = Inf)

    o[[i]] <- censtreg(y ~ x, data = d, censored = "censored", limit = d$th,
                       silent = TRUE, nu = Inf, method = "integrate")

    sm <- coef(summary(m[[i]], probs = c(.01, .99)))
    lo <- sm[, 4]
    hi <- sm[, 5]

    truth <- c(a, b, s)

    testthat::expect_true(all(between(truth, lo, hi)),
                          label = "Coefficients look ok: left censored")

    sm <- coef(summary(o[[i]], probs = c(.01, .99)))
    lo <- sm[, 4]
    hi <- sm[, 5]

    testthat::expect_true(all(between(truth, lo, hi)),
                          label = "Coefficients look ok: left censored")
  }
})







test_that("Simulated censored values look sane", {
  a <- sample(-10:10, size = 1)
  b <- runif(1, -.5, 1.5)
  s <- abs(rnorm(1, 1, 4))
  nu <- runif(1, 20, 100)

  x <- rt(1000, nu) * s
  y <- a +  b * x + rt(1000, nu) * s

  th <- quantile(y, prob = runif(1000, .6, .95))

  d <- data.frame(x = x, y = y, th = unname(th), censored = as.numeric(y > th))
  rownames(d) <- 1:nrow(d)

  m <- censtreg(y ~ x, data = d, censored = "censored", limit = d$th,
                silent = TRUE, upper = TRUE, nu = Inf)

  f <- fitted(m, method = "impute", m = 1000)

  wh <- apply(f, 2, min)
  plot(wh, d$th[d$censored == 1])

  expect_true(all(f > unique(th)))
  expect_true(min(f - th) < .001)

  # Need to add test of unimodality - I've read that there's an example in the
  # book by Efron and Tibshirani

  #hist(f[, 1])

})
