library(censtreg)
library(testthat)

if (FALSE){
  a <- sample(-10:10, size = 1)
  b <- runif(1, -.5, 1.5)
  s <- abs(rnorm(1, 1, 4))
  nu <- runif(1, 3, 7)

  x <- rt(1000, nu) * s
  y <- a +  b * x + rt(1000, nu) * s

  th <- quantile(y, prob = runif(1, .05, .4))

  d <- data.frame(x = x, y = ifelse (y >= th, y, th), realy = y)

  m <- censtreg(y ~ x, data = d, limit = th, silent = TRUE)

  coef(m)
  d$f <- fitted(m)

  f <- fitted(m, method = "summary")

  d <- cbind(d, f)

  ynames <- names(m$model)
  ynames <- ynames[startsWith(ynames, "y_cens")]
  ycens <- rstan::extract(m$model, pars = ynames)
  ycens <- sapply(ycens, what)

  d$ycens <- y
  d$ycens[d$ycens <= th] <- ycens

  d$ycens <- f

  f <- fitted(m, method = "impute")

  ggplot(d, aes(x, y)) +
    geom_point(color = "blue") +
    geom_point(aes(x, f), inherit.aes = FALSE) +
    geom_point(aes(x, realy), color = "orange", inherit.aes = FALSE) +
    geom_point(aes(x, ycens), color = "green", inherit.aes = FALSE) +
    geom_point(aes(x, `1`), color = "red", inherit.aes = FALSE) +
    geom_point(aes(x, `3`), color = "purple", inherit.aes = FALSE) +
    geom_point(aes(x, `2`), color = "steelblue", inherit.aes = FALSE)
}


test_that("Left-censored regression: posterior means match generating mechanism", {

  set.seed(20200330)

  between <- function(x, lower, upper){
    x < upper & x > lower
  }

  m <- list(); i <- 1

  for (i in 1:10){
    a <- sample(-10:10, size = 1)
    b <- runif(1, -.5, 1.5)
    s <- abs(rnorm(1, 1, 4))
    nu <- runif(1, 3, 7)

    x <- rt(1000, nu) * s
    y <- a +  b * x + rt(1000, nu) * s

    th <- quantile(y, prob = runif(1, .05, .4))

    d <- data.frame(x = x, y = ifelse (y >= th, y, th))

    m[[i]] <- censtreg(y ~ x, data = d, limit = th, silent = TRUE)

    sm <- coef(summary(m[[i]]))
    lo <- sm[, 4]
    hi <- sm[, 8]

    truth <- c(a, b, s, nu)

    testthat::expect_true(all(between(truth, lo, hi)),
                          label = "Coefficients look ok: left censored")
  }
})

test_that("Right-censored regression looks ok", {

  set.seed(20200330)

  between <- function(x, lower, upper){
    x < upper & x > lower
  }

  m <- list(); i <- 1

  for (i in 1:10){
    a <- sample(-10:10, size = 1)
    b <- runif(1, -.5, 1.5)
    s <- abs(rnorm(1, 1, 4))
    nu <- runif(1, 3, 7)

    x <- rt(1000, nu) * s
    y <- a +  b * x + rt(1000, nu) * s

    th <- quantile(y, prob = runif(1, .6, .95))

    d <- data.frame(x = x, y = ifelse (y <= th, y, th))

    m[[i]] <- censtreg(y ~ x, data = d, limit = th, silent = TRUE, upper = TRUE)

    s <- coef(summary(m[[i]]))
    lo <- s[, 1] - 2 * s[, 2]
    hi <- s[, 1] + 2 * s[, 2]

    testthat::expect_true(all(between(s[, 1], lo, hi)),
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


