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
    nu <- runif(1, 10, 1e12)

    x <- rt(1000, nu) * s
    y <- a +  b * x + rt(1000, nu) * s

    th <- quantile(y, prob = runif(1, .05, .4))

    d <- data.frame(x = x, y = ifelse (y >= th, y, -Inf))

    m[[i]] <- censtreg(y ~ x, data = d, limit = th, silent = TRUE,
                       nu = Inf)

    o[[i]] <- censtreg(y ~ x, data = d, limit = th, silent = TRUE,
                       method = "integrate", nu = Inf)

    sm <- coef(summary(m[[i]], probs = c(.01, .99)))
    lo <- sm[, 4]
    hi <- sm[, 5]

    truth <- c(a, b, s)

    testthat::expect_true(all(between(truth, lo, hi)),
                          label = "Coefficients look ok: left censored, estimated")

    sm <- coef(summary(o[[i]], probs = c(.01, .99)))
    lo <- sm[, 4]
    hi <- sm[, 5]

    testthat::expect_true(all(between(truth, lo, hi)),
                          label = "Coefficients look ok: left censored, integrated")

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

    d <- data.frame(x = x, y = ifelse (y <= th, y, Inf))

    m[[i]] <- censtreg(y ~ x, data = d, limit = th, silent = TRUE, upper = TRUE,
                       nu = Inf)
    o[[i]] <- censtreg(y ~ x, data = d, limit = th, method = "integrate",
                       silent = TRUE, upper = TRUE, nu = Inf)

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

  expect_error(censtreg(y ~ x, data = d, limit = th),
               label = "No values beneath threshold")

  th <- max(y) * 2

  expect_error(censtreg(y ~ x, data = d, limit = th, upper = TRUE),
               label = "No values above threshold")
})
