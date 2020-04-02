test_that("The fitted method returns imputed values correctly", {

  set.seed(20200331)

  a <- sample(-10:10, size = 1)
  b <- runif(1, -.5, 1.5)
  s <- abs(rnorm(1, 1, 4))
  nu <- runif(1, 3, 7)

  x <- rt(1000, nu) * s
  y <- a +  b * x + rt(1000, nu) * s

  th <- quantile(y, prob = runif(1, .05, .4))

  d <- data.frame(x = x, y = ifelse (y >= th, y, th))
  d <- d[order(d$y), ]

  censored <- sum(d$y) == th

  m <- censtreg(y ~ x, data = d, limit = th, silent = TRUE)

  n <- sample(5:15, size = 1)
  f <- fitted(m, method = "impute", m = n)

  expect_equal(n, ncol(f), label = "Imputation returns correct number of columns")
  expect_equal(nrow(d), nrow(f), label = "Imputation returns correct number of rows")

  for (j in 1:n){
    expect_true(all(f[censored, j] < d[censored, "y"]),
                label = "Imputed values are less than censoring level")
    expect_true(all(f[!censored, j] == d[censored, "y"]),
                label = "Observed values are the same as source data")
  }
})

test_that("The predict method returns output that makes sense", {
  set.seed(20200401)

  a <- sample(-10:10, size = 1)
  b <- runif(1, -.5, 1.5)
  s <- abs(rnorm(1, 1, 4))
  nu <- runif(1, 3, 7)

  x <- rt(100, nu) * s
  y <- a +  b * x + rt(100, nu) * s

  th <- quantile(y, prob = runif(1, .05, .4))

  d <- data.frame(x = x, y = ifelse (y >= th, y, th))

  censored <- sum(d$y) == th

  m <- censtreg(y ~ x, data = d, limit = th, silent = TRUE)

  p <- predict(m)
  f <- fitted(m)
  expect_true(all(p == f), label = "predict passes to fitted when no newdata")

  p <- predict(m, newdata = d)
  expect_true(all(p[, 1] == f),
              label = "predict returns same point estimates as fitted when appropriate")
  expect_true(all(p[, 2] < p[, 3]),
              label = "upper limit is always greater than lower limit")

  p <- predict(m, newdata = d, ci.fit = FALSE)
  expect_equal(ncol(p), 1, label = "predict returns a matrix with 1 column when se.fit and ci.fit are FALSE")

  p <- predict(m, newdata = d, se.fit = TRUE, ci.fit = FALSE)
  expect_equal(ncol(p), 2, label = "predict adds se when ci.fit is FALSE and se.fit is TRUE")

  p <- predict(m, newdata = d, se.fit = TRUE)
  expect_equal(ncol(p), 4, label = "predict adds se when ci.fit is TRUE and se.fit is TRUE")

  expect_equal(cor(d$x, p[, 1]), 1, label = "predictions have correlation 1 with x")

  p <- p[order(d$x), ]
  lo <- p[1, 4] - p[1, 3]
  mid <- p[nrow(p) / 2, 4] - p[nrow(p) / 2, 3]
  hi <- p[nrow(p), 4] - p[nrow(p), 3]

  expect_lt(mid, hi, label = "CI wider at upper end than middle")
  expect_lt(mid, lo, label = "CI wider at lower end than middle")
})
