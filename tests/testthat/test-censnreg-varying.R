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
    
    d <- data.frame(x = x, y = ifelse (y >= th, y, -Inf))
    
    m[[i]] <- censtreg(y ~ x, data = d, limit = th, silent = TRUE, nu = Inf)
    
    o[[i]] <- censtreg(y ~ x, data = d, limit = th, silent = TRUE, nu = Inf,
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
