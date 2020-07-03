#' @export
print.censtreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  print(x$call)
  cat("\n")

  b <- x$breaches

  cat(paste("There were", nrow(b), "values that breached the limit.\n\n"))

  o <- rstan::summary(x$model)$summary
  o <- o[substring(rownames(o), 1, 7) != "y_cens[", c(1, 3)]
  o <- o[substring(rownames(o), 1, 7) != "y_tmp[", c(1, 3)]
  o <- o[!(rownames(o) %in% c("lognu", "lp__")), ]

  rownames(o)[substring(rownames(o), 1, 5) == "beta["] <- x$names

  print(t(o), digits = digits)
  invisible()
}

#' @method summary censtreg
#' @export
summary.censtreg <- function(x, probs = c(.025, .25, .5, .75, .975),...){
  o <- rstan::summary(x$model, probs = probs)$summary
  o <- o[substring(rownames(o), 1, 7) != "y_cens[", ]
  o <- o[substring(rownames(o), 1, 6) != "y_tmp[", ]
  o <- o[!(rownames(o) %in% c("lognu", "lp__")), ]
  rownames(o)[substring(rownames(o), 1, 5) == "beta["] <- x$names

  o <- list(coefficients = o, breaches = x$breaches)

  class(o) <- "summary.censtreg"
  o
}

#' @method print summary.censtreg
#' @export
print.summary.censtreg <- function(object, digits = max(3L, getOption("digits") - 3L), ...){
  breaches <- object[[2]]

  cat("Values breaching limit:\n")
  print(head(breaches, 6))

  if (nrow(breaches) > 6){
    cat("... ", nrow(breaches) - 6, " other breaches.\n")
  }

  cat("\n")

  object <- t(object[[1]])

  print(object, digits = digits)
  invisible()
}

#' @method plot censtreg
#' @export
plot.censtreg <- function(x, y, what = "trace"){
  m <- x$model
  s <- rstan::summary(m)$summary

  what <- c(rownames(s)[startsWith(rownames(s), "beta[")], "sigma")
  if ("lognu" %in% rownames(s)){
    what <- c(what, "lognu")
  }

  rstan::traceplot(m, pars = what)
}

#' Fitted values from a censtreg object
#' @param object An object of class 'censtreg'.
#' @param method One of "lp" for the linear predictor, "summary" for a summary
#'   of the simulated censored values, or "impute" for multiple imputation.
#' @param what A function to summarize the coefficients and simulated censored
#'   values. Defaults to \code{what = mean} and the mean is used.
#' @param m The number of imputations of the censored response values.
#' @return If the method is "lp" or "summary", a vector; otherwise a matrix.
#' @details It's not very clear to me what use the returned values are when
#'   \code{method = "summary"}.
#' @method fitted censtreg
#' @export
fitted.censtreg <- function(object, method = "lp", what = mean, m = 10, ...){
  co <- coef(object)

  X <- model.matrix(object$formula, object$data)
  y <- model.response(model.frame(object$formula, object$data))

  if (method == "lp") {
    f <- c(X %*% co)
  } else {
    ynames <- names(object$model)
    ynames <- ynames[startsWith(ynames, "y_cens")]
    ycens <- rstan::extract(object$model, pars = ynames)

    if (method == "summary"){
      ycens <- sapply(ycens, what)
      f <- c(X %*% co)
      if (object$upper){
        f[y >= object$limit] <- ycens
      } else {
        f[y <= object$limit] <- ycens
      }
    } else if (method == "impute"){
      ycens <- t(sapply(ycens, function(X){
        sample(X, size = m, replace = FALSE)
      }))

      f <- t(ycens)

      if (FALSE){
        ycens <- matrix(ycens, ncol = m, byrow = FALSE)

        cy <- ifelse(object$upper & y >= object$limit, NA,
                     ifelse(!object$upper & y<= object$limit, NA, y))

        f <- matrix(rep(cy, m), ncol = m, byrow = FALSE)

        f[is.na(f)] <- ycens
      }
    }
  }

  f
}

#' @method coef censtreg
#' @export
coef.censtreg <- function(object, what = mean, ...){
  conames <- names(object$model)
  conames <- conames[startsWith(conames, "beta")]

  nms <- colnames(model.matrix(object$formula, object$data))

  co <- rstan::extract(object$model, conames)
  names(co) <- nms
  sapply(co, what)
}

#' Predicted values from a censtreg object
#' @param object An object of class 'censtreg'.
#' @param newdata A data frame with the predictors. If not specified, the
#'   function passes its arguments through to \code{fitted} unchanged.
#' @param se.fit Whether to include posterior sample deviations of the predictions
#'   in the output. Defaults to \code{se.fit = FALSE}.
#' @param ci.fit Whether to include confidence (credible) intervals in the output. Defaults
#'   to \code{ci.fit = TRUE} and percentiles of the posterior predictive
#'   distribution are included.
#' @param level The confidence (credible) level when \code{ci.fit = TRUE}.
#' @param method One of "lp" for the linear predictor, "summary" for a summary
#'   of the simulated censored values, or "impute" for multiple imputation.
#'   Only used if \code{newdata} is unspecified.
#' @param what A function to summarize the coefficients and simulated censored
#'   values. Defaults to \code{what = mean} and the mean is used.
#'   Only used if \code{newdata} is unspecified.
#' @param m The number of imputations of the censored response values.
#'   Only used if \code{newdata} is unspecified.
#' @method predict censtreg
#' @export
predict.censtreg <- function(object, newdata, se.fit = FALSE, ci.fit = TRUE,
                             level = .950, method = "lp", what = mean, m = 10, ...){
  if (missing(newdata)){
    fitted(object, method = method, what = what, m = m, ...)
  } else{
    nms <- names(object$model)
    betas <- nms[startsWith(nms, "beta")]
    beta <- rstan::extract(object$model, pars = betas)
    beta <- matrix(unlist(beta), ncol = length(beta), byrow = FALSE)

    co <- coef(object)

    X <- model.matrix(object$formula[-2], newdata)

    fit <- c(X %*% co)
    out <- matrix(fit, ncol = 1)
    colnames(out) <- "fit"

    p <- X %*% t(beta)

    if (se.fit){
      se <- apply(p, 1, sd)
      out <- cbind(out, se)
    }

    if (ci.fit){
      ci <- t(apply(p, 1, quantile, probs = c((1 - level) / 2, (1 + level) / 2)))
      out <- cbind(out, ci)
    }

    out
  }
}

#' @export
residuals.censtreg <- function(object, what = mean, ...){
  f <- fitted(object, method = "lp", what = what)
  y <- model.response(model.frame(object$formula, object$data))

  y - f
}


#' Compute WAIC from a censtreg object
#' @param object A censtreg model, as computed by \code{censtreg}.
#' @details The code has been copied and pasted from the Stan GitHub issues
#'   pages.
#' @export
waic.censtreg <- function(object){
  stanfit <- object$model
  log_lik <- extract (stanfit, "log_lik")$log_lik
  lppd <- sum (log (colMeans(exp(log_lik))))
  p_waic_1 <- 2*sum (log(colMeans(exp(log_lik))) - colMeans(log_lik))
  p_waic_2 <- sum (colVars(log_lik))
  waic_2 <- -2*lppd + 2*p_waic_2
  return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1))
}

ggplot.censtreg <- function(data, mapping = NULL, ..., environment = NULL){
  pd <- data.frame(residuals = resid(data), fitted = fitted(data))
  y <- model.response(model.frame(data$formula, data$data))
  pd$y <- y

  pd$limit <- data$limit
  pd$censored <- if (data$upper) { pd$y >= pd$limit } else { pd$y <= pd$limit }

  p1 <- ggplot(pd, aes(fitted, residuals, color = censored)) +
    geom_hline(yintercept = 0) +
    geom_point() +
    ylab("Residuals") + xlab("Fitted values")


  p2 <- ggplot(pd, aes(y, fitted, color = censored)) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point() +
    ylab("Fitted values") + xlab("Observed values")

  gridExtra::grid.arrange(p1, p2, ncol = 2)
}
