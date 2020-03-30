#' @export
print.censtreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  print(x$call)
  cat("\n")

  b <- x$breaches

  cat(paste("There were", nrow(b), "values that breached the limit.\n\n"))

  o <- rstan::summary(x$model)$summary
  o <- o[substring(rownames(o), 1, 7) != "y_cens[", c(1, 3)]
  o <- o[!(rownames(o) %in% c("lognu", "lp__")), ]

  rownames(o)[substring(rownames(o), 1, 5) == "beta["] <- x$names

  print(t(o), digits = digits)
  invisible()
}


#' @export
#' @aliases print.summary.censtreg
summary.censtreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  o <- rstan::summary(x$model)$summary
  o <- o[substring(rownames(o), 1, 7) != "y_cens[", ]
  o <- o[!(rownames(o) %in% c("lognu", "lp__")), ]
  rownames(o)[substring(rownames(o), 1, 5) == "beta["] <- x$names

  o <- list(coefficients = o, breaches = x$breaches)

  class(o) <- "summary.censtreg"
  o
}

#' @export
print.summary.censtreg <- function(object, digits = max(3L, getOption("digits") - 3L), ...){
  breaches <- object[[2]]

  cat("Values breaching limit:\n")
  print(head(breaches, 6))

  if (nrow(breaches) > 6){
    cat("... ", nrow(breaches) - 6, " other breaches.\n")
  }

  cat("\n")

  object <- object[[1]]

  print(t(object), digits = digits)
  invisible()
}

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
#' @export
#' @method fitted censtreg
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
      f <- matrix(rep(y, m), ncol = m, byrow = FALSE)
      if (object$upper){
        f[y >= object$limit, ] <- ycens
      } else {
        f[y <= object$limit, ] <- ycens
      }
    }
  }
  
  f
}

#' @export
#' @method coef censtreg
coef.censtreg <- function(object, what = mean, ...){
  conames <- names(object$model)
  conames <- conames[startsWith(conames, "beta")]
  
  nms <- colnames(model.matrix(object$formula, object$data))
  
  co <- rstan::extract(object$model, conames)
  names(co) <- nms
  sapply(co, what)
}