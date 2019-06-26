#' @export
print.censtreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  print(x$call)
  cat("\n")

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

  class(o) <- "summary.censtreg"
  o
}

#' @export
print.summary.censtreg <- function(object, digits = max(3L, getOption("digits") - 3L), ...){
  object <- unclass(object)
  print(t(object), digits = digits)
  invisible()
}

#' @export
plot.censtreg <- function(x, y, what = "trace"){
  m <- x$model
  s <- rstan::summary(m)$summary

  what <- c(rownames(s)[startsWith(rownames(s), "beta[")], "sigma", "lognu")

  rstan::traceplot(m, pars = what)
}
