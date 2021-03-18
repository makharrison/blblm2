#' Little bag of bootstraps on logistic regression model
#'
#' Preform the bag of little bootstraps method to compute the bootstrapped logistic regression model.
#'
#' @param formula A logistic regression formula to be fitted to the bootstrapped data.
#' @param data A data frame containing the variables in the model.
#' @param m An integer specifying the number of subsets to split the data into.
#' @param B An integer specifying the number of times for the bootstrap procedure to be repeated.
#'
#' @return Bag of little bootstrapped logistic regression model
#' @export
#' @examples
#' fit <- blblog(vs ~ wt + disp, data = mtcars, m = 3, B = 100)
blblog <- function(formula, data, m = 10, B = 5000) {
  data_list <- split_data(data, m)
  estimates <- map(
    data_list,
    ~ log_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblog"
  invisible(res)
}

#' Compute the estimates
#'
#' @param formula A logistic regression formula to be fitted to the bootstrapped data.
#' @param data A data frame to be used to compute the estimates.
#' @param n An integer
#' @param B An integer specifying the number of repeats.
log_each_subsample <- function(formula, data, n, B) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  m <- model.frame(formula, data)
  X <- model.matrix(formula, m)
  y <- model.response(m)
  replicate(B, log1(X, y, n), simplify = FALSE)
}

#' Compute the logistic regression estimates for a blb dataset
#'
#' @param X matrix
#' @param y vector
#' @param n integer
log1 <- function(X, y, n) {
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  fit <- suppressWarnings(glm.fit(X, y, freqs, family = binomial()))
  list(coef = blbcoef(fit))
}


#' Print the blblog model
#'
#' Prints blblog fitted model to the console.
#'
#' @param x A fitted blblog model
#' @param ... Additional arguments
#'
#' @return The blblog model
#'
#' @export
#' @method print blblog
print.blblog <- function(x, ...) {
  cat("blblog model:", capture.output(x$formula))
  cat("\n")
}


#' Calculate blblog coefficients
#'
#' Calculate the coefficients of the fitted blblog model.
#'
#' @param object A fitted blblog model
#' @param ... Additional arguments
#'
#' @return Coefficients of fitted blblog model.
#'
#' @export
#' @method coef blblog
coef.blblog <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' Confidence interval for blblog coefficients
#'
#' Compute confidence interval for the coefficients of a blblog model.
#'
#' @param object A fitted blblog model
#' @param parm The number of parameters in the fitted blblog model. Default is set to NULL.
#' @param level A numeric value indicating the desired level of confidence interval.
#' @param ... Additional parameters
#'
#' @return A confidence interval for the blblog coefficients.
#'
#' @export
#' @method confint blblog
confint.blblog <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' Compute predictions from fitted blblog model
#'
#' Compute predictions of coefficients of new data using a fitted blblog model.
#'
#' @param object A fitted blblog model.
#' @param new_data A data frame containing the the data to be predicted.
#' @param confidence A boolean parameter indicating if a confidence level is desired. Default is set to FALSE.
#' @param level A numeric value indicating the level of the confidence interval if desired.
#' @param ... Additional arguments
#'
#' @export
#' @method predict blblog
predict.blblog <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}
