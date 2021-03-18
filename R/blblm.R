#' @import purrr
#' @import stats
#' @import furrr
#' @importFrom utils capture.output
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))


#' Little bag of bootstraps on linear regression model
#'
#' Preform the bag of little of bootstraps method to compute the bootstrapped linear regression model.
#'
#' @param formula A linear regression formula to be fitted to the bootstrapped data.
#' @param data A data frame containing the variables in the model.
#' @param m An integer specifying the number of subsets to split the data into.
#' @param B An integer specifying the number of times for the bootstrap procedure to be repeated.
#' @param parallel A boolean value indicating whether to use parallelization to compute fitted blblm model.
#' Default is set to FALSE. If parallelization is desired user must use future::plan() before use of blblm.
#'
#' @return Bag of little bootstraps linear regression model
#' @export
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)
blblm <- function(formula, data, m = 10, B = 5000, parallel = FALSE) {
  data_list <- split_data(data, m)
  if (parallel == TRUE) {
    estimates <- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B),
      .options = furrr_options(seed = TRUE))
    res <- list(estimates = estimates, formula = formula)
  } else {
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
    res <- list(estimates = estimates, formula = formula)
  }
  class(res) <- "blblm"
  invisible(res)
}


#' Split data into m parts of approximated equal sizes.
#'
#' @param data A data frame to be split.
#' @param m An integer specifying the number of subsets to split the data into.
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' Compute the estimates
#'
#' @param formula A linear regression formula to be fitted to the bootstrapped data.
#' @param data A data frame to be used to compute the estimates.
#' @param n An integer
#' @param B An integer specifying the number of repeats.
lm_each_subsample <- function(formula, data, n, B) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  m <- model.frame(formula, data)
  X <- model.matrix(formula, m)
  y <- model.response(m)
  replicate(B, lm1(X, y, n), simplify = FALSE)
}


#' compute the linear regression estimates for a blb dataset
#'
#' @param X matrix
#' @param y vector
#' @param n integer
lm1 <- function(X, y, n) {
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  # lm.wfit can be written in C++ to improve performance
  fit <- lm.wfit(X, y, freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
#'
#' @param fit fitted model
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#'
#' @param fit fitted model
blbsigma <- function(fit) {
  p <- fit$rank
  e <- fit$residuals
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' Print the blblm model
#'
#' Prints blblm fitted model to the console.
#'
#' @param x A fitted blblm model
#' @param ... Additional arguments
#'
#' @return The blblm model
#'
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' Compute sigma of the blblm model
#'
#' Computes the estimated standard deviation of the blblm model.
#'
#' @param object A fitted blblm model.
#' @param confidence A boolean value indicating whether to return a confidince interval for the value of sigma.
#' @param level A numeric value indicating the confidence level. Used when confidence is set to TRUE.
#' @param ... Additional arguments
#'
#' @return A numeric value indicating the value of sigma calculated from the blblm model.
#'
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' Calculate blblm coefficients
#'
#' Calculate the coefficients of the fitted blblm model.
#'
#' @param object A fitted blblm model
#' @param ... Additional arguments
#'
#' @return Coefficients of fitted blblm model.
#'
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' Confidence interval for blblm coefficients
#'
#' Compute confidence interval for the coefficients of a blblm model.
#'
#' @param object A fitted blblm model
#' @param parm The number of parameters in the fitted blblm model. Default is set to NULL.
#' @param level A numeric value indicating the desired level of confidence interval.
#' @param ... Additional parameters
#'
#' @return A confidence interval for the blblm coefficients.
#'
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
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

#' Compute predictions from fitted blblm model
#'
#' Predicts intercept of the linear regression model based on fitted blblm model using new data.
#'
#' @param object A fitted blblm model.
#' @param new_data A data frame containing the the data to be predicted.
#' @param confidence A boolean parameter indicating if a confidence level is desired. Default is set to FALSE.
#' @param level A numeric value indicating the level of the confidence interval if desired.
#' @param ... Additional arguments
#'
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
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


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
