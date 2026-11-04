# Custom MICE imputation methods for zero-inflated counts ----------------------
#
# The thesis used the `countimp` package (methods "zip" and "zinb"). This file
# is a small reimplementation of the same idea so the demo runs on current R:
#
#   1. Bootstrap the observed cases (approximates draws from the posterior
#      of the model parameters, so between-imputation variance is kept).
#   2. Fit a zero-inflated Poisson or negative binomial model on the bootstrap sample.
#   3. For each missing case, draw a structural-zero indicator from the zero
#      part and, if it is not a structural zero, a count from the count part.
#
# mice finds these functions by name: method = "zip_boot" calls mice.impute.zip_boot().

impute_zi <- function(y, ry, x, wy, dist) {
  if (is.null(wy)) wy <- !ry
  obs <- data.frame(y = y[ry], x[ry, , drop = FALSE])
  boot <- obs[sample(nrow(obs), replace = TRUE), , drop = FALSE]
  fit <- pscl::zeroinfl(y ~ ., data = boot, dist = dist)

  new <- data.frame(x[wy, , drop = FALSE])
  p_zero <- predict(fit, newdata = new, type = "zero")
  mu     <- predict(fit, newdata = new, type = "count")
  n <- sum(wy)

  counts <- if (dist == "poisson") rpois(n, mu) else rnbinom(n, size = fit$theta, mu = mu)
  as.integer(ifelse(rbinom(n, 1, p_zero) == 1, 0L, counts))
}

mice.impute.zip_boot  <- function(y, ry, x, wy = NULL, ...) impute_zi(y, ry, x, wy, "poisson")
mice.impute.zinb_boot <- function(y, ry, x, wy = NULL, ...) impute_zi(y, ry, x, wy, "negbin")
