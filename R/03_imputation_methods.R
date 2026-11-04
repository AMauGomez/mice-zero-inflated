# Zero-inflated imputation methods for mice ------------------------------------
#
# The thesis imputed the outcome with the `countimp` package (Kleinke & Reinecke),
# methods "zip" and "zinb". countimp is not on CRAN; install it with
#   remotes::install_github("kkleinke/countimp")
# If it is installed, this demo uses countimp's own methods. If not, it falls
# back to a short stand-in that follows the same steps.
#
# Both read the predictor matrix with countimp's convention, as in the thesis:
#
#   predictor matrix code 1 = covariate in both parts of the model
#                         2 = covariate in the count part only
#                         3 = covariate in the zero-inflation part only
#
# so the count process and the zero process can have different covariates
# (in the thesis, chosen from correlograms, PCA, stepwise AIC and 10-fold CV).
#
# Each imputation:
#   1. fits a zero-inflated Poisson or negative binomial model (pscl::zeroinfl)
#      to the observed cases, with separate count and zero-part formulas;
#   2. draws the coefficients from their approximate posterior
#      N(beta_hat, Var(beta_hat)), so parameter uncertainty is propagated
#      (countimp's default, draw = "bayes");
#   3. for each missing case, draws "zero" yes/no and, if no, a count.

if (requireNamespace("countimp", quietly = TRUE)) {
  IMPUTATION_ENGINE <- "countimp"
  mice.impute.zip  <- getFromNamespace("mice.impute.zip",  "countimp")
  mice.impute.zinb <- getFromNamespace("mice.impute.zinb", "countimp")
} else {
  IMPUTATION_ENGINE <- "stand-in (countimp not installed)"

  impute_zi <- function(y, ry, x, wy, type, dist) {
    if (is.null(wy)) wy <- !ry
    if (is.null(type)) type <- rep(1, ncol(x))
    x <- as.data.frame(x)
    count_vars <- colnames(x)[type %in% c(1, 2)]
    zero_vars  <- colnames(x)[type %in% c(1, 3)]
    rhs <- function(v) if (length(v)) paste(v, collapse = " + ") else "1"
    f <- as.formula(paste("y ~", rhs(count_vars), "|", rhs(zero_vars)))

    obs <- data.frame(y = y[ry], x[ry, , drop = FALSE])
    fit <- pscl::zeroinfl(f, data = obs, dist = dist)

    # Draw the coefficients from N(beta_hat, V)
    beta <- MASS::mvrnorm(1, coef(fit), vcov(fit))
    fit$coefficients$count <- beta[startsWith(names(beta), "count_")]
    fit$coefficients$zero  <- beta[startsWith(names(beta), "zero_")]
    names(fit$coefficients$count) <- sub("^count_", "", names(fit$coefficients$count))
    names(fit$coefficients$zero)  <- sub("^zero_",  "", names(fit$coefficients$zero))

    new    <- x[wy, , drop = FALSE]
    p_zero <- predict(fit, newdata = new, type = "zero")
    mu     <- predict(fit, newdata = new, type = "count")
    n <- sum(wy)
    counts <- if (dist == "poisson") rpois(n, mu) else rnbinom(n, size = fit$theta, mu = mu)
    as.integer(ifelse(rbinom(n, 1, p_zero) == 1, 0L, counts))
  }

  mice.impute.zip  <- function(y, ry, x, wy = NULL, type = NULL, ...) impute_zi(y, ry, x, wy, type, "poisson")
  mice.impute.zinb <- function(y, ry, x, wy = NULL, type = NULL, ...) impute_zi(y, ry, x, wy, type, "negbin")
}

cat("Zero-inflated imputation engine:", IMPUTATION_ENGINE, "\n")
