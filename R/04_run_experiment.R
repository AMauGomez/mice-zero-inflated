# Masking experiment: PMM vs ZIP vs ZINB -----------------------------------------
#
# Design (same as the thesis):
#   * Start from a dataset where the outcome is fully known.
#   * Set 5% or 10% of `days_symptoms` to missing completely at random (MCAR).
#   * Impute with MICE using PMM, ZIP and ZINB for the outcome
#     (covariates with missing values are imputed with logistic regression).
#   * Compare the imputed values with the true hidden values.
#
# Unlike the thesis scripts, which scored the first completed dataset, this
# version scores every one of the M imputations and averages the metrics.

# Covariates in the real data also had missing values. Impute them once first,
# so every method starts from the same complete covariates and only the outcome
# imputation differs.
base_imp <- mice(cohort, m = 1, maxit = MAXIT, printFlag = FALSE, seed = SEED)
base <- complete(base_imp)

score <- function(truth, imputed) {
  zero <- truth == 0
  data.frame(
    MSE  = mean((truth - imputed)^2),
    RMSE = sqrt(mean((truth - imputed)^2)),
    MAE  = mean(abs(truth - imputed)),
    exact_match_pct = 100 * mean(truth == imputed),
    zeros_recovered_pct = 100 * mean(imputed[zero] == 0),
    nonzero_exact_pct = 100 * mean(truth[!zero] == imputed[!zero]),
    imputed_zero_share_pct = 100 * mean(imputed == 0),
    true_zero_share_pct = 100 * mean(zero),
    # Does the imputation keep the spread of the real distribution? (1 = yes)
    variance_ratio = var(imputed) / var(truth)
  )
}

run_scenario <- function(frac, seed) {
  set.seed(seed)
  masked <- base
  hide <- runif(nrow(base)) < frac
  truth <- base$days_symptoms[hide]
  masked$days_symptoms[hide] <- NA

  meth <- make.method(masked)
  meth[] <- ""
  out <- list()
  for (method in c("pmm", "zip_boot", "zinb_boot")) {
    meth["days_symptoms"] <- method
    t0 <- Sys.time()
    imp <- mice(masked, m = M, maxit = 1, method = meth, printFlag = FALSE, seed = seed)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    per_imp <- lapply(seq_len(M), function(i) {
      score(truth, complete(imp, i)$days_symptoms[hide])
    })
    res <- colMeans(do.call(rbind, per_imp))
    out[[method]] <- data.frame(missing_pct = 100 * frac, method = method,
                                t(res), seconds = round(secs, 1))
  }
  do.call(rbind, out)
}

# With only the outcome missing, a single MICE iteration is enough
# (nothing else is updated between iterations).
results <- do.call(rbind, Map(run_scenario, SCENARIOS, SEED + seq_along(SCENARIOS)))
results$method <- factor(results$method, levels = c("pmm", "zip_boot", "zinb_boot"),
                         labels = c("PMM", "ZIP", "ZINB"))
rownames(results) <- NULL
num <- vapply(results, is.numeric, logical(1))
results[num] <- lapply(results[num], round, 3)
readr::write_csv(results, "results/demo/imputation_metrics.csv")
print(results)
