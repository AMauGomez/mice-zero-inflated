# Masking experiment: PMM vs ZIP vs ZINB -----------------------------------------
#
# Design (same as the thesis):
#   * Start from a dataset where the outcome is known.
#   * Set 5% or 10% of `days_symptoms` to missing at random, independently of
#     everything else (in the thesis: caTools::sample.split).
#   * Covariates with their own missing values are imputed first with PMM, so
#     all three methods start from the same covariates.
#   * Impute the outcome with MICE using PMM, ZIP and ZINB, with a hand-built
#     predictor matrix: code 1 = both model parts, 2 = count part only,
#     3 = zero part only (see R/03_imputation_methods.R).
#   * Compare imputed values with the hidden true values.
#
# Difference from the thesis scripts: the thesis scored the first completed
# dataset (complete(imp)); this demo scores every one of the M imputations and
# averages the metrics.

# 1. Covariates: impute once with PMM (thesis section 5.7.1.1)
cov_meth <- make.method(cohort)
cov_meth[] <- ifelse(colSums(is.na(cohort)) > 0, "pmm", "")
cov_meth["days_symptoms"] <- ""
base_imp <- mice(cohort, m = 1, maxit = MAXIT, method = cov_meth,
                 printFlag = FALSE, seed = SEED)
base <- complete(base_imp)

# 2. Predictor matrix for the outcome (thesis: built from correlograms, PCA,
#    stepwise AIC and cross-validation; here: the variables used to simulate)
pred <- make.predictorMatrix(base)
pred[] <- 0
pred["days_symptoms", c("sex", "fever", "cough", "myalgia", "positive")] <- 2  # count part
pred["days_symptoms", c("intubated", "pneumonia", "sari")]              <- 3  # zero part
pred["days_symptoms", c("age", "dyspnea")]                              <- 1  # both

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
    variance_ratio = var(imputed) / var(truth),
    max_imputed = max(imputed)
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
  for (method in c("pmm", "zip", "zinb")) {
    meth["days_symptoms"] <- method
    t0 <- Sys.time()
    imp <- mice(masked, m = M, maxit = 1, method = meth, predictorMatrix = pred,
                printFlag = FALSE, seed = seed)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    per_imp <- lapply(seq_len(M), function(i) {
      score(truth, complete(imp, i)$days_symptoms[hide])
    })
    res <- colMeans(do.call(rbind, per_imp))
    out[[method]] <- data.frame(missing_pct = 100 * frac, method = toupper(method),
                                t(res), seconds = round(secs, 1))
  }
  do.call(rbind, out)
}

# With only the outcome missing at this stage, one MICE iteration is enough
# (nothing else is updated between iterations).
results <- do.call(rbind, Map(run_scenario, SCENARIOS, SEED + seq_along(SCENARIOS)))
rownames(results) <- NULL
num <- vapply(results, is.numeric, logical(1))
results[num] <- lapply(results[num], round, 3)
readr::write_csv(results, "results/demo/imputation_metrics.csv")
print(results)
