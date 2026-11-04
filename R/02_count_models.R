# Which count model fits the outcome? ------------------------------------------
# Same diagnostic sequence as the thesis: Poisson -> negative binomial -> ZIP -> ZINB,
# comparing overdispersion, predicted zeros and information criteria.

cc <- na.omit(cohort)
f_count <- days_symptoms ~ sex + age + fever + cough + myalgia + dyspnea + positive +
  intubated + pneumonia + sari
f_zi <- days_symptoms ~ sex + age + fever + cough + myalgia + dyspnea + positive +
  intubated + pneumonia + sari | intubated + pneumonia + sari + age

fits <- list(
  Poisson = glm(f_count, family = poisson, data = cc),
  NegBin  = MASS::glm.nb(f_count, data = cc),
  ZIP     = pscl::zeroinfl(f_zi, data = cc, dist = "poisson"),
  ZINB    = pscl::zeroinfl(f_zi, data = cc, dist = "negbin")
)

expected_zeros <- function(fit) {
  if (inherits(fit, "zeroinfl")) return(sum(predict(fit, type = "prob")[, 1]))
  mu <- fitted(fit)
  if (inherits(fit, "negbin")) return(sum(dnbinom(0, size = fit$theta, mu = mu)))
  sum(dpois(0, mu))
}
dispersion <- function(fit) {
  sum(residuals(fit, type = "pearson")^2) / (nrow(cc) - length(coef(fit)))
}

model_table <- data.frame(
  model = names(fits),
  AIC = vapply(fits, AIC, numeric(1)),
  dispersion_ratio = vapply(fits, dispersion, numeric(1)),
  zeros_predicted = vapply(fits, expected_zeros, numeric(1)),
  zeros_observed = sum(cc$days_symptoms == 0)
)
model_table$zero_ratio <- model_table$zeros_predicted / model_table$zeros_observed
model_table[, -1] <- round(model_table[, -1], 3)
readr::write_csv(model_table, "results/demo/count_model_comparison.csv")
print(model_table, row.names = FALSE)
