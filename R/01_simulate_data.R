# Simulate a dataset with the same structure as the thesis data ---------------
#
# The original data (SISVER, Mexican Ministry of Health, suspected COVID-19
# cases Jan-Jun 2020) are not public. This script builds a synthetic cohort
# with the same kind of variables and the same key feature: a zero-inflated,
# overdispersed count outcome, `days_symptoms` (admission date minus symptom
# onset date).
#
# In the real data the outcome has two kinds of zeros:
#   * true zeros (e.g. asymptomatic confirmed cases), and
#   * recording artefacts: when a severely ill patient arrived and the onset
#     date could not be obtained, clinicians recorded the admission date as the
#     onset date. These zeros are really missing values in disguise.
# The two cannot be told apart in the data. Here the excess zeros come from a
# zero process whose probability rises with severity, as in the real data.

simulate_cohort <- function(n, seed = SEED) {
  set.seed(seed)
  b <- function(p) rbinom(n, 1, p)

  d <- data.frame(
    sex           = b(0.48),
    age           = pmin(pmax(round(rnorm(n, 45, 17)), 0), 100),
    diabetes      = b(0.16),
    hypertension  = b(0.20),
    obesity       = b(0.17),
    copd          = b(0.03),
    ckd           = b(0.03),
    cardiovascular = b(0.03),
    vaccinated    = b(0.30),
    fever         = b(0.80),
    cough         = b(0.82),
    headache      = b(0.75),
    myalgia       = b(0.60),
    diarrhea      = b(0.20),
    rhinorrhea    = b(0.35),
    anosmia       = b(0.15)
  )

  # Severity-related variables depend on age and comorbidities
  lp_sev <- -3 + 0.04 * d$age + 0.5 * d$diabetes + 0.4 * d$obesity +
    0.6 * d$copd + 0.5 * d$ckd
  d$dyspnea   <- rbinom(n, 1, plogis(lp_sev + 1.2))
  d$pneumonia <- rbinom(n, 1, plogis(lp_sev + 0.8 * d$dyspnea))
  d$intubated <- rbinom(n, 1, plogis(lp_sev - 2 + 1.5 * d$pneumonia))
  d$sari      <- rbinom(n, 1, plogis(lp_sev + 1 + d$dyspnea))  # severe acute respiratory infection
  d$positive  <- rbinom(n, 1, plogis(-0.6 + 0.8 * d$anosmia + 0.5 * d$fever))

  # Zero-inflation part: probability of an excess zero rises with severity
  lp_zero <- -2.2 + 1.3 * d$intubated + 0.8 * d$pneumonia + 0.6 * d$sari +
    0.01 * d$age
  # Count part: negative binomial days with symptoms
  lp_count <- 1.4 + 0.004 * d$age + 0.15 * d$fever + 0.20 * d$cough +
    0.15 * d$myalgia - 0.20 * d$dyspnea + 0.10 * d$positive
  excess_zero <- rbinom(n, 1, plogis(lp_zero))
  counts <- rnbinom(n, size = 1.8, mu = exp(lp_count))
  d$days_symptoms <- ifelse(excess_zero == 1, 0L, counts)

  # Small MCAR missingness in some covariates, as in the real data
  for (v in c("intubated", "anosmia", "positive", "vaccinated")) {
    d[[v]][runif(n) < 0.04] <- NA
  }
  d
}

cohort <- simulate_cohort(N_PATIENTS)
readr::write_csv(cohort, "results/demo/simulated_cohort.csv")

cat(sprintf("Simulated %d patients; %.1f%% zeros in days_symptoms; mean %.2f, variance %.2f\n",
            nrow(cohort), 100 * mean(cohort$days_symptoms == 0),
            mean(cohort$days_symptoms), var(cohort$days_symptoms)))
