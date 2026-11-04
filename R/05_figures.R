# Figures ------------------------------------------------------------------------
pal <- c(PMM = "#9AA0A6", ZIP = "#1A73E8", ZINB = "#E37400")

# 1. Outcome distribution: the excess of zeros
p1 <- ggplot(cohort, aes(days_symptoms)) +
  geom_bar(fill = "#1A73E8", width = 0.8) +
  coord_cartesian(xlim = c(-0.5, 25)) +
  labs(x = "Days with symptoms before admission", y = "Patients",
       title = "Zero-inflated outcome (simulated data)") +
  theme_minimal(base_size = 12)
ggsave("figures/outcome_distribution.png", p1, width = 7, height = 4, dpi = 150)

# 2. Thesis results (original data): MSE and zeros recovered
thesis <- readr::read_csv("results/thesis/imputation_results_thesis.csv", show_col_types = FALSE)
long <- thesis |>
  dplyr::select(method, missing_pct, MSE, zeros_recovered_pct) |>
  tidyr::pivot_longer(c(MSE, zeros_recovered_pct), names_to = "metric") |>
  dplyr::mutate(metric = factor(metric, c("MSE", "zeros_recovered_pct"),
                                c("MSE (lower is better)", "Zeros recovered, % (higher is better)")),
                method = factor(method, c("PMM", "ZIP", "ZINB")),
                missing = factor(paste0(missing_pct, "% missing"), c("5% missing", "10% missing")))

p2 <- ggplot(long, aes(missing, value, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Thesis results: 583,678 patients, SISVER COVID-19 surveillance data") +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave("figures/thesis_results.png", p2, width = 8, height = 4.2, dpi = 150)

# 3. Demo: point accuracy (MSE) vs distributional fidelity (variance ratio)
demo <- readr::read_csv("results/demo/imputation_metrics.csv", show_col_types = FALSE) |>
  dplyr::mutate(method = factor(method, c("PMM", "ZIP", "ZINB")),
                missing = factor(paste0(missing_pct, "% missing"), c("5% missing", "10% missing")))
p3 <- ggplot(demo, aes(variance_ratio, MSE, colour = method, shape = missing)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_point(size = 4) +
  scale_colour_manual(values = pal) +
  labs(x = "Variance of imputed / variance of true values (1 = spread preserved)",
       y = "MSE", colour = NULL, shape = NULL,
       title = "Lowest MSE is not the same as the best imputation",
       subtitle = "Simulated data where ZINB is the true model") +
  theme_minimal(base_size = 12)
ggsave("figures/mse_vs_variance.png", p3, width = 7, height = 4.5, dpi = 150)
