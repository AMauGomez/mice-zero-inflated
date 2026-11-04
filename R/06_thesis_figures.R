# Figures built from thesis artefacts (no microdata needed) -----------------------
#
#  * zinb_coefficients.png: final ZINB imputation model (thesis Tables 5.9-5.10),
#    count part as rate ratios and zero part as odds ratios.
#  * predictor_matrix.png: the hand-built MICE predictor matrix used for the ZIP
#    run (thesis_code/predictor_matrices/), with countimp's 1/2/3 codes.
library(ggplot2)

zi <- readr::read_csv("results/thesis/zi_models_thesis.csv", show_col_types = FALSE)

coefs <- zi |>
  dplyr::filter(variable != "CONSTANTE") |>
  dplyr::select(label, count = zinb_count, zero = zinb_zero) |>
  tidyr::pivot_longer(c(count, zero), names_to = "part", values_to = "beta") |>
  dplyr::filter(!is.na(beta)) |>
  dplyr::mutate(
    ratio = exp(beta),
    part = factor(part, c("zero", "count"),
                  c("Zero part: odds of an excess zero (OR)",
                    "Count part: expected days (rate ratio)")),
    direction = ifelse(beta > 0, "increases", "decreases")
  )
ord <- coefs |> dplyr::group_by(label) |> dplyr::summarise(k = max(abs(beta)))
coefs$label <- factor(coefs$label, ord$label[order(ord$k)])

p_coef <- ggplot(coefs, aes(ratio, label, colour = direction)) +
  geom_vline(xintercept = 1, colour = "grey55") +
  geom_segment(aes(x = 1, xend = ratio, yend = label), linewidth = 0.8) +
  geom_point(size = 2.6) +
  scale_x_log10() +
  scale_colour_manual(values = c(increases = "#C5221F", decreases = "#1A73E8"), guide = "none") +
  facet_wrap(~part, scales = "free_x") +
  labs(x = "exp(coefficient), log scale", y = NULL,
       title = "What drives days with symptoms? Final ZINB imputation model",
       subtitle = "Thesis Tables 5.9-5.10 | 29,373 complete cases | all effects p < 0.05") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
ggsave("figures/zinb_coefficients.png", p_coef, width = 10, height = 7, dpi = 150)

pm <- as.matrix(read.csv("thesis_code/predictor_matrices/predictor_matrix_zip.csv",
                         row.names = 1, check.names = FALSE))
pm <- pm[rowSums(pm != 0) > 0, , drop = FALSE]
grid <- expand.grid(target = rownames(pm), predictor = colnames(pm), stringsAsFactors = FALSE)
grid$code <- pm[cbind(grid$target, grid$predictor)]
grid$code <- factor(grid$code, 0:3, c("not used", "1: predictor (both parts)",
                                      "2: count part only", "3: zero part only"))
lev <- colnames(pm)
grid$target <- factor(grid$target, rev(intersect(lev, rownames(pm))))
grid$predictor <- factor(grid$predictor, lev)

p_pm <- ggplot(grid, aes(predictor, target, fill = code)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("not used" = "#F1F3F4", "1: predictor (both parts)" = "#5F6368",
                               "2: count part only" = "#1A73E8", "3: zero part only" = "#E37400")) +
  labs(x = "Predictor", y = "Variable being imputed", fill = NULL,
       title = "Hand-built MICE predictor matrix (44 variables)",
       subtitle = "Row DIASSINTOMAS: separate predictors for the count and zero parts") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top", panel.grid = element_blank())
ggsave("figures/predictor_matrix.png", p_pm, width = 10, height = 9, dpi = 150)
