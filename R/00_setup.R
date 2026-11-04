# Packages and shared settings ------------------------------------------------
pkgs <- c("mice", "pscl", "MASS", "dplyr", "tidyr", "ggplot2", "readr")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) install.packages(missing_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

SEED <- 2025
N_PATIENTS <- 20000     # demo size (the thesis used 583,678 patients)
SCENARIOS <- c(0.05, 0.10)  # fraction of the outcome set to missing (MCAR)
M <- 5                   # number of imputations
MAXIT <- 5               # MICE iterations

dir.create("results/demo", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
