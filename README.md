# Multiple imputation of zero-inflated counts with MICE

**When the variable to impute has far more zeros than a count distribution allows, can zero-inflated models (ZIP, ZINB) beat Predictive Mean Matching inside MICE?**

BSc thesis in Mathematics, UNAM (2025): *"Métodos de Estimación e Imputación de Variables Cero Infladas con el Modelo de Imputación Multivariada por Ecuaciones Encadenadas"*. Defended 8 August 2025. Supervisor: Dra. Ruth Selene Fuentes García.

📄 [Thesis (TESIUNAM)](https://tesiunamdocumentos.dgb.unam.mx/ptd2025/abr_jun/0872591/Index.html) · 📊 [Full project report](https://amaugomez.github.io/mice-zero-inflated/) · 🗂️ [Original analysis code](thesis_code/) · 🇲🇽 [Resumen en español](https://amaugomez.github.io/mice-zero-inflated/es.html)

---

## At a glance

| | |
|---|---|
| **Data** | Mexico's COVID-19 respiratory-disease surveillance system (SISVER), Jan–Jun 2020: **583,678 patients, 113 variables** |
| **Variable imputed** | Days with symptoms before admission: mean 3.7, variance 10.9, **15.9% zeros** |
| **Imputation** | Multivariate MICE: **42 of 44 variables** had missing values and all were imputed together |
| **Methods compared** | For days with symptoms: PMM, zero-inflated Poisson (ZIP) and zero-inflated negative binomial (ZINB) |
| **Design** | Hide 5% and 10% of the known values, impute them, compare with the truth |
| **Main result** | ZIP had the lowest error (**−29% MSE vs PMM**); ZINB recovered the most zeros (**+38% vs PMM**); both took **5–7× longer** |
| **Second look** | ZIP's lower MSE comes from imputations that are **too narrow**: they keep only ~50% of the real variance. Lowest MSE ≠ best imputation |

## The problem

Days with symptoms is the admission date minus the symptom-onset date. It has more zeros than any standard count model predicts, for two reasons that cannot be told apart in the data:

- **True zeros:** confirmed patients without symptoms.
- **Zeros that are really missing values:** when a patient arrived too ill to answer the intake questionnaire, the onset date was often filled in with the admission date.

Other variables were also incomplete, so the problem calls for **multivariate** multiple imputation. The hypothesis: a model built for excess zeros should impute this variable better than a general-purpose method like PMM.

## How the work was done

The thesis covers the theory first (GLMs and zero-inflated mixture models, missing-data mechanisms and Little's test, EM, FIML, Gibbs sampling and MCMC convergence, Bayesian regression imputation, PMM, MICE as a Gibbs sampler), then this applied study. The [report](https://amaugomez.github.io/mice-zero-inflated/) summarises the theory.

**1. Data preparation** · [`thesis_code/00`](thesis_code/00_data_processing.R), [`01`](thesis_code/01_exploratory_analysis.R)
- 113 → 61 → 44 numeric variables: dropped variables that did not apply to every patient (pregnancy, travel, death dates, hospital coordinates) and converted yes/no items to binary.
- Recoded "SE IGNORA" (unknown) and unusable lab results (rejected, not received, did not amplify…) as missing, so MICE could impute them.
- 1.65 million missing cells before any masking: intubation status 78% missing, anosmia and dysgeusia 42%, antimicrobial treatment 29%, lab result 14%.

![Missing-data pattern](figures/thesis/missingness_pattern.png)

**2. Which distribution?** · [`thesis_code/02`](thesis_code/02_distribution_and_dispersion.R)
- Chi-square goodness-of-fit rejects both Poisson (λ = 3.68) and negative binomial (size = 1.80) for the raw distribution.
- Overdispersion tests and zero-inflation tests on simulated residuals (DHARMa) point to a zero-inflated model.

![Empirical vs Poisson and negative binomial](figures/thesis/empirical_vs_poisson_nb.png)

**3. Model building** · [`thesis_code/03`](thesis_code/03_model_selection.R)
- Candidate predictors from correlograms, PCA and subject-matter judgement.
- Stepwise selection by AIC for a logistic model of "zero days" and for Poisson / NB count models; checks for multicollinearity, residuals and deviance; 10-fold cross-validation.
- The logistic model became the zero part and the count model the count part of ZIP and ZINB. Each part ends up with its own set of predictors.

| Model (listwise-deleted data) | Variables | CV RMSE | Dispersion | Predicted / observed zeros |
|---|---|---|---|---|
| Poisson | 25 | 3.02 | 2.14 (over) | 0.27 |
| Negative binomial | 23 | 3.02 | 0.96 | 0.75 |
| ZIP | 20 count / 24 zero | 3.00 | 1.68 (over) | 0.98 |
| **ZINB** | 20 count / 19 zero | **3.00** | **1.00** | **1.04** |

ZINB was the only model that fixed both the dispersion and the zeros, and had the lowest AIC.

**ZIP and ZINB.** Both mix a logit model for excess zeros (probability π) with a count model: Poisson (ZIP) or negative binomial (ZINB). ZIP: Var(Y) = (1−π)λ(1+πλ). ZINB adds a dispersion parameter α: Var(Y) = (1−π)μ(1+μ(π+α)), which absorbs the overdispersion that ZIP leaves (1.68 → 1.00).

**Why a mixture and not a hurdle model?** In a hurdle model every zero comes from one binary process and the count part cannot produce zeros. Here a zero can be genuine (admitted the day symptoms began) or an artefact (admission date copied as onset date). Only a zero-inflated mixture lets a zero come from either source.

![Final ZINB model](figures/zinb_coefficients.png)

Intubated, immunosuppressed and cardiovascular patients were about twice as likely to have an excess zero; patients reporting symptoms much less. Sudden onset shortened the delay to admission by ~17%.

**4. MICE design** · [`thesis_code/04`](thesis_code/04_masking_and_covariate_imputation.R), [`05a–c`](thesis_code/), [`predictor_matrices/`](thesis_code/predictor_matrices/)
- Multivariate imputation: the 41 other incomplete variables (symptoms, comorbidities, treatments, intubation, lab result) imputed with PMM in the same MICE framework, using a hand-built predictor matrix.
- Days with symptoms imputed with `mice` (PMM) or `countimp` (ZIP, ZINB). The predictor matrices encode the model-building results: **1** = predictor in both parts, **2** = count part only, **3** = zero part only.
- 5% and 10% of the observed values hidden by simple random sampling (about 29,200 and 58,400 values); m = iterations = 5 and 10.
- About 28 hours of imputation in total on a laptop (Core i5, 16 GB RAM).

**5. Convergence** — trace plots of the mean and SD of the imputed values across iterations, for every method and scenario, plus imputed-vs-observed densities:

![ZINB trace plot, 10% missing](figures/thesis/trace_zinb_10.png)

## Results

![Thesis results](figures/thesis_results.png)

| 10% missing (m = 10) | PMM | ZIP | ZINB |
|---|---|---|---|
| MSE | 23.68 | **16.86** | 20.74 |
| MAE | 3.46 | **2.94** | 3.24 |
| Exact matches | 12.3% | **14.2%** | 13.8% |
| Zeros recovered | 22.9% | 29.9% | **31.7%** |
| Imputation time | **1 h 36 min** | 8 h 21 min | 10 h 45 min |

- **ZIP:** lowest MSE, RMSE and MAE in both scenarios; best on non-zero values.
- **ZINB:** best at recovering zeros in both scenarios.
- **PMM:** 4–7× faster, but the least accurate.
- Full table: [`results/thesis/imputation_results_thesis.csv`](results/thesis/imputation_results_thesis.csv)

## A second look: lowest MSE ≠ best imputation

*Post-thesis analysis (2026), using the imputations saved during the thesis.*

ZINB fitted the data best (table above), yet ZIP had the lowest MSE. Going back to the saved imputations explains why:

![MSE vs variance](figures/mse_vs_variance.png)

| Real data, 10% missing | PMM | ZIP | ZINB |
|---|---|---|---|
| Variance of imputed / true values | 1.01 | **0.53** | 0.85 |
| Largest imputed value (true max = 115) | 36 | **18** | 35 |

- An imputation is a **random draw** from a predictive distribution, not a best guess.
- Poisson draws vary less than negative binomial draws, so ZIP's draws sit closer to the mean and score a lower MSE.
- The price: ZIP imputations keep only about **half of the real variance** and never produce long symptom durations. Any analysis on those data would get standard errors that are too small.
- A simulation in which ZINB is the true model reproduces the same pattern (right panel).

**Lesson:** judge an imputation method by whether it preserves the distribution and supports valid inference, not by point-error metrics alone.

## What this means in practice

| Situation | Recommendation |
|---|---|
| First pass, large data, little time | **PMM**: fast, but the least accurate |
| Zeros matter | **ZINB**: most zeros recovered, keeps most of the real spread |
| Overdispersion persists after modelling the zeros | **ZINB** over ZIP |
| Choosing between methods | Compare the variance and range of the imputations, not only the MSE |
| Building the imputation model | Model the zero and count parts separately and pass that to MICE through the predictor matrix |

**Where this applies:** any count with excess zeros and gaps: days in hospital, emergency visits, insurance claims per policy, missed payments, days of work lost in surveys.

## Repository

```
thesis_code/            original analysis scripts (as run; need the non-public data)
  predictor_matrices/   MICE predictor matrices (0/1/2/3 codes)
R/                      runnable demo on simulated data
results/thesis/         aggregated results from the thesis
results/demo/           demo output
figures/                figures (figures/thesis/ = taken from the thesis)
index.qmd → docs/       Quarto report published with GitHub Pages
```

## Run the demo

The SISVER microdata are not public, so the demo runs the same pipeline on a **simulated cohort** with the same structure: binary symptoms and comorbidities, age, severity markers, incomplete covariates, and a zero-inflated, overdispersed outcome.

```r
# R >= 4.3
remotes::install_github("kkleinke/countimp")   # optional, the method used in the thesis
Rscript run_all.R                              # about 2 minutes
```

| Script | What it does |
|---|---|
| `R/01_simulate_data.R` | Simulated cohort (20,000 patients) with the SISVER structure |
| `R/02_count_models.R` | Poisson → NB → ZIP → ZINB: dispersion, predicted zeros, AIC |
| `R/03_imputation_methods.R` | `zip` and `zinb` methods: countimp's own if installed, otherwise a stand-in |
| `R/04_run_experiment.R` | Covariate PMM, masking, PMM vs ZIP vs ZINB, scored over all *m* |
| `R/05_figures.R` | Figures |
| `R/06_thesis_figures.R` | Coefficient plot and predictor-matrix heatmap from thesis artefacts |

**About `R/03`:** the thesis used the `countimp` package (Kleinke & Reinecke), which is on GitHub but not on CRAN. When it is installed, the demo uses countimp's `zip` and `zinb` methods directly. Otherwise it uses a short stand-in that follows the same steps (fit `pscl::zeroinfl` with separate count and zero predictors, draw the coefficients from their approximate posterior, then draw zero / count), so the demo runs anywhere. Demo numbers come from simulated data and are not meant to match the thesis.

## How to cite

See [`CITATION.cff`](CITATION.cff) (GitHub's *Cite this repository* button) or the BibTeX entry in the [report](https://amaugomez.github.io/mice-zero-inflated/#how-to-cite).

## Skills shown

Multiple imputation (MICE, PMM, custom `mice` methods) · count and zero-inflated models (Poisson, NB, ZIP, ZINB) · model selection and diagnostics (stepwise AIC, cross-validation, DHARMa) · MCMC convergence checks · simulation studies · large-data workflow in R (583k patients) · reproducible reports with Quarto

---

Aaron Mauricio Gómez Jiménez · [LinkedIn](https://www.linkedin.com/in/amaugomezjimenez) · maugomez@ciencias.unam.mx
