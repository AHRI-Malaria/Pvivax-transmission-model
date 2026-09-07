###############################################################################

DESCRIPTIVE_PSEUDOCOUNT <- 0.1

# Final microscopy-based joint Bayesian latent-variable analysis 

###############################################################################

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/")))
  }
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      return(dirname(normalizePath(frames[[i]]$ofile, winslash = "/")))
    }
  }
  normalizePath(getwd(), winslash = "/")
}

required_packages <- c(
  "cmdstanr", "posterior", "dplyr", "tidyr", "purrr", "readr",
  "stringr", "tibble", "ggplot2", "bayesplot", "loo", "scales",
  "patchwork"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install these R packages first: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(bayesplot)
  library(loo)
  library(scales)
  library(patchwork)
})

ROOT <- get_script_dir()
source(file.path(ROOT, "R", "helpers.R"))

if (is.na(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  stop("CmdStan is not installed. Run cmdstanr::install_cmdstan() once.")
}

# --------------------------- Analysis settings -------------------------------
DATA_FILE <- file.path(ROOT, "data", "data_pvivax.csv")
STAN_FILE <- file.path(ROOT, "stan", "joint_pvivax_model.stan")
OUT_DIR <- file.path(ROOT, "results")

# Age 

AGE_CENTER_YEARS <- 15
AGE_INTERNAL_SCALE_YEARS <- 10
AGE_REPORT_INCREMENT_YEARS <- 15
AGE_REPORT_MULTIPLIER <- AGE_REPORT_INCREMENT_YEARS / AGE_INTERNAL_SCALE_YEARS

# Fit objects are cached only inside this analysis directory.
# Set LOAD_EXISTING_FITS below to TRUE to resume an interrupted run.
for (subdir in c("audit", "fits", "tables", "figures", "diagnostics", "loo")) {
  dir.create(file.path(OUT_DIR, subdir), recursive = TRUE, showWarnings = FALSE)
}

SEED <- 20260710
CHAINS <- 4L
PARALLEL_CHAINS <- min(CHAINS, max(1L, parallel::detectCores(logical = TRUE) - 1L))
ADAPT_DELTA <- 0.99
MAX_TREEDEPTH <- 15L


RUN_MODE <- Sys.getenv("PVIVAX_RUN_MODE", unset = "final")
if (!RUN_MODE %in% c("final", "smoke")) stop("PVIVAX_RUN_MODE must be final or smoke.")
if (RUN_MODE == "smoke") {
  PRIMARY_WARMUP <- 300L
  PRIMARY_SAMPLING <- 300L
  SENS_WARMUP <- 250L
  SENS_SAMPLING <- 250L
} else {
  PRIMARY_WARMUP <- 2000L
  PRIMARY_SAMPLING <- 2000L
  SENS_WARMUP <- 1000L
  SENS_SAMPLING <- 1000L
}

LOAD_EXISTING_FITS <- tolower(Sys.getenv("PVIVAX_LOAD_EXISTING_FITS", unset = "false")) %in% c("1", "true", "yes")
RUN_PRIOR_PREDICTIVE <- TRUE
RUN_REDUCED_MODEL <- TRUE
RUN_SENSITIVITY <- TRUE


# study is used as a transparent proxy and varied by 0.5x and 2x.
USER_LOD_P <- c("Study 1" = NA_real_, "Study 2" = NA_real_, "Study 3" = NA_real_)
USER_LOD_G <- c("Study 1" = NA_real_, "Study 2" = NA_real_, "Study 3" = NA_real_)

# Reference values used only for centring regression terms.
CENTER_P <- log10(3000)
CENTER_G <- log10(200)

# Fixed microscopy measurement-error SDs on the log10 scale. They are not
# estimated simultaneously with latent biological SDs because one observation
# per participant cannot identify both components without external information.
PRIMARY_SIGMA_ME_P <- 0.50
PRIMARY_SIGMA_ME_G <- 0.50

# Broad priors. Ethiopian transmission studies inform plausible density scales
# and support a positive biological relationship, but the primary effect priors
# remain centred at zero so the posterior is not forced to be positive.
PRIMARY_PRIORS <- list(
  prior_mu_P_mean = log10(3000),
  prior_mu_P_sd = 1.50,
  prior_mu_G_mean = log10(200),
  prior_mu_G_sd = 1.50,
  prior_latent_sd_scale = 1.25,
  prior_gamma_GP_mean = 0.00,
  prior_gamma_GP_sd = 0.75,
  prior_age_density_sd = 0.50,
  prior_density_study_sd = 1.25,
  prior_alpha_mean = qlogis(0.25),
  prior_alpha_sd = 1.50,
  prior_beta_density_sd = 1.00,
  prior_beta_age_sd = 0.75,
  prior_beta_fever_sd = 1.00,
  prior_beta_study_sd = 1.25,
  prior_log_phi_mean = log(10),
  prior_log_phi_sd = 1.50
)

set.seed(SEED)
options(mc.cores = PARALLEL_CHAINS)

# ----------------------------- Data preparation ------------------------------
if (!file.exists(DATA_FILE)) stop("Data file not found: ", DATA_FILE)
raw <- readr::read_csv(DATA_FILE, show_col_types = FALSE)
required_columns <- c(
  "id", "study", "par", "gam", "number_dissect", "number_inf", "age", "fever"
)
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

primary_dat <- prepare_analysis_data(raw, studies = c(1L, 2L, 3L))

analysis_audit <- tibble::tibble(
  raw_rows = nrow(raw),
  analysed_unique_participants = nrow(primary_dat),
  studies = paste(levels(primary_dat$study_label), collapse = ", "),
  age_reference_used_for_centering = unique(primary_dat$analysis_age_center)[1],
  primary_sigma_me_P_log10 = PRIMARY_SIGMA_ME_P,
  primary_sigma_me_G_log10 = PRIMARY_SIGMA_ME_G,
  retained_zero_gametocyte_records = sum(primary_dat$gam == 0),
  reason_for_retaining_zeros = paste(
    "Microscopy zeros are treated as censored measurements; excluding them",
    "would condition the analysis on microscopy detection."
  )
)
readr::write_csv(analysis_audit, file.path(OUT_DIR, "audit", "analysis_audit.csv"))

make_descriptive_table(primary_dat)
make_observed_data_figure(primary_dat)
make_observed_density_relationship(primary_dat)

primary_data <- build_stan_data(
  primary_dat,
  priors = PRIMARY_PRIORS,
  sigma_me_P = PRIMARY_SIGMA_ME_P,
  sigma_me_G = PRIMARY_SIGMA_ME_G,
  lod_multiplier = 1,
  include_beta_P = 1L,
  prior_only = 0L
)
readr::write_csv(attr(primary_data, "lodP"), file.path(OUT_DIR, "audit", "asexual_thresholds.csv"))
readr::write_csv(attr(primary_data, "lodG"), file.path(OUT_DIR, "audit", "gametocyte_thresholds.csv"))

prior_table <- tibble::tibble(
  prior_parameter = names(PRIMARY_PRIORS),
  value = as.numeric(unlist(PRIMARY_PRIORS, use.names = FALSE))
)
write_table(prior_table, "TableS_prior_hyperparameters.csv")

# ------------------------------- Model fitting -------------------------------
model <- cmdstanr::cmdstan_model(STAN_FILE, force_recompile = FALSE)

if (RUN_PRIOR_PREDICTIVE) {
  prior_data <- build_stan_data(
    primary_dat,
    priors = PRIMARY_PRIORS,
    sigma_me_P = PRIMARY_SIGMA_ME_P,
    sigma_me_G = PRIMARY_SIGMA_ME_G,
    include_beta_P = 1L,
    prior_only = 1L
  )
  prior_fit <- run_or_load_fit(
    model, prior_data, "prior_predictive",
    warmup = 750L, sampling = 1000L, seed = SEED + 1L
  )
  make_prior_predictive_plot(prior_fit, primary_dat)
}

primary_fit <- run_or_load_fit(
  model, primary_data, "primary_all_three_studies",
  warmup = PRIMARY_WARMUP,
  sampling = PRIMARY_SAMPLING,
  seed = SEED + 2L
)

save_hmc_diagnostics(primary_fit, "primary_all_three_studies")
posterior_parameter_table(
  primary_fit, primary_dat, primary_data, "primary_all_three_studies"
)
make_posterior_associations_figure(primary_fit)
make_latent_observed_diagnostic(
  primary_fit, primary_dat, "primary_all_three_studies"
)
make_relationship_figure(
  primary_fit, primary_dat, "primary_all_three_studies"
)
make_posterior_predictive_outputs(
  primary_fit, primary_dat, "primary_all_three_studies"
)
make_probability_contrasts(
  primary_fit, primary_dat, "primary_all_three_studies"
)
make_manuscript_ready_tables(
  primary_fit, primary_dat, primary_data, "primary_all_three_studies"
)

# ---------------------- Residual asexual-stage comparison --------------------
loo_models <- list(full = compute_loo(primary_fit, "full_primary"))

if (RUN_REDUCED_MODEL) {
  reduced_data <- build_stan_data(
    primary_dat,
    priors = PRIMARY_PRIORS,
    sigma_me_P = PRIMARY_SIGMA_ME_P,
    sigma_me_G = PRIMARY_SIGMA_ME_G,
    lod_multiplier = 1,
    include_beta_P = 0L,
    prior_only = 0L
  )
  reduced_fit <- run_or_load_fit(
    model, reduced_data, "without_residual_asexual_term",
    warmup = PRIMARY_WARMUP,
    sampling = PRIMARY_SAMPLING,
    seed = SEED + 3L
  )
  save_hmc_diagnostics(reduced_fit, "without_residual_asexual_term")
  loo_models$no_residual_asexual <- compute_loo(
    reduced_fit, "without_residual_asexual_term"
  )
  compare_loo(loo_models)
}

# ------------------------------ Sensitivity ----------------------------------
sensitivity_results <- sensitivity_summary(
  primary_fit, "Primary: all studies; ME SD 0.50; zeros censored"
)

if (RUN_SENSITIVITY) {
  scenarios <- list(
    list(
      name = "Lower measurement error: SD 0.25",
      studies = c(1L, 2L, 3L),       sigmaP = 0.25, sigmaG = 0.25, lod_mult = 1,
      priors = PRIMARY_PRIORS
    ),
    list(
      name = "Higher measurement error: SD 0.75",
      studies = c(1L, 2L, 3L),       sigmaP = 0.75, sigmaG = 0.75, lod_mult = 1,
      priors = PRIMARY_PRIORS
    ),
    list(
      name = "Lower microscopy threshold: 0.5x",
      studies = c(1L, 2L, 3L),       sigmaP = 0.50, sigmaG = 0.50, lod_mult = 0.5,
      priors = PRIMARY_PRIORS
    ),
    list(
      name = "Higher microscopy threshold: 2x",
      studies = c(1L, 2L, 3L),       sigmaP = 0.50, sigmaG = 0.50, lod_mult = 2,
      priors = PRIMARY_PRIORS
    ),
    list(
      name = "Literature-positive asexual-to-gametocyte prior",
      studies = c(1L, 2L, 3L),       sigmaP = 0.50, sigmaG = 0.50, lod_mult = 1,
      priors = modifyList(
        PRIMARY_PRIORS,
        list(prior_gamma_GP_mean = 0.30, prior_gamma_GP_sd = 0.50)
      )
    ),
    list(
      name = "Wider infectivity coefficient priors",
      studies = c(1L, 2L, 3L),       sigmaP = 0.50, sigmaG = 0.50, lod_mult = 1,
      priors = modifyList(
        PRIMARY_PRIORS,
        list(prior_beta_density_sd = 1.50, prior_beta_fever_sd = 1.50)
      )
    ),
    list(
      name = "Wider beta-binomial precision prior",
      studies = c(1L, 2L, 3L),       sigmaP = 0.50, sigmaG = 0.50, lod_mult = 1,
      priors = modifyList(PRIMARY_PRIORS, list(prior_log_phi_sd = 2.25))
    ),
    list(
      name = "Studies 1 and 2 only",
      studies = c(1L, 2L),       sigmaP = 0.50, sigmaG = 0.50, lod_mult = 1,
      priors = PRIMARY_PRIORS
    )
  )

  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    sc_dat <- prepare_analysis_data(
      raw,
      studies = sc$studies,
      positive_gametocytes_only = FALSE
    )
    sc_data <- build_stan_data(
      sc_dat,
      priors = sc$priors,
      sigma_me_P = sc$sigmaP,
      sigma_me_G = sc$sigmaG,
      lod_multiplier = sc$lod_mult,
      include_beta_P = 1L,
      prior_only = 0L
    )
    fit_name <- paste0(
      "sensitivity_", sprintf("%02d", i), "_",
      stringr::str_replace_all(stringr::str_to_lower(sc$name), "[^a-z0-9]+", "_")
    )
    sc_fit <- run_or_load_fit(
      model, sc_data, fit_name,
      warmup = SENS_WARMUP,
      sampling = SENS_SAMPLING,
      seed = SEED + 100L + i
    )
    save_hmc_diagnostics(sc_fit, fit_name)
    sensitivity_results <- dplyr::bind_rows(
      sensitivity_results,
      sensitivity_summary(sc_fit, sc$name)
    )
    rm(sc_fit)
    invisible(gc())
  }
}

write_table(sensitivity_results, "TableS_sensitivity_analysis.csv")
make_sensitivity_plot(sensitivity_results)

# This is run only after model comparison and sensitivity analyses are complete.
export_complete_manuscript_fields(
  fit = primary_fit,
  dat = primary_dat,
  stan_data = primary_data,
  fit_name = "primary_all_three_studies",
  sensitivity_data = sensitivity_results,
  loo_full = loo_models$full,
  loo_reduced = loo_models$no_residual_asexual %||% NULL
)


capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
writeLines(
  c(
    paste0("Run mode: ", RUN_MODE),
    paste0("CmdStan version: ", cmdstanr::cmdstan_version()),
    paste0("Seed: ", SEED),
    paste0("Chains: ", CHAINS),
    paste0("Primary warm-up iterations per chain: ", PRIMARY_WARMUP),
    paste0("Primary retained iterations per chain: ", PRIMARY_SAMPLING),
    paste0("adapt_delta: ", ADAPT_DELTA),
    paste0("max_treedepth: ", MAX_TREEDEPTH),
    "Thinning: none",
    "Primary zero handling: left-censored microscopy observations",
    "Primary analysis retains microscopy zeros as left-censored observations",
    paste0("Age reference value: ", AGE_CENTER_YEARS, " years"),
    paste0("Age reporting contrast: ", AGE_REPORT_INCREMENT_YEARS, " years")
  ),
  file.path(OUT_DIR, "analysis_settings.txt")
)

