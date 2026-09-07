# Helper functions for the final P. vivax joint Bayesian analysis.
#
# This file contains only the functions used by the GitHub reproduction
# workflow. Figure 1 is intentionally excluded. Final Figures 2–6 are
# reproduced from posterior draws or the supplied microscopy data.

`%||%` <- function(x, y) if (is.null(x)) y else x


# =============================================================================
# General output and formatting helpers
# =============================================================================
manuscript_theme <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.8, "lines")
    )
}

save_plot <- function(plot, stem, width = 7.3, height = 5.4, dpi = 600) {
  dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file.path(OUT_DIR, "figures", paste0(stem, ".tiff")),
    plot = plot, device = "tiff", width = width, height = height,
    units = "in", dpi = dpi, compression = "lzw", bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path(OUT_DIR, "figures", paste0(stem, ".pdf")),
    plot = plot, device = grDevices::cairo_pdf, width = width, height = height,
    units = "in", bg = "white"
  )
  invisible(plot)
}

write_table <- function(x, filename) {
  dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(as.data.frame(x), file.path(OUT_DIR, "tables", filename))
  invisible(x)
}

median_iqr <- function(x, digits = 1) {
  q <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"),
    q[[2]], q[[1]], q[[3]]
  )
}

n_percent <- function(x, digits = 1) {
  n <- sum(x, na.rm = TRUE)
  d <- sum(!is.na(x))
  sprintf(paste0("%d/%d (%.", digits, "f%%)"), n, d, 100 * n / d)
}

make_descriptive_table <- function(dat) {
  one_group <- function(d, label) {
    tibble::tibble(
      study = label,
      participants = nrow(d),
      age_years_median_IQR = median_iqr(d$age, 1),
      fever_n_percent = n_percent(d$fever_bin == 1L, 1),
      asexual_density_per_uL_median_IQR = median_iqr(d$par, 1),
      gametocyte_density_per_uL_median_IQR = median_iqr(d$gam, 1),
      mosquitoes_dissected_per_feed_median_IQR = median_iqr(d$number_dissect, 0),
      infected_mosquitoes = sum(d$number_inf),
      dissected_mosquitoes = sum(d$number_dissect),
      overall_infected_percent = 100 * sum(d$number_inf) / sum(d$number_dissect)
    )
  }

  split_dat <- split(dat, dat$study_label, drop = TRUE)
  by_study <- purrr::imap_dfr(
    split_dat,
    function(d, nm) one_group(d, as.character(nm))
  )
  out <- dplyr::bind_rows(by_study, one_group(dat, "Overall"))
  write_table(out, "Table1_descriptive_microscopy_data.csv")
  out
}


# =============================================================================
# Data preparation and Stan-data construction
# =============================================================================
derive_lod_table <- function(dat, variable, user_values, multiplier = 1) {
  out <- dat |>
    dplyr::group_by(study_id, study_label) |>
    dplyr::summarise(
      minimum_positive = {
        z <- .data[[variable]][.data[[variable]] > 0]
        if (length(z) == 0L) NA_real_ else min(z)
      },
      .groups = "drop"
    ) |>
    dplyr::arrange(study_id)

  supplied <- unname(user_values[as.character(out$study_label)])
  supplied_ok <- is.finite(supplied) & supplied > 0
  out$lod_unscaled <- ifelse(supplied_ok, supplied, out$minimum_positive)
  out$lod <- multiplier * out$lod_unscaled
  out$source <- ifelse(
    supplied_ok,
    "Protocol-specified",
    "Smallest positive study-specific microscopy density"
  )

  if (any(!is.finite(out$lod)) || any(out$lod <= 0)) {
    stop("A positive microscopy threshold could not be derived for every study.")
  }
  out
}

prepare_analysis_data <- function(raw, studies = c(1L, 2L, 3L),
                                  positive_gametocytes_only = FALSE) {
  dat <- raw |>
    dplyr::mutate(
      id = as.character(id),
      study_original = as.integer(study),
      fever_text = stringr::str_to_lower(stringr::str_trim(as.character(fever))),
      fever_bin = dplyr::case_when(
        fever_text %in% c("yes", "y", "1", "true") ~ 1L,
        fever_text %in% c("no", "n", "0", "false") ~ 0L,
        TRUE ~ NA_integer_
      )
    ) |>
    dplyr::filter(study_original %in% studies)

  key_vars <- c(
    "par", "gam", "number_dissect", "number_inf", "age", "fever_bin"
  )

  duplicated_keys <- dat |>
    dplyr::count(study_original, id, name = "n") |>
    dplyr::filter(n > 1L)

  if (nrow(duplicated_keys) > 0L) {
    duplicate_rows <- dat |>
      dplyr::semi_join(duplicated_keys, by = c("study_original", "id"))
    readr::write_csv(
      duplicate_rows,
      file.path(OUT_DIR, "audit", "duplicate_participant_rows.csv")
    )
    disagreement <- duplicate_rows |>
      dplyr::group_by(study_original, id) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(key_vars), ~dplyr::n_distinct(.x, na.rm = FALSE)),
        .groups = "drop"
      )
    if (any(as.matrix(disagreement[key_vars]) > 1L)) {
      stop("Conflicting duplicate study-participant records were found. See audit file.")
    }
    dat <- dat |>
      dplyr::distinct(study_original, id, .keep_all = TRUE)
  }

  missing_counts <- vapply(dat[key_vars], function(x) sum(is.na(x)), numeric(1))
  if (any(missing_counts > 0)) {
    readr::write_csv(
      tibble::tibble(variable = names(missing_counts), n_missing = missing_counts),
      file.path(OUT_DIR, "audit", "missing_core_variables.csv")
    )
    stop("Missing values remain in core model variables. See audit/missing_core_variables.csv.")
  }

  if (any(dat$par < 0 | dat$gam < 0)) stop("Microscopy densities cannot be negative.")
  if (any(dat$number_dissect < 1)) stop("number_dissect must be at least one.")
  if (any(dat$number_inf < 0 | dat$number_inf > dat$number_dissect)) {
    stop("number_inf must be between zero and number_dissect.")
  }

  if (isTRUE(positive_gametocytes_only)) {
    dat <- dat |>
      dplyr::filter(gam > 0)
  }

  age_center <- AGE_CENTER_YEARS
  dat |>
    dplyr::mutate(
      study_label = factor(
        paste0("Study ", study_original),
        levels = paste0("Study ", sort(unique(study_original)))
      ),
      study_id = as.integer(study_label),
      age10 = (age - age_center) / AGE_INTERNAL_SCALE_YEARS,
      observed_proportion = number_inf / number_dissect,
      any_infected = as.integer(number_inf > 0),
      analysis_age_center = age_center,
      analysis_age_scale = AGE_INTERNAL_SCALE_YEARS
    ) |>
    dplyr::arrange(study_id, id)
}

build_stan_data <- function(dat, priors, sigma_me_P, sigma_me_G,
                            lod_multiplier = 1, include_beta_P = 1L,
                            prior_only = 0L,
                            use_y = rep(1L, nrow(dat))) {
  lodP <- derive_lod_table(dat, "par", USER_LOD_P, lod_multiplier)
  lodG <- derive_lod_table(dat, "gam", USER_LOD_G, lod_multiplier)

  ans <- c(
    list(
      N = nrow(dat),
      J = nlevels(dat$study_label),
      m = as.integer(dat$number_dissect),
      y = as.integer(dat$number_inf),
      age10 = as.numeric(dat$age10),
      fever = as.integer(dat$fever_bin),
      study = as.integer(dat$study_id),
      P_positive = as.integer(dat$par > 0),
      G_positive = as.integer(dat$gam > 0),
      logP_obs = ifelse(dat$par > 0, log10(dat$par), 0),
      logG_obs = ifelse(dat$gam > 0, log10(dat$gam), 0),
      log_lod_P = log10(lodP$lod),
      log_lod_G = log10(lodG$lod),
      sigma_me_P = sigma_me_P,
      sigma_me_G = sigma_me_G,
      center_P = CENTER_P,
      center_G = CENTER_G,
      include_beta_P = as.integer(include_beta_P),
      prior_only = as.integer(prior_only),
      use_y = as.integer(use_y)
    ),
    priors
  )
  attr(ans, "lodP") <- lodP
  attr(ans, "lodG") <- lodG
  ans
}


# =============================================================================
# Model fitting and posterior extraction
# =============================================================================
fit_storage_id <- function(fit_name) {
  if (identical(fit_name, "prior_predictive")) return("prior")
  if (identical(fit_name, "primary_all_three_studies")) return("primary")
  if (identical(fit_name, "without_residual_asexual_term")) return("reduced")

  sensitivity_id <- stringr::str_match(
    fit_name,
    "^sensitivity_([0-9]{2})_"
  )[, 2]
  if (!is.na(sensitivity_id)) return(paste0("s", sensitivity_id))

  cleaned <- stringr::str_replace_all(
    stringr::str_to_lower(fit_name),
    "[^a-z0-9]+",
    "_"
  )
  substr(cleaned, 1L, 18L)
}

run_or_load_fit <- function(model, stan_data, fit_name, warmup, sampling, seed) {
  # Keep filesystem names deliberately short. This avoids the classic Windows
  # MAX_PATH problem when a long project directory is combined with long
  # sensitivity-scenario names and CmdStan's default CSV basename.
  storage_id <- fit_storage_id(fit_name)
  fit_dir <- file.path(OUT_DIR, "fits", storage_id)
  rds_path <- file.path(fit_dir, paste0(storage_id, ".rds"))

  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(LOAD_EXISTING_FITS) && file.exists(rds_path)) {
    message("Loading existing fit: ", rds_path)
    return(readRDS(rds_path))
  }

  message("\nFitting: ", fit_name, " [cache id: ", storage_id, "]")

  # Sample in a short temporary directory. The permanent project path can be
  # very long on Windows, while CmdStan creates one CSV per chain and otherwise
  # adds the model name, timestamp, and random characters to each filename.
  stage_dir <- file.path(
    tempdir(),
    paste0("pntd_", storage_id, "_", as.integer(seed))
  )
  if (dir.exists(stage_dir)) unlink(stage_dir, recursive = TRUE, force = TRUE)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

  fit <- model$sample(
    data = stan_data,
    seed = seed,
    chains = CHAINS,
    parallel_chains = PARALLEL_CHAINS,
    iter_warmup = warmup,
    iter_sampling = sampling,
    adapt_delta = ADAPT_DELTA,
    max_treedepth = MAX_TREEDEPTH,
    refresh = 200,
    output_dir = stage_dir,
    output_basename = storage_id,
    save_warmup = FALSE
  )

  # Save a self-contained fit first. CmdStanR's save_object() forces the draws
  # and diagnostics to be read before serialization, so the RDS remains usable
  # even if the raw CSV files are later moved or removed.
  fit$save_object(rds_path)

  # Preserve raw chain CSV files under short deterministic names. If this copy
  # step fails, the self-contained RDS above is still valid.
  old_csv <- Sys.glob(file.path(fit_dir, "draws-*.csv"))
  if (length(old_csv) > 0L) unlink(old_csv, force = TRUE)

  moved <- try(
    fit$save_output_files(
      dir = fit_dir,
      basename = "draws",
      timestamp = FALSE,
      random = FALSE
    ),
    silent = TRUE
  )

  if (inherits(moved, "try-error") || anyNA(moved)) {
    warning(
      "The fitted model was saved successfully as ", rds_path,
      ", but one or more raw CmdStan CSV files could not be moved. ",
      "The RDS is self-contained and can still be used for all analyses."
    )
  } else {
    # Update the saved object's internal CSV paths after moving the files.
    fit$save_object(rds_path)
  }

  fit
}

indexed_draw_matrix <- function(fit, variable, expected_columns = NULL) {
  mat <- as.matrix(fit$draws(variable, format = "matrix"))
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1L)
  cn <- colnames(mat)
  if (ncol(mat) > 1L && !is.null(cn)) {
    pattern <- paste0("^", variable, "\\[(\\d+)\\]$")
    idx <- suppressWarnings(as.integer(sub(pattern, "\\1", cn)))
    if (all(is.finite(idx))) {
      mat <- mat[, order(idx), drop = FALSE]
    }
  }
  if (!is.null(expected_columns) && ncol(mat) != expected_columns) {
    stop(variable, " has ", ncol(mat), " columns; expected ", expected_columns, ".")
  }
  mat
}

scalar_draw <- function(fit, variable) {
  mat <- as.matrix(fit$draws(variable, format = "matrix"))
  if (ncol(mat) != 1L) stop(variable, " is not scalar.")
  as.numeric(mat[, 1])
}

subsample_rows <- function(n, max_draws = 2000L, seed = SEED) {
  if (n <= max_draws) return(seq_len(n))
  set.seed(seed)
  sort(sample.int(n, max_draws))
}

save_hmc_diagnostics <- function(fit, fit_name) {
  diag <- fit$diagnostic_summary()
  readr::write_csv(
    as.data.frame(diag),
    file.path(OUT_DIR, "diagnostics", paste0("HMC_", fit_name, ".csv"))
  )

  key_vars <- c(
    "mu_P", "mu_G", "sigma_P", "sigma_G", "gamma_G_P",
    "gamma_P_age", "gamma_G_age", "alpha", "beta_G", "beta_P",
    "beta_age", "beta_fever", "phi"
  )
  # Exclude parameters fixed to a constant in a sensitivity model (for example,
  # beta_P in the reduced model). Rank plots are not meaningful for constants.
  draw_matrix <- as.matrix(fit$draws(variables = key_vars, format = "matrix"))
  variable_sd <- apply(draw_matrix, 2, stats::sd)
  plot_vars <- names(variable_sd)[is.finite(variable_sd) & variable_sd > 1e-12]
  if (length(plot_vars) == 0L) {
    warning("No non-constant variables available for trace/rank plots: ", fit_name)
    return(invisible(diag))
  }
  draws <- fit$draws(variables = plot_vars)
  trace <- bayesplot::mcmc_trace(draws, facet_args = list(ncol = 3)) +
    ggplot2::labs(x = "Post-warm-up iteration", y = NULL) +
    manuscript_theme(8)
  rank <- bayesplot::mcmc_rank_overlay(draws) + manuscript_theme(8)
  save_plot(trace, paste0("FigS_trace_", fit_name), 9.5, 8)
  save_plot(rank, paste0("FigS_rank_", fit_name), 9.5, 8)
  invisible(diag)
}


# =============================================================================
# Tables, manuscript fields, model comparison, and sensitivity summaries
# =============================================================================
posterior_parameter_table <- function(fit, dat, stan_data, fit_name) {
  J <- nlevels(dat$study_label)
  vars <- c(
    "mu_P", "mu_G", "sigma_P", "sigma_G", "gamma_P_age",
    "gamma_G_age", "gamma_G_P", "alpha", "beta_G", "beta_P",
    "beta_age", "beta_fever", "phi", "assoc_P_via_G",
    "assoc_P_residual", "assoc_P_total", "assoc_age_total"
  )
  if (J > 1L) {
    vars <- c(
      vars,
      paste0("delta_P_study[", 2:J, "]"),
      paste0("delta_G_study[", 2:J, "]"),
      paste0("beta_study[", 2:J, "]")
    )
  }

  sm <- posterior::summarise_draws(
    fit$draws(variables = vars),
    median = stats::median,
    lower_95 = function(x) stats::quantile(x, 0.025),
    upper_95 = function(x) stats::quantile(x, 0.975),
    mean = base::mean,
    sd = stats::sd,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  ) |>
    dplyr::rename(parameter = variable)

  fixed_me <- tibble::tibble(
    parameter = c("sigma_me_P_fixed", "sigma_me_G_fixed"),
    median = c(stan_data$sigma_me_P, stan_data$sigma_me_G),
    lower_95 = c(stan_data$sigma_me_P, stan_data$sigma_me_G),
    upper_95 = c(stan_data$sigma_me_P, stan_data$sigma_me_G),
    mean = c(stan_data$sigma_me_P, stan_data$sigma_me_G),
    sd = c(NA_real_, NA_real_),
    rhat = c(NA_real_, NA_real_),
    ess_bulk = c(NA_real_, NA_real_),
    ess_tail = c(NA_real_, NA_real_)
  )
  sm <- dplyr::bind_rows(sm, fixed_me)
  write_table(sm, paste0("Table2_parameter_estimates_", fit_name, ".csv"))

  or_variables <- c("beta_G", "beta_P", "beta_age", "beta_fever")
  if (J > 1L) or_variables <- c(or_variables, paste0("beta_study[", 2:J, "]"))
  or_sm <- purrr::map_dfr(or_variables, function(v) {
    x <- as.numeric(as.matrix(fit$draws(v, format = "matrix"))[, 1])
    multiplier <- if (identical(v, "beta_age")) AGE_REPORT_MULTIPLIER else 1
    z <- exp(multiplier * x)
    tibble::tibble(
      parameter = v,
      odds_ratio = stats::median(z),
      lower_95 = stats::quantile(z, 0.025),
      upper_95 = stats::quantile(z, 0.975),
      posterior_probability_OR_above_1 = mean(z > 1)
    )
  })
  write_table(or_sm, paste0("Table2_odds_ratios_", fit_name, ".csv"))
  list(parameters = sm, odds_ratios = or_sm)
}

make_manuscript_ready_tables <- function(fit, dat, stan_data, fit_name) {
  get_summary <- function(variable, transform = identity, digits = 2) {
    x <- as.numeric(as.matrix(fit$draws(variable, format = "matrix"))[, 1])
    z <- transform(x)
    c(
      estimate = stats::median(z),
      lower = stats::quantile(z, 0.025),
      upper = stats::quantile(z, 0.975),
      prob_positive = mean(x > 0)
    )
  }
  rows <- list(
    c("Asexual-stage to gametocyte association", "Change in latent log10 gametocyte density per 1-log10 higher latent asexual-stage density", get_summary("gamma_G_P")),
    c("Gametocyte density and mosquito infection", "Odds ratio per 10-fold higher gametocyte density", get_summary("beta_G", exp)),
    c("Residual conditional asexual-stage association", "Odds ratio per 10-fold higher asexual-stage density", get_summary("beta_P", exp)),
    c("Host age and mosquito infection", "Odds ratio per 15-year higher age", get_summary("beta_age", function(x) exp(AGE_REPORT_MULTIPLIER * x))),
    c("Fever and mosquito infection", "Odds ratio: fever versus no fever", get_summary("beta_fever", exp)),
    c("Beta-binomial precision", "Precision parameter; lower values indicate greater overdispersion", get_summary("phi")),
    c("Latent asexual-stage SD", "Log10 scale", get_summary("sigma_P")),
    c("Latent gametocyte SD", "Log10 scale", get_summary("sigma_G"))
  )
  out <- purrr::map_dfr(rows, function(r) {
    vals <- as.numeric(r[3:6])
    tibble::tibble(
      parameter = r[[1]], interpretation = r[[2]],
      posterior_median = vals[1], lower_95_CrI = vals[2], upper_95_CrI = vals[3],
      posterior_probability_positive = vals[4]
    )
  })
  out <- dplyr::bind_rows(
    out,
    tibble::tibble(
      parameter = c("Fixed asexual-stage measurement-error SD", "Fixed gametocyte measurement-error SD"),
      interpretation = "Fixed value on the log10 scale; varied in sensitivity analyses",
      posterior_median = c(stan_data$sigma_me_P, stan_data$sigma_me_G),
      lower_95_CrI = NA_real_, upper_95_CrI = NA_real_, posterior_probability_positive = NA_real_
    )
  )
  write_table(out, "Table2_main_Bayesian_estimates.csv")

  contrasts <- readr::read_csv(
    file.path(OUT_DIR, "tables", paste0("standardized_probability_contrasts_", fit_name, ".csv")),
    show_col_types = FALSE
  )
  write_table(contrasts, "Table3_probability_contrasts.csv")

  manuscript_fields <- tibble::tibble(
    field = c(
      "N_ANALYSIS", "N_INFECTIOUS", "PCT_INFECTIOUS",
      "OR_GAM", "OR_GAM_LO", "OR_GAM_HI",
      "OR_ASEX_RES", "OR_ASEX_RES_LO", "OR_ASEX_RES_HI",
      "PHI", "PHI_LO", "PHI_HI",
      "SIGMA_P", "SIGMA_G"
    ),
    value = c(
      nrow(dat), sum(dat$number_inf > 0), 100 * mean(dat$number_inf > 0),
      get_summary("beta_G", exp)[1:3],
      get_summary("beta_P", exp)[1:3],
      get_summary("phi")[1:3],
      get_summary("sigma_P")[1], get_summary("sigma_G")[1]
    )
  )
  write_table(manuscript_fields, "manuscript_result_fields.csv")
  invisible(list(table2 = out, table3 = contrasts, fields = manuscript_fields))
}

export_complete_manuscript_fields <- function(
    fit, dat, stan_data, fit_name,
    sensitivity_data = NULL,
    loo_full = NULL,
    loo_reduced = NULL) {

  qsum <- function(variable, transform = identity) {
    x <- scalar_draw(fit, variable)
    z <- transform(x)
    c(
      estimate = stats::median(z),
      lower = unname(stats::quantile(z, 0.025)),
      upper = unname(stats::quantile(z, 0.975))
    )
  }
  fmt <- function(x, digits = 2) formatC(as.numeric(x), digits = digits, format = "f")
  fmt3 <- function(x) formatC(as.numeric(x), digits = 3, format = "f")
  pct <- function(x, digits = 1) formatC(100 * as.numeric(x), digits = digits, format = "f")

  gam <- qsum("beta_G", exp)
  asex <- qsum("beta_P", exp)
  age <- qsum("beta_age", function(x) exp(AGE_REPORT_MULTIPLIER * x))
  fever <- qsum("beta_fever", exp)
  gp <- qsum("gamma_G_P")
  sigp <- qsum("sigma_P")
  sigg <- qsum("sigma_G")
  phi <- qsum("phi")

  J <- nlevels(dat$study_label)
  study2 <- if (J >= 2) qsum("beta_study[2]", exp) else rep(NA_real_, 3)
  study3 <- if (J >= 3) qsum("beta_study[3]", exp) else rep(NA_real_, 3)

  contrast_file <- file.path(
    OUT_DIR, "tables",
    paste0("standardized_probability_contrasts_", fit_name, ".csv")
  )
  contrasts <- if (file.exists(contrast_file)) {
    readr::read_csv(contrast_file, show_col_types = FALSE)
  } else {
    make_probability_contrasts(fit, dat, fit_name)
  }
  get_contrast <- function(pattern) {
    z <- contrasts[stringr::str_detect(contrasts$contrast, stringr::fixed(pattern)), , drop = FALSE]
    if (nrow(z) != 1L) return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
    c(estimate = z$risk_difference[[1]], lower = z$RD_lower[[1]], upper = z$RD_upper[[1]])
  }
  rd_g <- get_contrast("Tenfold higher gametocyte density")
  rd_pr <- get_contrast("residual conditional association")
  rd_pt <- get_contrast("total model-implied association")

  contrast_row <- function(pattern) {
    z <- contrasts[stringr::str_detect(contrasts$contrast, stringr::fixed(pattern)), , drop = FALSE]
    if (nrow(z) != 1L) return(NULL)
    z
  }
  cg <- contrast_row("Tenfold higher gametocyte density")
  cpr <- contrast_row("residual conditional association")
  cpt <- contrast_row("total model-implied association")
  cage <- contrast_row("15-year higher host age")
  cfev <- contrast_row("Fever versus no fever")
  prob_ci <- function(z, stem) {
    if (is.null(z)) return("NA")
    paste0(
      pct(z[[paste0(stem, "_median")]], 1), "% (95% CrI ",
      pct(z[[paste0(stem, "_lower")]], 1), "%–",
      pct(z[[paste0(stem, "_upper")]], 1), "%)"
    )
  }
  rd_ci <- function(z) {
    if (is.null(z)) return("NA")
    paste0(
      pct(z$risk_difference, 1), " (95% CrI ",
      pct(z$RD_lower, 1), " to ", pct(z$RD_upper, 1), ")"
    )
  }
  rr_ci <- function(z) {
    if (is.null(z)) return("NA")
    paste0(
      fmt(z$risk_ratio), " (95% CrI ", fmt(z$RR_lower), "–", fmt(z$RR_upper), ")"
    )
  }

  # Rank-normalised convergence summaries.
  sm <- fit$summary()
  max_rhat <- suppressWarnings(max(sm$rhat[is.finite(sm$rhat)], na.rm = TRUE))
  min_bulk <- suppressWarnings(min(sm$ess_bulk[is.finite(sm$ess_bulk)], na.rm = TRUE))
  min_tail <- suppressWarnings(min(sm$ess_tail[is.finite(sm$ess_tail)], na.rm = TRUE))

  ds <- as.data.frame(fit$diagnostic_summary())
  n_div <- 0
  n_td <- 0
  min_ebfmi <- NA_real_
  div_col <- grep("diverg", names(ds), value = TRUE, ignore.case = TRUE)
  td_col <- grep("treedepth", names(ds), value = TRUE, ignore.case = TRUE)
  eb_col <- grep("ebfmi|e_bfmi", names(ds), value = TRUE, ignore.case = TRUE)
  if (length(div_col)) n_div <- sum(ds[[div_col[[1]]]], na.rm = TRUE)
  if (length(td_col)) n_td <- sum(ds[[td_col[[1]]]], na.rm = TRUE)
  if (length(eb_col)) min_ebfmi <- min(ds[[eb_col[[1]]]], na.rm = TRUE)

  # Posterior-predictive coverage and mean calibration summary.
  pred_file <- file.path(
    OUT_DIR, "tables", paste0("posterior_predictions_", fit_name, ".csv")
  )
  if (file.exists(pred_file)) {
    pred <- readr::read_csv(pred_file, show_col_types = FALSE)
    coverage <- mean(
      pred$observed >= pred$predictive_lower &
        pred$observed <= pred$predictive_upper,
      na.rm = TRUE
    )
    mean_obs <- mean(pred$observed, na.rm = TRUE)
    mean_pred <- mean(pred$predicted_median, na.rm = TRUE)
    ppc_summary <- paste0(
      pct(coverage, 1), "% of participant-level observed proportions lay within the 95% ",
      "beta-binomial predictive intervals; the overall observed and posterior-median mean ",
      "infection proportions were ", pct(mean_obs, 1), "% and ", pct(mean_pred, 1), "%"
    )
  } else {
    ppc_summary <- "posterior predictive results are provided in Fig 6 and the accompanying output table"
  }

  phi_interpretation <- if (phi[["estimate"]] < 5) {
    "substantial"
  } else if (phi[["estimate"]] < 20) {
    "moderate"
  } else {
    "limited"
  }

  age_fever_summary <- paste0(
    "an age odds ratio of ", fmt(age[["estimate"]]), " (95% CrI ",
    fmt(age[["lower"]]), "–", fmt(age[["upper"]]), ") per 15 years and a fever odds ratio of ",
    fmt(fever[["estimate"]]), " (95% CrI ", fmt(fever[["lower"]]), "–",
    fmt(fever[["upper"]]), ")"
  )

  sensitivity_summary_text <- "not evaluated"
  if (!is.null(sensitivity_data) && nrow(sensitivity_data) > 0L) {
    main_par <- sensitivity_data |>
      dplyr::filter(parameter %in% c("beta_G", "beta_P", "gamma_G_P", "phi"))
    stable_direction <- main_par |>
      dplyr::filter(parameter %in% c("beta_G", "beta_P")) |>
      dplyr::group_by(parameter) |>
      dplyr::summarise(all_above = all(estimate > 1), all_below = all(estimate < 1), .groups = "drop") |>
      dplyr::summarise(stable = all(all_above | all_below)) |>
      dplyr::pull(stable)
    sensitivity_summary_text <- if (isTRUE(stable_direction)) {
      "qualitatively consistent"
    } else {
      "not fully directionally consistent"
    }
  }

  elpd_diff <- NA_real_
  elpd_se <- NA_real_
  if (!is.null(loo_full) && !is.null(loo_reduced)) {
    pf <- loo_full$pointwise[, "elpd_loo"]
    pr <- loo_reduced$pointwise[, "elpd_loo"]
    dd <- pf - pr
    elpd_diff <- sum(dd)
    elpd_se <- sqrt(length(dd) * stats::var(dd))
  }

  values <- c(
    N_ANALYSIS = as.character(nrow(dat)),
    N_INFECTIOUS = as.character(sum(dat$number_inf > 0)),
    PCT_INFECTIOUS = pct(mean(dat$number_inf > 0), 1),
    OR_GAM = fmt(gam[["estimate"]]), OR_GAM_LO = fmt(gam[["lower"]]), OR_GAM_HI = fmt(gam[["upper"]]),
    OR_ASEX_RES = fmt(asex[["estimate"]]), OR_ASEX_RES_LO = fmt(asex[["lower"]]), OR_ASEX_RES_HI = fmt(asex[["upper"]]),
    OR_AGE = fmt(age[["estimate"]]), OR_AGE_LO = fmt(age[["lower"]]), OR_AGE_HI = fmt(age[["upper"]]),
    OR_FEVER = fmt(fever[["estimate"]]), OR_FEVER_LO = fmt(fever[["lower"]]), OR_FEVER_HI = fmt(fever[["upper"]]),
    OR_STUDY2 = fmt(study2[["estimate"]]), OR_STUDY2_LO = fmt(study2[["lower"]]), OR_STUDY2_HI = fmt(study2[["upper"]]),
    OR_STUDY3 = fmt(study3[["estimate"]]), OR_STUDY3_LO = fmt(study3[["lower"]]), OR_STUDY3_HI = fmt(study3[["upper"]]),
    GAMMA_GP = fmt(gp[["estimate"]]), GAMMA_GP_LO = fmt(gp[["lower"]]), GAMMA_GP_HI = fmt(gp[["upper"]]),
    PHI = fmt(phi[["estimate"]]), PHI_LO = fmt(phi[["lower"]]), PHI_HI = fmt(phi[["upper"]]),
    PHI_INTERPRETATION = phi_interpretation,
    SIGMA_P = fmt(sigp[["estimate"]]), SIGMA_P_LO = fmt(sigp[["lower"]]), SIGMA_P_HI = fmt(sigp[["upper"]]),
    SIGMA_G = fmt(sigg[["estimate"]]), SIGMA_G_LO = fmt(sigg[["lower"]]), SIGMA_G_HI = fmt(sigg[["upper"]]),
    RD_GAM = pct(rd_g[["estimate"]], 1), RD_GAM_LO = pct(rd_g[["lower"]], 1), RD_GAM_HI = pct(rd_g[["upper"]], 1),
    RD_ASEX_RES = pct(rd_pr[["estimate"]], 1), RD_ASEX_RES_LO = pct(rd_pr[["lower"]], 1), RD_ASEX_RES_HI = pct(rd_pr[["upper"]], 1),
    RD_ASEX_TOTAL = pct(rd_pt[["estimate"]], 1), RD_ASEX_TOTAL_LO = pct(rd_pt[["lower"]], 1), RD_ASEX_TOTAL_HI = pct(rd_pt[["upper"]], 1),
    P0_GAM = prob_ci(cg, "p0"), P1_GAM = prob_ci(cg, "p1"), RR_GAM = rr_ci(cg),
    P0_ASEX_RES = prob_ci(cpr, "p0"), P1_ASEX_RES = prob_ci(cpr, "p1"), RR_ASEX_RES = rr_ci(cpr),
    P0_ASEX_TOTAL = prob_ci(cpt, "p0"), P1_ASEX_TOTAL = prob_ci(cpt, "p1"), RR_ASEX_TOTAL = rr_ci(cpt),
    P0_AGE = prob_ci(cage, "p0"), P1_AGE = prob_ci(cage, "p1"), RD_AGE = rd_ci(cage), RR_AGE = rr_ci(cage),
    P0_FEVER = prob_ci(cfev, "p0"), P1_FEVER = prob_ci(cfev, "p1"), RD_FEVER = rd_ci(cfev), RR_FEVER = rr_ci(cfev),
    MAX_RHAT = fmt3(max_rhat), MIN_ESS_BULK = format(round(min_bulk), scientific = FALSE),
    MIN_ESS_TAIL = format(round(min_tail), scientific = FALSE),
    N_DIVERGENCES = as.character(n_div), N_TREEDEPTH = as.character(n_td),
    MIN_EBFMI = ifelse(is.finite(min_ebfmi), fmt3(min_ebfmi), "NA"),
    PPC_SUMMARY = ppc_summary,
    AGE_FEVER_SUMMARY = age_fever_summary,
    SENSITIVITY_SUMMARY = sensitivity_summary_text,
    ELPD_DIFF = ifelse(is.finite(elpd_diff), fmt(elpd_diff), "NA"),
    ELPD_SE = ifelse(is.finite(elpd_se), fmt(elpd_se), "NA")
  )
  out <- tibble::tibble(field = names(values), value = unname(values))
  write_table(out, "manuscript_result_fields_complete.csv")
  invisible(out)
}

make_probability_contrasts <- function(fit, dat, fit_name, max_draws = 2000L) {
  N <- nrow(dat)
  logP <- indexed_draw_matrix(fit, "logP", N)
  logG <- indexed_draw_matrix(fit, "logG", N)
  nd_all <- nrow(logP)
  keep <- subsample_rows(nd_all, max_draws, SEED + 30L)
  logP <- logP[keep, , drop = FALSE]
  logG <- logG[keep, , drop = FALSE]

  alpha <- scalar_draw(fit, "alpha")[keep]
  beta_G <- scalar_draw(fit, "beta_G")[keep]
  beta_P <- scalar_draw(fit, "beta_P")[keep]
  beta_age <- scalar_draw(fit, "beta_age")[keep]
  beta_fever <- scalar_draw(fit, "beta_fever")[keep]
  gamma_P_age <- scalar_draw(fit, "gamma_P_age")[keep]
  gamma_G_age <- scalar_draw(fit, "gamma_G_age")[keep]
  gamma_G_P <- scalar_draw(fit, "gamma_G_P")[keep]

  nd <- length(keep)
  study_effect <- matrix(0, nrow = nd, ncol = N)
  if (nlevels(dat$study_label) > 1L) {
    for (j in 2:nlevels(dat$study_label)) {
      bj <- indexed_draw_matrix(fit, "beta_study", nlevels(dat$study_label))[keep, j]
      cols <- which(dat$study_id == j)
      study_effect[, cols] <- matrix(bj, nrow = nd, ncol = length(cols))
    }
  }

  rep_draw <- function(x) matrix(x, nrow = nd, ncol = N)
  rep_person <- function(x) matrix(x, nrow = nd, ncol = N, byrow = TRUE)

  eta0 <- rep_draw(alpha) + study_effect +
    sweep(logG - CENTER_G, 1, beta_G, "*") +
    sweep(logP - CENTER_P, 1, beta_P, "*") +
    outer(beta_age, dat$age10, "*") +
    outer(beta_fever, dat$fever_bin, "*")

  summarise_contrast <- function(label, eta_base, eta_new) {
    p0 <- rowMeans(plogis(eta_base))
    p1 <- rowMeans(plogis(eta_new))
    rd <- p1 - p0
    rr <- p1 / pmax(p0, 1e-10)
    tibble::tibble(
      contrast = label,
      p0_median = stats::median(p0),
      p0_lower = stats::quantile(p0, 0.025),
      p0_upper = stats::quantile(p0, 0.975),
      p1_median = stats::median(p1),
      p1_lower = stats::quantile(p1, 0.025),
      p1_upper = stats::quantile(p1, 0.975),
      risk_difference = stats::median(rd),
      RD_lower = stats::quantile(rd, 0.025),
      RD_upper = stats::quantile(rd, 0.975),
      risk_ratio = stats::median(rr),
      RR_lower = stats::quantile(rr, 0.025),
      RR_upper = stats::quantile(rr, 0.975)
    )
  }

  eta_g <- eta0 + rep_draw(beta_G)
  eta_p_residual <- eta0 + rep_draw(beta_P)
  eta_p_total <- eta0 + rep_draw(beta_P + gamma_G_P * beta_G)
  age_increment <- AGE_REPORT_MULTIPLIER * (
    beta_age + gamma_P_age * beta_P + gamma_G_age * beta_G +
      gamma_P_age * gamma_G_P * beta_G
  )
  eta_age <- eta0 + rep_draw(age_increment)

  eta_no_fever <- eta0 - outer(beta_fever, dat$fever_bin, "*")
  eta_fever <- eta_no_fever + rep_draw(beta_fever)

  out <- dplyr::bind_rows(
    summarise_contrast("Tenfold higher gametocyte density", eta0, eta_g),
    summarise_contrast(
      "Tenfold higher asexual-stage density: residual conditional association",
      eta0, eta_p_residual
    ),
    summarise_contrast(
      "Tenfold higher asexual-stage density: total model-implied association",
      eta0, eta_p_total
    ),
    summarise_contrast("15-year higher host age: total model-implied association", eta0, eta_age),
    summarise_contrast("Fever versus no fever", eta_no_fever, eta_fever)
  )
  write_table(out, paste0("standardized_probability_contrasts_", fit_name, ".csv"))

  plot_dat <- out |>
    dplyr::mutate(
      contrast = factor(contrast, levels = rev(contrast))
    )
  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(risk_difference, contrast)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = RD_lower, xmax = RD_upper), height = 0.12
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(
      x = "Standardized difference in predicted mosquito infection probability",
      y = NULL
    ) + manuscript_theme()
  save_plot(p, "FigS_standardized_probability_contrasts", 7.5, 4.5)
  out
}

compute_loo <- function(fit, name) {
  ll <- as.matrix(fit$draws("log_lik", format = "matrix"))
  out <- loo::loo(ll)
  saveRDS(out, file.path(OUT_DIR, "loo", paste0("loo_", name, ".rds")))
  out
}

compare_loo <- function(named_loo_list) {
  cmp <- loo::loo_compare(named_loo_list)
  out <- as.data.frame(cmp) |>
    tibble::rownames_to_column("model")
  write_table(out, "TableS_LOO_model_comparison.csv")
  out
}

sensitivity_summary <- function(fit, scenario) {
  variables <- c("beta_G", "beta_P", "gamma_G_P", "beta_fever", "phi")
  purrr::map_dfr(variables, function(v) {
    x <- scalar_draw(fit, v)
    transform_to_or <- v %in% c("beta_G", "beta_P", "beta_fever")
    z <- if (transform_to_or) exp(x) else x
    tibble::tibble(
      scenario = scenario,
      parameter = v,
      scale = ifelse(transform_to_or, "Odds ratio", "Model scale"),
      estimate = stats::median(z),
      lower = stats::quantile(z, 0.025),
      upper = stats::quantile(z, 0.975)
    )
  })
}


# =============================================================================
# Final manuscript figure palette
# =============================================================================
PNTD_STUDY_COLOURS <- c(
  "Study 1" = "#A9C7DE",
  "Study 2" = "#E9C9A5",
  "Study 3" = "#B8D9C7"
)
PNTD_LINE_COLOURS <- c(
  "Study 1" = "#4E7FA3",
  "Study 2" = "#B47B4D",
  "Study 3" = "#5C9476"
)


# =============================================================================
# Observed-data and posterior manuscript figures generated during analysis
# =============================================================================
make_observed_data_figure <- function(dat) {
  common_violin <- function(yvar, ylab, percent = FALSE) {
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = study_label, y = .data[[yvar]], fill = study_label)
    ) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.62, colour = "grey55", linewidth = 0.35) +
      ggplot2::geom_boxplot(
        width = 0.15, outlier.shape = NA, fill = "white", colour = "grey35",
        linewidth = 0.35
      ) +
      ggplot2::scale_fill_manual(values = PNTD_STUDY_COLOURS, guide = "none") +
      ggplot2::labs(x = NULL, y = ylab) +
      manuscript_theme(9)
    if (percent) p <- p + ggplot2::scale_y_continuous(labels = scales::label_percent())
    p
  }

  plot_dat <- dat |>
    dplyr::mutate(
      observed_logP = log10(par + 0.1),
      observed_logG = log10(gam + 0.1)
    )
  old_dat <- dat
  dat <- plot_dat
  p_a <- common_violin(
    "observed_logP",
    "Asexual parasite density"
  )
  p_b <- common_violin(
    "observed_logG",
    "Gametocyte density"
  )
  p_c <- common_violin(
    "observed_proportion",
    "Proportion of infected mosquitoes (%)",
    percent = TRUE
  )
  p_d <- common_violin(
    "number_dissect",
    "Number of mosquitoes dissected"
  )
  dat <- old_dat

  combined <- (p_a + p_b) / (p_c + p_d) +
    patchwork::plot_annotation(tag_levels = "A")
  save_plot(combined, "Fig2_observed_microscopy_and_feeding_data", 8.0, 7.0)
  invisible(combined)
}

make_observed_density_relationship <- function(dat) {
  plot_dat <- dat |>
    dplyr::mutate(
      observed_logP = log10(par + 0.1),
      observed_logG = log10(gam + 0.1)
    )
  p <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(observed_logP, observed_logG, colour = observed_proportion)
  ) +
    ggplot2::geom_point(alpha = 0.82, size = 1.45) +
    ggplot2::facet_wrap(~study_label, nrow = 1) +
    ggplot2::scale_colour_gradientn(
      colours = c("#4F81A8", "#E8DCA8", "#C98273"),
      labels = scales::label_percent(),
      limits = c(0, 1),
      name = "Observed proportion\nof mosquitoes infected"
    ) +
    ggplot2::labs(
      x = "Asexual parasite density",
      y = "Gametocyte density"
    ) +
    manuscript_theme(9) +
    ggplot2::theme(legend.position = "right")
  save_plot(p, "Fig3_observed_density_relationship", 8.2, 4.5)
  invisible(p)
}

make_posterior_associations_figure <- function(fit) {
  summarise_or <- function(variable, label, multiplier = 1) {
    x <- exp(multiplier * scalar_draw(fit, variable))
    tibble::tibble(
      label = label,
      estimate = stats::median(x),
      lower = stats::quantile(x, 0.025),
      upper = stats::quantile(x, 0.975)
    )
  }

  panel_a_dat <- dplyr::bind_rows(
    summarise_or("beta_age", "Host age (15-year increase)", AGE_REPORT_MULTIPLIER),
    summarise_or("beta_fever", "Fever (yes vs no)"),
    summarise_or("beta_G", "Gametocyte density\n(10-fold increase)"),
    summarise_or("beta_P", "Asexual parasite density\n(10-fold increase; residual conditional)")
  ) |>
    dplyr::mutate(label = factor(label, levels = rev(label)))

  panel_b_dat <- dplyr::bind_rows(
    summarise_or("assoc_age_total", "Age: total model-implied association\n(15-year increase)", AGE_REPORT_MULTIPLIER),
    summarise_or("assoc_P_via_G", "Asexual parasite density:\npathway via gametocytes"),
    summarise_or("assoc_P_total", "Asexual parasite density:\ntotal model-implied association")
  ) |>
    dplyr::mutate(label = factor(label, levels = rev(label)))

  forest <- function(x, xlab) {
    ggplot2::ggplot(x, ggplot2::aes(estimate, label)) +
      ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = lower, xmax = upper), height = 0.12,
        colour = "#82A9C4", linewidth = 0.60
      ) +
      ggplot2::geom_point(colour = "#365F91", size = 2.1) +
      ggplot2::scale_x_log10(
        breaks = c(0.5, 1, 2, 4),
        labels = c("0.5", "1.0", "2.0", "4.0"),
        limits = c(0.45, 4.7)
      ) +
      ggplot2::labs(x = xlab, y = NULL) +
      manuscript_theme(8) +
      ggplot2::theme(legend.position = "none")
  }

  combined <- (
    forest(panel_a_dat, "Odds ratio (95% CrI)") |
      forest(panel_b_dat, "Model-implied odds ratio (95% CrI)")
  ) +
    patchwork::plot_annotation(tag_levels = "A")
  save_plot(combined, "Fig4_posterior_associations", 8.2, 4.4)
  invisible(list(conditional = panel_a_dat, pathway = panel_b_dat))
}

make_relationship_figure <- function(fit, dat, fit_name, max_draws = 2000L) {
  d <- posterior::as_draws_df(fit$draws(
    variables = c(
      "alpha", "beta_G", "beta_P", "beta_age", "beta_fever",
      "beta_study", "phi"
    )
  ))
  keep <- subsample_rows(nrow(d), max_draws, SEED + 20L)
  d <- d[keep, , drop = FALSE]
  nd <- nrow(d)
  J <- nlevels(dat$study_label)

  study_effect <- function(j) {
    if (j == 1L) rep(0, nd) else d[[paste0("beta_study[", j, "]")]]
  }

  summarize_grid <- function(x, eta_matrix, study_name, x_name, m_ref = 30L) {
    pmat <- plogis(eta_matrix)
    out <- tibble::tibble(
      x = x,
      study_label = factor(study_name, levels = levels(dat$study_label)),
      median = apply(pmat, 2, stats::median),
      lower = apply(pmat, 2, stats::quantile, probs = 0.025),
      upper = apply(pmat, 2, stats::quantile, probs = 0.975),
      x_name = x_name
    )

    set.seed(SEED + 21L)
    out$predictive_lower <- vapply(seq_along(x), function(k) {
      p <- pmat[, k]
      theta <- stats::rbeta(nd, p * d$phi, (1 - p) * d$phi)
      stats::quantile(stats::rbinom(nd, size = m_ref, prob = theta) / m_ref, 0.025)
    }, numeric(1))
    out$predictive_upper <- vapply(seq_along(x), function(k) {
      p <- pmat[, k]
      theta <- stats::rbeta(nd, p * d$phi, (1 - p) * d$phi)
      stats::quantile(stats::rbinom(nd, size = m_ref, prob = theta) / m_ref, 0.975)
    }, numeric(1))
    out
  }

  p_grid <- seq(min(log10(dat$par[dat$par > 0])), max(log10(dat$par[dat$par > 0])), length.out = 80)
  g_grid <- seq(min(log10(dat$gam[dat$gam > 0])), max(log10(dat$gam[dat$gam > 0])), length.out = 80)
  age_grid <- seq(stats::quantile(dat$age, 0.02), stats::quantile(dat$age, 0.98), length.out = 80)
  age_center <- unique(dat$analysis_age_center)[1]

  g_curves <- purrr::map_dfr(seq_len(J), function(j) {
    eta <- outer(d$beta_G, g_grid - CENTER_G, "*") +
      matrix(d$alpha + study_effect(j), nrow = nd, ncol = length(g_grid))
    summarize_grid(g_grid, eta, levels(dat$study_label)[j], "gametocyte")
  })

  p_curves <- purrr::map_dfr(seq_len(J), function(j) {
    eta <- outer(d$beta_P, p_grid - CENTER_P, "*") +
      matrix(d$alpha + study_effect(j), nrow = nd, ncol = length(p_grid))
    summarize_grid(p_grid, eta, levels(dat$study_label)[j], "asexual")
  })

  age_curves <- purrr::map_dfr(seq_len(J), function(j) {
    age10_grid <- (age_grid - age_center) / AGE_INTERNAL_SCALE_YEARS
    eta <- outer(d$beta_age, age10_grid, "*") +
      matrix(d$alpha + study_effect(j), nrow = nd, ncol = length(age_grid))
    summarize_grid(age_grid, eta, levels(dat$study_label)[j], "age")
  })

  fever_curves <- purrr::map_dfr(seq_len(J), function(j) {
    x <- c(0, 1)
    eta <- outer(d$beta_fever, x, "*") +
      matrix(d$alpha + study_effect(j), nrow = nd, ncol = 2L)
    summarize_grid(x, eta, levels(dat$study_label)[j], "fever")
  })

  plot_curve <- function(curve, xlab) {
    ggplot2::ggplot(curve, ggplot2::aes(x, median)) +
      # The main predictor figure displays uncertainty in the mean probability only.
      # Participant-level beta-binomial predictive uncertainty is shown in Fig 6.
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        fill = "#A9C7DE", alpha = 0.42
      ) +
      ggplot2::geom_line(colour = "#3F6F96", linewidth = 0.85) +
      ggplot2::facet_wrap(~study_label, nrow = 1) +
      ggplot2::scale_y_continuous(
        labels = scales::label_percent(accuracy = 1),
        breaks = seq(0, 1, 0.25)
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        x = xlab,
        y = "Predicted proportion of infected mosquitoes (%)"
      ) +
      manuscript_theme(9) +
      ggplot2::theme(legend.position = "none")
  }

  panel_a <- plot_curve(g_curves, "Gametocyte density")
  panel_b <- plot_curve(p_curves, "Asexual parasite density")
  panel_c <- plot_curve(age_curves, "Host age (years)") +
    ggplot2::geom_vline(
      xintercept = AGE_CENTER_YEARS,
      linetype = 3, colour = "grey55", linewidth = 0.45
    )

  fever_plot_dat <- fever_curves |>
    dplyr::mutate(
      fever_label = factor(ifelse(x == 1, "Fever", "No fever"),
                           levels = c("No fever", "Fever"))
    )

  panel_d <- ggplot2::ggplot(fever_plot_dat, ggplot2::aes(fever_label, median)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 0.10, colour = "#3F6F96", linewidth = 0.75
    ) +
    ggplot2::geom_point(colour = "#3F6F96", size = 2.1) +
    ggplot2::facet_wrap(~study_label, nrow = 1) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      breaks = seq(0, 1, 0.25)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Fever status",
      y = "Predicted proportion of infected mosquitoes (%)"
    ) +
    manuscript_theme(9) +
    ggplot2::theme(legend.position = "none")

  combined <- panel_a / panel_b / panel_c / panel_d +
    patchwork::plot_annotation(tag_levels = "A")
  save_plot(combined, "Fig5_modelled_predictor_relationships", 7.9, 10.5)

  write_table(g_curves, paste0("relationship_gametocyte_", fit_name, ".csv"))
  write_table(p_curves, paste0("relationship_asexual_", fit_name, ".csv"))
  write_table(age_curves, paste0("relationship_age_", fit_name, ".csv"))
  write_table(fever_curves, paste0("relationship_fever_", fit_name, ".csv"))
  invisible(list(g = g_curves, p = p_curves, age = age_curves, fever = fever_curves))
}

make_posterior_predictive_outputs <- function(fit, dat, fit_name) {
  N <- nrow(dat)
  yrep <- indexed_draw_matrix(fit, "y_rep", N)
  pmat <- indexed_draw_matrix(fit, "p_infect", N)
  prop_rep <- sweep(yrep, 2, dat$number_dissect, "/")

  pred <- tibble::tibble(
    id = dat$id,
    study_label = dat$study_label,
    observed = dat$observed_proportion,
    predicted_median = apply(pmat, 2, stats::median),
    predicted_mean_lower = apply(pmat, 2, stats::quantile, probs = 0.025),
    predicted_mean_upper = apply(pmat, 2, stats::quantile, probs = 0.975),
    predictive_lower = apply(prop_rep, 2, stats::quantile, probs = 0.025),
    predictive_upper = apply(prop_rep, 2, stats::quantile, probs = 0.975)
  ) |>
    dplyr::mutate(
      predictive_width = predictive_upper - predictive_lower,
      covered = observed >= predictive_lower & observed <= predictive_upper
    )
  write_table(pred, paste0("posterior_predictions_", fit_name, ".csv"))

  panel_a <- ggplot2::ggplot(
    pred,
    ggplot2::aes(observed, predicted_median, colour = study_label)
  ) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = predictive_lower, ymax = predictive_upper),
      alpha = 0.12, linewidth = 0.24
    ) +
    ggplot2::geom_point(alpha = 0.76, size = 1.55) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = 2,
      colour = "grey50", linewidth = 0.45
    ) +
    ggplot2::facet_wrap(~study_label, nrow = 1) +
    ggplot2::scale_colour_manual(values = PNTD_LINE_COLOURS, guide = "none") +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, 0.25)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, 0.25)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Observed proportion of mosquitoes infected (%)",
      y = "Posterior median infection probability (%)"
    ) +
    manuscript_theme(8)

  # Five within-study groups are used instead of deciles so Study 3 retains
  # enough participants per group for an interpretable calibration display.
  cal <- pred |>
    dplyr::mutate(y = dat$number_inf, m = dat$number_dissect) |>
    dplyr::group_by(study_label) |>
    dplyr::mutate(prediction_group = dplyr::ntile(predicted_median, 5L)) |>
    dplyr::group_by(study_label, prediction_group) |>
    dplyr::summarise(
      predicted = mean(predicted_median),
      observed = sum(y) / sum(m),
      n = dplyr::n(),
      .groups = "drop"
    )

  panel_b <- ggplot2::ggplot(
    cal,
    ggplot2::aes(predicted, observed, colour = study_label)
  ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = 2,
      colour = "grey50", linewidth = 0.45
    ) +
    ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.86) +
    ggplot2::scale_colour_manual(values = PNTD_LINE_COLOURS, guide = "none") +
    ggplot2::scale_size_continuous(range = c(2.0, 5.0)) +
    ggplot2::facet_wrap(~study_label, nrow = 1) +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, 0.25)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, 0.25)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Mean predicted infection probability (%)",
      y = "Observed proportion of mosquitoes infected (%)",
      size = "Participants"
    ) +
    manuscript_theme(8)

  coverage <- pred |>
    dplyr::group_by(study_label) |>
    dplyr::summarise(
      coverage = mean(covered),
      median_predictive_width = stats::median(predictive_width),
      .groups = "drop"
    )
  write_table(cal, paste0("calibration_groups_", fit_name, ".csv"))
  write_table(coverage, paste0("predictive_coverage_", fit_name, ".csv"))

  combined <- panel_a / panel_b + patchwork::plot_annotation(tag_levels = "A")
  save_plot(combined, "Fig6_model_fit_and_calibration", 8.3, 6.0)
  invisible(list(predictions = pred, calibration = cal, coverage = coverage))
}


# =============================================================================
# Supplementary diagnostic figures
# =============================================================================
make_prior_predictive_plot <- function(prior_fit, dat) {
  yrep <- indexed_draw_matrix(prior_fit, "y_rep", nrow(dat))
  keep <- subsample_rows(nrow(yrep), 1500L, SEED + 10L)
  proportions <- sweep(yrep[keep, , drop = FALSE], 2, dat$number_dissect, "/")
  prior_means <- rowMeans(proportions)
  observed_mean <- mean(dat$observed_proportion)
  p <- ggplot2::ggplot(tibble::tibble(prior_mean = prior_means), ggplot2::aes(prior_mean)) +
    ggplot2::geom_histogram(bins = 40, fill = "grey75", color = "white") +
    ggplot2::geom_vline(xintercept = observed_mean, linewidth = 0.8, linetype = 2) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(
      x = "Prior-predictive mean proportion of mosquitoes infected",
      y = "Prior-predictive draws"
    ) +
    manuscript_theme()
  save_plot(p, "FigS1_prior_predictive", 6.2, 4.5)
}

make_latent_observed_diagnostic <- function(fit, dat, fit_name) {
  N <- nrow(dat)
  latentP <- indexed_draw_matrix(fit, "logP", N)
  latentG <- indexed_draw_matrix(fit, "logG", N)

  make_stage <- function(mat, observed, positive, stage_name) {
    keep <- which(positive)
    tibble::tibble(
      id = dat$id[keep],
      study_label = dat$study_label[keep],
      stage = stage_name,
      observed_log10 = log10(observed[keep]),
      latent_median = apply(mat[, keep, drop = FALSE], 2, stats::median),
      latent_lower = apply(mat[, keep, drop = FALSE], 2, stats::quantile, probs = 0.025),
      latent_upper = apply(mat[, keep, drop = FALSE], 2, stats::quantile, probs = 0.975)
    )
  }

  plot_dat <- dplyr::bind_rows(
    make_stage(latentP, dat$par, dat$par > 0, "Asexual-stage density"),
    make_stage(latentG, dat$gam, dat$gam > 0, "Gametocyte density")
  )
  write_table(plot_dat, paste0("latent_vs_observed_", fit_name, ".csv"))

  diagnostics <- plot_dat |>
    dplyr::group_by(stage, study_label) |>
    dplyr::summarise(
      n = dplyr::n(),
      spearman_correlation = stats::cor(
        observed_log10, latent_median, method = "spearman"
      ),
      median_absolute_difference = stats::median(
        abs(observed_log10 - latent_median)
      ),
      .groups = "drop"
    )
  write_table(diagnostics, paste0("latent_anchoring_diagnostics_", fit_name, ".csv"))

  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(observed_log10, latent_median)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = latent_lower, ymax = latent_upper), alpha = 0.18
    ) +
    ggplot2::geom_point(alpha = 0.55, size = 1.1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::facet_grid(stage ~ study_label, scales = "free") +
    ggplot2::labs(
      x = "Observed density (log10 scale)",
      y = "Posterior latent density (log10 scale)"
    ) + manuscript_theme(9)
  save_plot(p, "FigS2_latent_vs_observed", 8.2, 6.2)
  invisible(list(data = plot_dat, diagnostics = diagnostics))
}

make_sensitivity_plot <- function(sensitivity_data) {
  plot_data <- sensitivity_data |>
    dplyr::filter(parameter %in% c("beta_G", "beta_P", "gamma_G_P", "phi")) |>
    dplyr::mutate(
      parameter_label = dplyr::recode(
        parameter,
        beta_G = "Gametocyte–infectivity association (OR)",
        beta_P = "Residual asexual-stage–infectivity association (OR)",
        gamma_G_P = "Asexual-stage–gametocyte coefficient",
        phi = "Beta-binomial precision"
      ),
      scenario = factor(scenario, levels = rev(unique(scenario)))
    )

  reference_data <- tibble::tibble(
    parameter_label = c(
      "Gametocyte–infectivity association (OR)",
      "Residual asexual-stage–infectivity association (OR)",
      "Asexual-stage–gametocyte coefficient"
    ),
    reference = c(1, 1, 0)
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(estimate, scenario)) +
    ggplot2::geom_vline(
      data = reference_data,
      ggplot2::aes(xintercept = reference),
      linetype = 2, colour = "grey60", linewidth = 0.4,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = lower, xmax = upper), height = 0.12
    ) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~parameter_label, scales = "free_x", ncol = 2) +
    ggplot2::labs(x = "Posterior median and 95% credible interval", y = NULL) +
    manuscript_theme(8)
  save_plot(p, "FigS4_sensitivity_analysis", 9.2, 7.2)
  invisible(p)
}
