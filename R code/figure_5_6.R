###############################################################################
#  FIGURES 5 AND 6
# Microscopy-based joint Bayesian latent-variable model
#
# 
#   Fig 5: observed participant points + model-estimated mean relationship +
#          95% CrI for the mean. Posterior predictive limits are intentionally
#          omitted from Fig 5 to keep the biological relationships clear.
#   Fig 6: participant-level model fit, separation of mean uncertainty from
#          posterior predictive uncertainty, and study-specific calibration.
#

###############################################################################

# =============================================================================
# 1. PACKAGES
# =============================================================================
required_packages <- c(
  "cmdstanr", "posterior", "dplyr", "readr", "stringr", "tibble",
  "ggplot2", "patchwork", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install these packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# =============================================================================
# 2. USER SETTINGS
# =============================================================================
# Recommended use:
#   A. Run this script after the main analysis while primary_fit and primary_dat
#      are still in memory; OR
#   B. Edit PROJECT_DIR below so the script loads the saved fit and raw data.

get_project_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/")
    script_dir <- dirname(script_path)
    # Figure scripts live one directory below the repository root.
    if (basename(script_dir) == "figures") return(dirname(script_dir))
    return(script_dir)
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
PROJECT_DIR <- get_project_dir()


# Saved primary CmdStanR fit produced by the final analysis pipeline.
FIT_RDS <- file.path(
  PROJECT_DIR,
  "results", "fits", "primary", "primary.rds"
)

# Raw microscopy analysis dataset.
DATA_CSV <- file.path(
  PROJECT_DIR,
  "data", "data_pvivax.csv"
)

# Output folder for final submission figures and supporting tables.
FIGURE_OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "results", "final_figures"
)
TABLE_OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "results", "final_tables"
)

dir.create(FIGURE_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Model centring constants used in the final Stan model.
CENTER_P <- log10(3000)
CENTER_G <- log10(200)
AGE_CENTER_YEARS <- 15
AGE_INTERNAL_SCALE_YEARS <- 10

# Figure settings.
SEED <- 20260713
MAX_POSTERIOR_DRAWS <- 2000L
GRID_POINTS <- 100L

# Low-density limits displayed in Fig 5.
# Gametocyte -1 shows the descriptive location of recorded microscopy zeros
# under log10(x + 0.1). The model curve itself remains a latent-density curve.
GAMETOCYTE_X_MIN <- -1
ASEXUAL_X_MIN <- 0

# Gentle PLOS-compatible colours.
LINE_COLOUR <- "#356F9F"
RIBBON_COLOUR <- "#B9D4E8"
OBSERVED_POINT_COLOUR <- "#595959"
PREDICTIVE_LINE_COLOUR <- "#8F8F8F"
STUDY_COLOURS <- c(
  "Study 1" = "#5B8DB8",
  "Study 2" = "#C68B57",
  "Study 3" = "#6FA08A"
)

# =============================================================================
# 3. LOAD OR PREPARE THE FIT AND DATA
# =============================================================================
prepare_plot_data <- function(raw) {
  dat <- raw |>
    dplyr::mutate(
      id = as.character(id),
      study_original = as.integer(study),
      fever_text = stringr::str_to_lower(
        stringr::str_trim(as.character(fever))
      ),
      fever_bin = dplyr::case_when(
        fever_text %in% c("yes", "y", "1", "true") ~ 1L,
        fever_text %in% c("no", "n", "0", "false") ~ 0L,
        TRUE ~ NA_integer_
      )
    ) |>
    dplyr::filter(study_original %in% c(1L, 2L, 3L))

  required <- c(
    "id", "study_original", "par", "gam", "number_dissect",
    "number_inf", "age", "fever_bin"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop("Missing required variables: ", paste(missing, collapse = ", "))
  }

  # Remove exact duplicate participant records only. Conflicting duplicates stop.
  duplicate_check <- dat |>
    dplyr::group_by(study_original, id) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_par = dplyr::n_distinct(par, na.rm = FALSE),
      n_gam = dplyr::n_distinct(gam, na.rm = FALSE),
      n_m = dplyr::n_distinct(number_dissect, na.rm = FALSE),
      n_y = dplyr::n_distinct(number_inf, na.rm = FALSE),
      n_age = dplyr::n_distinct(age, na.rm = FALSE),
      n_fever = dplyr::n_distinct(fever_bin, na.rm = FALSE),
      .groups = "drop"
    ) |>
    dplyr::filter(n > 1L)

  if (nrow(duplicate_check) > 0L) {
    conflict_columns <- c(
      "n_par", "n_gam", "n_m", "n_y", "n_age", "n_fever"
    )
    if (any(as.matrix(duplicate_check[conflict_columns]) > 1L)) {
      stop("Conflicting duplicate participant records were detected.")
    }
    dat <- dat |>
      dplyr::distinct(study_original, id, .keep_all = TRUE)
  }

  if (anyNA(dat[c(
    "par", "gam", "number_dissect", "number_inf", "age", "fever_bin"
  )])) {
    stop("Missing values remain in the core plotting variables.")
  }
  if (any(dat$par < 0 | dat$gam < 0)) {
    stop("Microscopy densities cannot be negative.")
  }
  if (any(dat$number_dissect < 1L)) {
    stop("number_dissect must be at least one.")
  }
  if (any(dat$number_inf < 0L | dat$number_inf > dat$number_dissect)) {
    stop("number_inf must be between zero and number_dissect.")
  }

  dat |>
    dplyr::mutate(
      study_label = factor(
        paste0("Study ", study_original),
        levels = c("Study 1", "Study 2", "Study 3")
      ),
      study_id = as.integer(study_label),
      age10 = (age - AGE_CENTER_YEARS) / AGE_INTERNAL_SCALE_YEARS,
      observed_proportion = number_inf / number_dissect,
      observed_logP_display = log10(par + 0.1),
      observed_logG_display = log10(gam + 0.1),
      fever_label = factor(
        ifelse(fever_bin == 1L, "Fever", "No fever"),
        levels = c("No fever", "Fever")
      )
    ) |>
    dplyr::arrange(study_id, id)
}

if (!exists("primary_fit", inherits = TRUE)) {
  if (!file.exists(FIT_RDS)) {
    stop(
      "primary_fit is not in memory and the fit file was not found:\n",
      FIT_RDS,
      "\nRun the main model first or edit FIT_RDS."
    )
  }
  primary_fit <- readRDS(FIT_RDS)
}

if (!exists("primary_dat", inherits = TRUE)) {
  if (!file.exists(DATA_CSV)) {
    stop(
      "primary_dat is not in memory and the data file was not found:\n",
      DATA_CSV,
      "\nRun the main model first or edit DATA_CSV."
    )
  }
  primary_dat <- prepare_plot_data(
    readr::read_csv(DATA_CSV, show_col_types = FALSE)
  )
} else {
  # Ensure all display variables exist when primary_dat came from the main code.
  primary_dat <- primary_dat |>
    dplyr::mutate(
      observed_proportion = number_inf / number_dissect,
      observed_logP_display = log10(par + 0.1),
      observed_logG_display = log10(gam + 0.1),
      fever_label = factor(
        ifelse(fever_bin == 1L, "Fever", "No fever"),
        levels = c("No fever", "Fever")
      )
    )
}

# =============================================================================
# 4. GENERAL HELPERS
# =============================================================================
manuscript_theme <- function(base_size = 9) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_rect(
        fill = "grey96", colour = "grey82", linewidth = 0.35
      ),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      axis.title = ggplot2::element_text(face = "plain"),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.35),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(5.5, 7, 5.5, 7)
    )
}

save_manuscript_plot <- function(plot, stem, width, height) {
  pdf_file <- file.path(FIGURE_OUTPUT_DIR, paste0(stem, ".pdf"))
  tif_file <- file.path(FIGURE_OUTPUT_DIR, paste0(stem, ".tiff"))

  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = tif_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )

  message("Saved: ", pdf_file)
  message("Saved: ", tif_file)
}

indexed_draw_matrix <- function(fit, variable, expected_columns = NULL) {
  mat <- as.matrix(fit$draws(variable, format = "matrix"))
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1L)

  cn <- colnames(mat)
  if (!is.null(cn) && ncol(mat) > 1L) {
    pattern <- paste0("^", variable, "\\[(\\d+)\\]$")
    idx <- suppressWarnings(as.integer(sub(pattern, "\\1", cn)))
    if (all(is.finite(idx))) {
      mat <- mat[, order(idx), drop = FALSE]
    }
  }

  if (!is.null(expected_columns) && ncol(mat) != expected_columns) {
    stop(
      variable, " has ", ncol(mat),
      " columns; expected ", expected_columns, "."
    )
  }
  mat
}

scalar_draw <- function(fit, variable) {
  mat <- indexed_draw_matrix(fit, variable, expected_columns = 1L)
  as.numeric(mat[, 1L])
}

subsample_draws <- function(n, max_draws = MAX_POSTERIOR_DRAWS, seed = SEED) {
  if (n <= max_draws) return(seq_len(n))
  set.seed(seed)
  sort(sample.int(n, max_draws))
}

clip_probability <- function(p) {
  pmin(pmax(p, 1e-8), 1 - 1e-8)
}

summarise_curve <- function(
    x,
    eta_matrix,
    phi_draws,
    study_label,
    m_reference,
    predictive_seed) {

  p_matrix <- clip_probability(stats::plogis(eta_matrix))
  n_draws <- nrow(p_matrix)
  n_grid <- ncol(p_matrix)

  result <- tibble::tibble(
    x = x,
    study_label = factor(
      study_label,
      levels = levels(primary_dat$study_label)
    ),
    median = apply(p_matrix, 2, stats::median),
    mean_lower = apply(p_matrix, 2, stats::quantile, probs = 0.025),
    mean_upper = apply(p_matrix, 2, stats::quantile, probs = 0.975),
    m_reference = as.integer(m_reference)
  )

  set.seed(predictive_seed)
  predictive_lower <- numeric(n_grid)
  predictive_upper <- numeric(n_grid)

  for (k in seq_len(n_grid)) {
    p_k <- p_matrix[, k]
    theta_k <- stats::rbeta(
      n_draws,
      shape1 = p_k * phi_draws,
      shape2 = (1 - p_k) * phi_draws
    )
    replicated_proportion <- stats::rbinom(
      n_draws,
      size = as.integer(m_reference),
      prob = theta_k
    ) / as.integer(m_reference)

    predictive_lower[k] <- stats::quantile(
      replicated_proportion, 0.025, names = FALSE
    )
    predictive_upper[k] <- stats::quantile(
      replicated_proportion, 0.975, names = FALSE
    )
  }

  result$predictive_lower <- predictive_lower
  result$predictive_upper <- predictive_upper
  result
}

# =============================================================================
# 5. BUILD STUDY-SPECIFIC PREDICTION CURVES FOR FIGURE 5
# =============================================================================
build_fig5_prediction_data <- function(
    fit,
    dat,
    max_draws = MAX_POSTERIOR_DRAWS) {

  J <- nlevels(dat$study_label)

  alpha <- scalar_draw(fit, "alpha")
  beta_G <- scalar_draw(fit, "beta_G")
  beta_P <- scalar_draw(fit, "beta_P")
  beta_age <- scalar_draw(fit, "beta_age")
  beta_fever <- scalar_draw(fit, "beta_fever")
  phi <- scalar_draw(fit, "phi")

  mu_P <- scalar_draw(fit, "mu_P")
  mu_G <- scalar_draw(fit, "mu_G")
  gamma_P_age <- scalar_draw(fit, "gamma_P_age")
  gamma_G_age <- scalar_draw(fit, "gamma_G_age")
  gamma_G_P <- scalar_draw(fit, "gamma_G_P")

  beta_study <- indexed_draw_matrix(fit, "beta_study", J)
  delta_P_study <- indexed_draw_matrix(fit, "delta_P_study", J)
  delta_G_study <- indexed_draw_matrix(fit, "delta_G_study", J)

  n_all <- length(alpha)
  keep <- subsample_draws(n_all, max_draws, SEED + 10L)

  alpha <- alpha[keep]
  beta_G <- beta_G[keep]
  beta_P <- beta_P[keep]
  beta_age <- beta_age[keep]
  beta_fever <- beta_fever[keep]
  phi <- phi[keep]
  mu_P <- mu_P[keep]
  mu_G <- mu_G[keep]
  gamma_P_age <- gamma_P_age[keep]
  gamma_G_age <- gamma_G_age[keep]
  gamma_G_P <- gamma_G_P[keep]
  beta_study <- beta_study[keep, , drop = FALSE]
  delta_P_study <- delta_P_study[keep, , drop = FALSE]
  delta_G_study <- delta_G_study[keep, , drop = FALSE]

  n_draws <- length(alpha)

  # Density ranges. The y-axis starts at 0%; curves are not constrained to 0%.
  positive_g <- dat$gam[dat$gam > 0]
  positive_p <- dat$par[dat$par > 0]
  if (length(positive_g) == 0L || length(positive_p) == 0L) {
    stop("Positive microscopy densities are required to construct density grids.")
  }

  g_max <- ceiling(max(log10(positive_g), na.rm = TRUE) * 2) / 2
  p_max <- ceiling(max(log10(positive_p), na.rm = TRUE) * 2) / 2

  g_grid <- seq(GAMETOCYTE_X_MIN, g_max, length.out = GRID_POINTS)
  p_grid <- seq(ASEXUAL_X_MIN, p_max, length.out = GRID_POINTS)
  age_grid <- seq(
    max(0, floor(min(dat$age, na.rm = TRUE))),
    ceiling(max(dat$age, na.rm = TRUE)),
    length.out = GRID_POINTS
  )

  g_curves <- list()
  p_curves <- list()
  age_curves <- list()
  fever_curves <- list()

  for (j in seq_len(J)) {
    study_name <- levels(dat$study_label)[j]
    m_reference <- as.integer(round(stats::median(
      dat$number_dissect[dat$study_id == j],
      na.rm = TRUE
    )))

    # Study-specific expected latent densities for a participant aged 15 years.
    logP_reference <- mu_P + delta_P_study[, j]
    logG_reference <- mu_G + delta_G_study[, j] +
      gamma_G_P * (logP_reference - CENTER_P)

    # Panel A: vary latent gametocyte density; hold latent asexual density at
    # its study-specific expected value, age at 15, and fever at no fever.
    eta_G <- matrix(
      alpha + beta_study[, j] +
        beta_P * (logP_reference - CENTER_P),
      nrow = n_draws,
      ncol = length(g_grid)
    ) + outer(beta_G, g_grid - CENTER_G, "*")

    g_curves[[j]] <- summarise_curve(
      x = g_grid,
      eta_matrix = eta_G,
      phi_draws = phi,
      study_label = study_name,
      m_reference = m_reference,
      predictive_seed = SEED + 100L + j
    )

    # Panel B: residual conditional association. Vary latent asexual density
    # while holding latent gametocyte density at its study-specific reference.
    eta_P <- matrix(
      alpha + beta_study[, j] +
        beta_G * (logG_reference - CENTER_G),
      nrow = n_draws,
      ncol = length(p_grid)
    ) + outer(beta_P, p_grid - CENTER_P, "*")

    p_curves[[j]] <- summarise_curve(
      x = p_grid,
      eta_matrix = eta_P,
      phi_draws = phi,
      study_label = study_name,
      m_reference = m_reference,
      predictive_seed = SEED + 200L + j
    )

    # Panel C: total model-implied age relationship. Age is allowed to change
    # the expected latent asexual and gametocyte densities according to the
    # fitted latent-process equations. Age remains continuous.
    age10_grid <- (age_grid - AGE_CENTER_YEARS) / AGE_INTERNAL_SCALE_YEARS

    logP_age <- outer(gamma_P_age, age10_grid, "*") +
      matrix(
        mu_P + delta_P_study[, j],
        nrow = n_draws,
        ncol = length(age_grid)
      )

    logG_age <- matrix(
      mu_G + delta_G_study[, j],
      nrow = n_draws,
      ncol = length(age_grid)
    ) +
      gamma_G_P * (logP_age - CENTER_P) +
      outer(gamma_G_age, age10_grid, "*")

    eta_age <- matrix(
      alpha + beta_study[, j],
      nrow = n_draws,
      ncol = length(age_grid)
    ) +
      beta_G * (logG_age - CENTER_G) +
      beta_P * (logP_age - CENTER_P) +
      outer(beta_age, age10_grid, "*")

    age_curves[[j]] <- summarise_curve(
      x = age_grid,
      eta_matrix = eta_age,
      phi_draws = phi,
      study_label = study_name,
      m_reference = m_reference,
      predictive_seed = SEED + 300L + j
    )

    # Panel D: fever contrast at study-specific expected latent densities and
    # age 15 years.
    eta_no_fever <- alpha + beta_study[, j] +
      beta_G * (logG_reference - CENTER_G) +
      beta_P * (logP_reference - CENTER_P)
    eta_fever <- eta_no_fever + beta_fever
    eta_fever_matrix <- cbind(eta_no_fever, eta_fever)

    fever_curves[[j]] <- summarise_curve(
      x = c(0, 1),
      eta_matrix = eta_fever_matrix,
      phi_draws = phi,
      study_label = study_name,
      m_reference = m_reference,
      predictive_seed = SEED + 400L + j
    ) |>
      dplyr::mutate(
        fever_label = factor(
          ifelse(x == 1, "Fever", "No fever"),
          levels = c("No fever", "Fever")
        )
      )
  }

  list(
    gametocyte = dplyr::bind_rows(g_curves),
    asexual = dplyr::bind_rows(p_curves),
    age = dplyr::bind_rows(age_curves),
    fever = dplyr::bind_rows(fever_curves)
  )
}

# =============================================================================
# 6. FINAL FIGURE 5
# =============================================================================
make_final_fig5 <- function(
    fit = primary_fit,
    dat = primary_dat,
    max_draws = MAX_POSTERIOR_DRAWS) {

  curves <- build_fig5_prediction_data(fit, dat, max_draws)

  # Fig 5 is intentionally restricted to the observed participant points,
  # the model-estimated mean relationship, and its 95% credible interval.
  # Participant-level posterior predictive uncertainty is shown in Fig 6.
  continuous_panel <- function(
      curve_data,
      observed_data,
      observed_x,
      x_label,
      x_limits,
      x_breaks = waiver(),
      add_age_reference = FALSE) {

    p <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(
        data = curve_data,
        ggplot2::aes(x = x, ymin = mean_lower, ymax = mean_upper),
        fill = RIBBON_COLOUR,
        alpha = 0.55
      ) +
      ggplot2::geom_point(
        data = observed_data,
        mapping = ggplot2::aes(
          x = .data[[observed_x]],
          y = observed_proportion
        ),
        colour = OBSERVED_POINT_COLOUR,
        size = 1.15,
        alpha = 0.24
      ) +
      ggplot2::geom_line(
        data = curve_data,
        ggplot2::aes(x = x, y = median),
        colour = LINE_COLOUR,
        linewidth = 0.90
      ) +
      ggplot2::facet_wrap(~study_label, nrow = 1) +
      ggplot2::scale_x_continuous(
        limits = x_limits,
        breaks = x_breaks,
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25),
        labels = scales::label_percent(accuracy = 1),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::labs(
        x = x_label,
        y = "Proportion of mosquitoes infected (%)"
      ) +
      manuscript_theme(9) +
      ggplot2::theme(legend.position = "none")

    if (isTRUE(add_age_reference)) {
      p <- p +
        ggplot2::geom_vline(
          xintercept = AGE_CENTER_YEARS,
          linetype = 2,
          colour = "grey45",
          linewidth = 0.45
        )
    }
    p
  }

  g_limits <- range(curves$gametocyte$x)
  p_limits <- range(curves$asexual$x)
  age_limits <- range(curves$age$x)

  panel_a <- continuous_panel(
    curve_data = curves$gametocyte,
    observed_data = dat,
    observed_x = "observed_logG_display",
    x_label = "Gametocyte density (log10 scale)",
    x_limits = g_limits,
    x_breaks = seq(floor(g_limits[1]), ceiling(g_limits[2]), by = 1)
  )

  panel_b <- continuous_panel(
    curve_data = curves$asexual,
    observed_data = dat,
    observed_x = "observed_logP_display",
    x_label = "Asexual parasite density (log10 scale)",
    x_limits = p_limits,
    x_breaks = seq(floor(p_limits[1]), ceiling(p_limits[2]), by = 1)
  )

  panel_c <- continuous_panel(
    curve_data = curves$age,
    observed_data = dat,
    observed_x = "age",
    x_label = "Host age (years)",
    x_limits = age_limits,
    x_breaks = scales::breaks_pretty(n = 6),
    add_age_reference = TRUE
  )

  panel_d <- ggplot2::ggplot() +
    ggplot2::geom_jitter(
      data = dat,
      ggplot2::aes(x = fever_label, y = observed_proportion),
      width = 0.08,
      height = 0,
      colour = OBSERVED_POINT_COLOUR,
      size = 1.15,
      alpha = 0.25
    ) +
    ggplot2::geom_errorbar(
      data = curves$fever,
      ggplot2::aes(
        x = fever_label,
        ymin = mean_lower,
        ymax = mean_upper
      ),
      width = 0.10,
      colour = LINE_COLOUR,
      linewidth = 0.75
    ) +
    ggplot2::geom_point(
      data = curves$fever,
      ggplot2::aes(x = fever_label, y = median),
      colour = LINE_COLOUR,
      size = 2.25
    ) +
    ggplot2::facet_wrap(~study_label, nrow = 1) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      x = "Fever status",
      y = "Proportion of mosquitoes infected (%)"
    ) +
    manuscript_theme(9) +
    ggplot2::theme(legend.position = "none")

  fig5 <- panel_a / panel_b / panel_c / panel_d +
    patchwork::plot_annotation(tag_levels = "A")

  save_manuscript_plot(
    fig5,
    "Fig5_modelled_predictor_relationships_FINAL",
    width = 8.4,
    height = 10.8
  )

  readr::write_csv(
    curves$gametocyte,
    file.path(TABLE_OUTPUT_DIR, "Fig5_gametocyte_curve.csv")
  )
  readr::write_csv(
    curves$asexual,
    file.path(TABLE_OUTPUT_DIR, "Fig5_asexual_curve.csv")
  )
  readr::write_csv(
    curves$age,
    file.path(TABLE_OUTPUT_DIR, "Fig5_age_curve.csv")
  )
  readr::write_csv(
    curves$fever,
    file.path(TABLE_OUTPUT_DIR, "Fig5_fever_contrast.csv")
  )

  invisible(list(plot = fig5, curves = curves))
}

# =============================================================================
# 7. FINAL FIGURE 6
# =============================================================================
make_final_fig6 <- function(fit = primary_fit, dat = primary_dat) {
  N <- nrow(dat)

  p_matrix <- indexed_draw_matrix(fit, "p_infect", N)
  yrep_matrix <- indexed_draw_matrix(fit, "y_rep", N)
  replicated_proportion <- sweep(
    yrep_matrix,
    2,
    dat$number_dissect,
    "/"
  )

  pred <- tibble::tibble(
    id = dat$id,
    study_label = dat$study_label,
    y = dat$number_inf,
    m = dat$number_dissect,
    observed = dat$observed_proportion,
    predicted_median = apply(p_matrix, 2, stats::median),
    predicted_mean_lower = apply(
      p_matrix, 2, stats::quantile, probs = 0.025
    ),
    predicted_mean_upper = apply(
      p_matrix, 2, stats::quantile, probs = 0.975
    ),
    predictive_lower = apply(
      replicated_proportion, 2, stats::quantile, probs = 0.025
    ),
    predictive_upper = apply(
      replicated_proportion, 2, stats::quantile, probs = 0.975
    )
  ) |>
    dplyr::mutate(
      covered = observed >= predictive_lower & observed <= predictive_upper,
      predictive_width = predictive_upper - predictive_lower
    )

  # Panel A uses the model-estimated mean probability on x and the observed
  # proportion on y. Horizontal coloured segments show the 95% credible
  # interval for the mean probability; vertical grey segments show the 95%
  # posterior predictive interval for an observed feeding-assay proportion.
  panel_a <- ggplot2::ggplot(
    pred,
    ggplot2::aes(
      x = predicted_median,
      y = observed,
      colour = study_label
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = predicted_median,
        xend = predicted_median,
        y = predictive_lower,
        yend = predictive_upper
      ),
      inherit.aes = TRUE,
      colour = PREDICTIVE_LINE_COLOUR,
      linewidth = 0.28,
      alpha = 0.17
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = predicted_mean_lower,
        xend = predicted_mean_upper,
        y = observed,
        yend = observed
      ),
      linewidth = 0.35,
      alpha = 0.30
    ) +
    ggplot2::geom_point(size = 1.55, alpha = 0.72) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2,
      colour = "grey45",
      linewidth = 0.50
    ) +
    ggplot2::facet_wrap(~study_label, nrow = 1) +
    ggplot2::scale_colour_manual(values = STUDY_COLOURS, guide = "none") +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::labs(
      x = "Model-estimated mean infection probability (%)",
      y = "Observed proportion of mosquitoes infected (%)"
    ) +
    manuscript_theme(8)

  # Five groups within each study. Predicted means are weighted by the number
  # of mosquitoes dissected, matching the observed pooled mosquito proportion.
  calibration <- pred |>
    dplyr::group_by(study_label) |>
    dplyr::mutate(
      prediction_group = dplyr::ntile(predicted_median, 5L)
    ) |>
    dplyr::group_by(study_label, prediction_group) |>
    dplyr::summarise(
      predicted = stats::weighted.mean(predicted_median, w = m),
      observed = sum(y) / sum(m),
      participants = dplyr::n(),
      mosquitoes = sum(m),
      .groups = "drop"
    )

  panel_b <- ggplot2::ggplot(
    calibration,
    ggplot2::aes(
      x = predicted,
      y = observed,
      colour = study_label,
      size = participants
    )
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2,
      colour = "grey45",
      linewidth = 0.50
    ) +
    ggplot2::geom_point(alpha = 0.88) +
    ggplot2::facet_wrap(~study_label, nrow = 1) +
    ggplot2::scale_colour_manual(values = STUDY_COLOURS, guide = "none") +
    ggplot2::scale_size_continuous(range = c(2.4, 5.2)) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::labs(
      x = "Mean predicted infection probability (%)",
      y = "Observed proportion of mosquitoes infected (%)",
      size = "Participants"
    ) +
    manuscript_theme(8) +
    ggplot2::theme(legend.position = "top")

  coverage <- pred |>
    dplyr::group_by(study_label) |>
    dplyr::summarise(
      participants = dplyr::n(),
      predictive_coverage = mean(covered),
      median_predictive_width = stats::median(predictive_width),
      .groups = "drop"
    )

  fig6 <- panel_a / panel_b +
    patchwork::plot_annotation(tag_levels = "A")

  save_manuscript_plot(
    fig6,
    "Fig6_model_fit_calibration_predictive_uncertainty_FINAL",
    width = 8.4,
    height = 6.4
  )

  readr::write_csv(
    pred,
    file.path(TABLE_OUTPUT_DIR, "Fig6_participant_predictions.csv")
  )
  readr::write_csv(
    calibration,
    file.path(TABLE_OUTPUT_DIR, "Fig6_calibration_groups.csv")
  )
  readr::write_csv(
    coverage,
    file.path(TABLE_OUTPUT_DIR, "Fig6_predictive_coverage.csv")
  )

  invisible(list(
    plot = fig6,
    predictions = pred,
    calibration = calibration,
    coverage = coverage
  ))
}

# =============================================================================
# 8. GENERATE BOTH FINAL FIGURES
# =============================================================================
fig5_result <- make_final_fig5(
  fit = primary_fit,
  dat = primary_dat,
  max_draws = MAX_POSTERIOR_DRAWS
)

fig6_result <- make_final_fig6(
  fit = primary_fit,
  dat = primary_dat
)

