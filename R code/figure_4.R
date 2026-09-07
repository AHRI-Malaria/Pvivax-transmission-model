# Final Figure 4: posterior conditional and pathway-specific associations.
# All plotted estimates are calculated from the saved primary posterior draws.

required_packages <- c("cmdstanr", "posterior", "dplyr", "tibble", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install these packages first: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(cmdstanr); library(posterior); library(dplyr); library(tibble)
  library(ggplot2); library(patchwork)
})

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

FIT_RDS <- file.path(PROJECT_DIR, "results", "fits", "primary", "primary.rds")
OUT_DIR <- file.path(PROJECT_DIR, "results", "final_figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(FIT_RDS)) stop("Primary fit not found: ", FIT_RDS, "\nRun 01_run_analysis.R first.")
fit <- readRDS(FIT_RDS)
AGE_REPORT_MULTIPLIER <- 15 / 10

scalar_draw <- function(variable) {
  x <- as.matrix(fit$draws(variable, format = "matrix"))
  if (ncol(x) != 1L) stop(variable, " is not scalar.")
  as.numeric(x[, 1L])
}

summarise_or <- function(variable, label, multiplier = 1) {
  x <- exp(multiplier * scalar_draw(variable))
  tibble(
    term = label,
    OR = median(x),
    lower = quantile(x, 0.025, names = FALSE),
    upper = quantile(x, 0.975, names = FALSE)
  )
}

panel_A <- bind_rows(
  summarise_or("beta_age", "Host age", AGE_REPORT_MULTIPLIER),
  summarise_or("beta_fever", "Fever\n(yes vs no)"),
  summarise_or("beta_G", "Gametocyte density"),
  summarise_or("beta_P", "Asexual parasite density\n(residual conditional)")
) |>
  mutate(term = factor(term, levels = rev(term)))

panel_B <- bind_rows(
  summarise_or("assoc_age_total", "Age: total model-implied\nassociation", AGE_REPORT_MULTIPLIER),
  summarise_or("assoc_P_via_G", "Asexual parasite density:\npathway via gametocytes"),
  summarise_or("assoc_P_total", "Asexual parasite density:\ntotal model-implied\nassociation")
) |>
  mutate(term = factor(term, levels = rev(term)))

forest_theme <- theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "sans", colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 11, margin = margin(t = 7)),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 10.2, colour = "black", lineheight = 0.95, margin = margin(r = 6)),
    axis.line = element_line(colour = "black", linewidth = 0.55),
    axis.ticks = element_line(colour = "black", linewidth = 0.45),
    panel.grid = element_blank(),
    plot.margin = margin(t = 8, r = 8, b = 6, l = 5)
  )

make_panel <- function(dat, xlab) {
  ggplot(dat, aes(x = OR, y = term)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.55) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.11,
                  linewidth = 0.8, colour = "#6FA3CF") +
    geom_point(size = 3.0, shape = 16, colour = "#1558A0") +
    scale_x_log10(breaks = c(0.5, 1, 2, 4), labels = c("0.5", "1.0", "2.0", "4.0"),
                  limits = c(0.40, 5.0), expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = xlab) + forest_theme
}

final_figure <- make_panel(panel_A, "Odds ratio (95% CrI)") +
  make_panel(panel_B, "Model-implied odds ratio (95% CrI)") +
  plot_layout(widths = c(1, 1.03)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 14)))

write.csv(panel_A, file.path(OUT_DIR, "Fig4_panelA_values.csv"), row.names = FALSE)
write.csv(panel_B, file.path(OUT_DIR, "Fig4_panelB_values.csv"), row.names = FALSE)

ggsave(file.path(OUT_DIR, "Fig4_posterior_associations_FINAL.pdf"), final_figure,
       width = 9, height = 5.2, units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(OUT_DIR, "Fig4_posterior_associations_FINAL.tiff"), final_figure,
       width = 9, height = 5.2, units = "in", dpi = 600, device = "tiff",
       compression = "lzw", bg = "white")

print(final_figure)
