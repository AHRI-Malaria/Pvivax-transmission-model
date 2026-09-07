
# =============================================================================
#  FIGURE 3 

# =============================================================================

library(dplyr)
library(ggplot2)
library(scales)

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


# -----------------------------------------------------------------------------
# 0. Load analysis data 
# -----------------------------------------------------------------------------
if (!exists("primary_dat", inherits = TRUE)) {
  data_file <- file.path(PROJECT_DIR, "data", "data_pvivax.csv")
  if (!file.exists(data_file)) stop("Data file not found: ", data_file)
  raw_fig3 <- read.csv(data_file, stringsAsFactors = FALSE, check.names = FALSE)
  required_raw <- c("id", "study", "par", "gam", "number_inf", "number_dissect")
  missing_raw <- setdiff(required_raw, names(raw_fig3))
  if (length(missing_raw)) stop("Missing required data columns: ", paste(missing_raw, collapse = ", "))
  primary_dat <- raw_fig3 |>
    dplyr::transmute(
      id = as.character(id),
      study_original = as.integer(study),
      par = as.numeric(par),
      gam = as.numeric(gam),
      number_inf = as.integer(number_inf),
      number_dissect = as.integer(number_dissect)
    ) |>
    dplyr::filter(study_original %in% c(1L, 2L, 3L)) |>
    dplyr::distinct(study_original, id, .keep_all = TRUE) |>
    dplyr::mutate(
      study_label = factor(paste0("Study ", study_original), levels = c("Study 1", "Study 2", "Study 3"))
    ) |>
    dplyr::arrange(study_original, id)
}

# -----------------------------------------------------------------------------
# 1.  theme
# -----------------------------------------------------------------------------

manuscript_theme <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      strip.background = element_rect(
        fill = "grey96",
        colour = "grey80",
        linewidth = 0.4
      ),
      strip.text = element_text(face = "bold", size = base_size),
      axis.title = element_text(size = base_size + 0.5),
      axis.text = element_text(size = base_size - 0.5),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.5),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = "white", colour = NA),
      plot.margin = margin(7, 9, 7, 7)
    )
}

# -----------------------------------------------------------------------------
# 2. Format Spearman 
# -----------------------------------------------------------------------------

format_spearman_label <- function(rho, p_value) {
  p_text <- if (is.na(p_value)) {
    "P = NA"
  } else if (p_value < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", formatC(p_value, format = "f", digits = 3))
  }

  paste0(
    " \u03C1 = ",
    formatC(rho, format = "f", digits = 2),
    ", ",
    p_text
  )
}

# -----------------------------------------------------------------------------
# 3. Create Figure 3
# -----------------------------------------------------------------------------

make_figure3_observed_density_relationship <- function(
    dat,
    output_dir = file.path(PROJECT_DIR, "results", "final_figures"),
    file_stem = "Fig3_observed_density_relationship_FINAL") {

  required_columns <- c(
    "study_label",
    "par",
    "gam",
    "number_inf",
    "number_dissect"
  )

  missing_columns <- setdiff(required_columns, names(dat))

  if (length(missing_columns) > 0L) {
    stop(
      "The following required columns are missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------------------
  # Observed plotting data
  # ---------------------------------------------------------------------------

  plot_dat <- dat |>
    mutate(
      study_label = factor(
        as.character(study_label),
        levels = c("Study 1", "Study 2", "Study 3")
      ),
      observed_proportion = number_inf / number_dissect,

      # Descriptive transformation only; microscopy zeros are retained.
      observed_logP = log10(par + 0.1),
      observed_logG = log10(gam + 0.1)
    ) |>
    filter(
      !is.na(study_label),
      is.finite(observed_logP),
      is.finite(observed_logG),
      is.finite(observed_proportion),
      number_dissect > 0
    )

  # ---------------------------------------------------------------------------
  # Study-specific Spearman correlations
  # ---------------------------------------------------------------------------

  correlation_table <- plot_dat |>
    group_by(study_label) |>
    group_modify(
      ~ {
        test <- suppressWarnings(
          cor.test(
            x = .x$par,
            y = .x$gam,
            method = "spearman",
            exact = FALSE
          )
        )

        tibble(
          n = sum(complete.cases(.x$par, .x$gam)),
          rho = unname(test$estimate),
          p_value = test$p.value
        )
      }
    ) |>
    ungroup() |>
    mutate(
      annotation = mapply(
        format_spearman_label,
        rho,
        p_value,
        USE.NAMES = FALSE
      ),

      # -Inf/Inf places the label consistently in the upper-left corner
      # of every facet, independent of panel-specific data ranges.
      x_position = -Inf,
      y_position = Inf
    )

  # Export exact correlation results used in the figure.
  write.csv(
    correlation_table |>
      select(study_label, n, rho, p_value),
    file = file.path(
      output_dir,
      paste0(file_stem, "_Spearman_results.csv")
    ),
    row.names = FALSE
  )

  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------

  p <- ggplot(
    plot_dat,
    aes(
      x = observed_logP,
      y = observed_logG,
      colour = observed_proportion
    )
  ) +
    geom_point(
      alpha = 0.78,
      size = 1.55
    ) +

    # Correlation annotation in the upper-left corner.
    geom_label(
      data = correlation_table,
      aes(
        x = x_position,
        y = y_position,
        label = annotation
      ),
      inherit.aes = FALSE,
      hjust = -0.08,
      vjust = 1.12,
      size = 3.15,
      lineheight = 1.05,
      colour = "black",
      fill = alpha("white", 0.85),
      label.size = 0,
      label.padding = unit(0.12, "lines")
    ) +

    facet_wrap(
      ~study_label,
      nrow = 1,
      drop = FALSE
    ) +

    scale_colour_gradientn(
      colours = c(
        "#4F81A8",
        "#A7C4D7",
        "#E8DCA8",
        "#D9A17E",
        "#C98273"
      ),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 1),
      labels = label_percent(accuracy = 1),
      name = paste0(
        "Observed proportion\n",
        "of mosquitoes infected"
      )
    ) +

    labs(
      x = "Asexual parasite density (log10 scale)",
      y = "Gametocyte density (log10 scale)"
    ) +

    manuscript_theme(base_size = 9) +

    theme(
      legend.position = "right"
    )

  # ---------------------------------------------------------------------------
  # Export
  # ---------------------------------------------------------------------------

  ggsave(
    filename = file.path(output_dir, paste0(file_stem, ".pdf")),
    plot = p,
    width = 8.2,
    height = 4.5,
    units = "in",
    device = cairo_pdf
  )

  ggsave(
    filename = file.path(output_dir, paste0(file_stem, ".tiff")),
    plot = p,
    width = 8.2,
    height = 4.5,
    units = "in",
    dpi = 600,
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )

  invisible(
    list(
      plot = p,
      spearman_results = correlation_table |>
        select(study_label, n, rho, p_value)
    )
  )
}

# -----------------------------------------------------------------------------
# 4. Run Figure 3
# -----------------------------------------------------------------------------

fig3_result <- make_figure3_observed_density_relationship(
  dat = primary_dat
)

print(fig3_result$plot)
print(fig3_result$spearman_results)

