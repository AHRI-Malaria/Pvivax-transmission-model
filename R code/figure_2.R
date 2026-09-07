
# =============================================================================
#  FIGURE 2 
# =============================================================================

library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(readr)

get_project_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/")
    script_dir <- dirname(script_path)
    if (basename(script_dir) == "figures") return(dirname(script_dir))
    return(script_dir)
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
PROJECT_DIR <- get_project_dir()

# -----------------------------------------------------------------------------
# 1. Data location
# -----------------------------------------------------------------------------
DATA_FILE <- file.path(PROJECT_DIR, "data", "data_pvivax.csv")

if (!exists("primary_dat", inherits = TRUE)) {

  if (!file.exists(DATA_FILE)) {
    stop("Data file not found: ", DATA_FILE)
  }

  message(
    "Reading microscopy data from: ",
    normalizePath(DATA_FILE, winslash = "/")
  )

  raw_fig2 <- readr::read_csv(
    DATA_FILE,
    show_col_types = FALSE
  )

  required_raw_columns <- c(
    "id",
    "study",
    "par",
    "gam",
    "number_inf",
    "number_dissect"
  )

  missing_raw_columns <- setdiff(
    required_raw_columns,
    names(raw_fig2)
  )

  if (length(missing_raw_columns) > 0L) {
    stop(
      "Missing required columns in the microscopy CSV: ",
      paste(missing_raw_columns, collapse = ", ")
    )
  }

  primary_dat <- raw_fig2 |>
    dplyr::transmute(
      id = as.character(id),
      study_original = as.integer(study),
      par = as.numeric(par),
      gam = as.numeric(gam),
      number_inf = as.integer(number_inf),
      number_dissect = as.integer(number_dissect)
    ) |>
    dplyr::filter(
      study_original %in% c(1L, 2L, 3L)
    ) |>
    dplyr::distinct(
      study_original,
      id,
      .keep_all = TRUE
    ) |>
    dplyr::mutate(
      study_label = factor(
        paste0("Study ", study_original),
        levels = c("Study 1", "Study 2", "Study 3")
      )
    ) |>
    dplyr::arrange(
      study_original,
      id
    )
}

# -----------------------------------------------------------------------------
# 2. Validate plotting data
# -----------------------------------------------------------------------------

required_columns <- c(
  "study_label",
  "par",
  "gam",
  "number_inf",
  "number_dissect"
)

missing_columns <- setdiff(
  required_columns,
  names(primary_dat)
)

if (length(missing_columns) > 0L) {
  stop(
    "The plotting dataset is missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# 3.  colour palette
# -----------------------------------------------------------------------------

study_palette <- c(
  "Study 1" = "#5D88A8",
  "Study 2" = "#C88B56",
  "Study 3" = "#6F9B7A"
)

# -----------------------------------------------------------------------------
# 4.theme
# -----------------------------------------------------------------------------

manuscript_theme <- function(base_size = 9) {

  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      axis.title = ggplot2::element_text(
        size = base_size + 0.4,
        colour = "black"
      ),
      axis.text = ggplot2::element_text(
        size = base_size - 0.2,
        colour = "black"
      ),
      axis.text.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 4)
      ),
      axis.line = ggplot2::element_line(
        colour = "black",
        linewidth = 0.45
      ),
      axis.ticks = ggplot2::element_line(
        colour = "black",
        linewidth = 0.35
      ),
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = base_size + 2
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(7, 8, 7, 7)
    )
}

# -----------------------------------------------------------------------------
# 5. plot
# -----------------------------------------------------------------------------

make_gentle_violin_panel <- function(
    dat,
    y_variable,
    y_label,
    percent_axis = FALSE,
    y_lower = NA_real_) {

  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = study_label,
      y = .data[[y_variable]],
      colour = study_label
    )
  ) +
    ggplot2::geom_violin(
      fill = "white",
      trim = TRUE,
      scale = "width",
      width = 0.82,
      linewidth = 0.70,
      alpha = 1,
      na.rm = TRUE
    ) +
    ggplot2::geom_boxplot(
      fill = "white",
      width = 0.16,
      linewidth = 0.70,
      outlier.shape = NA,
      fatten = 1.8,
      alpha = 1,
      na.rm = TRUE
    ) +
    ggplot2::geom_jitter(
      shape = 16,
      width = 0.085,
      height = 0,
      size = 1.15,
      alpha = 0.42,
      na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(
      values = study_palette,
      drop = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = y_label
    ) +
    manuscript_theme(base_size = 9)

  if (percent_axis) {

    p <- p +
      ggplot2::scale_y_continuous(
        breaks = c(0, 0.25, 0.50, 0.75, 1),
        labels = scales::label_percent(accuracy = 1),
        expand = ggplot2::expansion(mult = c(0.01, 0.05))
      ) +
      ggplot2::coord_cartesian(
        ylim = c(0, 1),
        clip = "on"
      )

  } else {

    p <- p +
      ggplot2::scale_y_continuous(
        breaks = scales::breaks_pretty(n = 5),
        expand = ggplot2::expansion(mult = c(0.03, 0.07))
      )

    if (!is.na(y_lower)) {
      p <- p +
        ggplot2::coord_cartesian(
          ylim = c(y_lower, NA),
          clip = "on"
        )
    }
  }

  p
}

# -----------------------------------------------------------------------------
# 6. Prepare observed values
# -----------------------------------------------------------------------------

plot_dat <- primary_dat |>
  dplyr::mutate(
    study_label = factor(
      as.character(study_label),
      levels = c("Study 1", "Study 2", "Study 3")
    ),

    # Retain microscopy-zero observations for descriptive visualization.
    observed_logP = log10(par + 0.1),
    observed_logG = log10(gam + 0.1),

    observed_proportion =
      number_inf / number_dissect
  ) |>
  dplyr::filter(
    !is.na(study_label),
    number_dissect > 0,
    is.finite(observed_logP),
    is.finite(observed_logG),
    is.finite(observed_proportion)
  )

# -----------------------------------------------------------------------------
# 7. Create Figure 2 panels
# -----------------------------------------------------------------------------

# Panel A: all observed asexual parasite densities
panel_a <- make_gentle_violin_panel(
  dat = plot_dat,
  y_variable = "observed_logP",
  y_label = "Asexual parasite density (log10 scale)"
)

# Panel B: POSITIVE gametocyte densities only, for clearer visualization
panel_b_dat <- plot_dat |>
  dplyr::filter(gam > 0)

panel_b <- make_gentle_violin_panel(
  dat = panel_b_dat,
  y_variable = "observed_logG",
  y_label = "Gametocyte density (log10 scale)",
  y_lower = 0
)

# Panel C: all observed proportions infected
panel_c <- make_gentle_violin_panel(
  dat = plot_dat,
  y_variable = "observed_proportion",
  y_label = "Proportion of mosquitoes infected (%)",
  percent_axis = TRUE
)

# Panel D: all numbers dissected
panel_d <- make_gentle_violin_panel(
  dat = plot_dat,
  y_variable = "number_dissect",
  y_label = "Number of mosquitoes dissected"
)

# -----------------------------------------------------------------------------
# 8. Combine panels
# -----------------------------------------------------------------------------

fig2_final <- (
  panel_a + panel_b
) / (
  panel_c + panel_d
) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 11
    )
  )

# -----------------------------------------------------------------------------
# 9. Export
# -----------------------------------------------------------------------------

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "final_figures"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Fig2_observed_microscopy_and_feeding_data_FINAL.pdf"
  ),
  plot = fig2_final,
  width = 8.0,
  height = 7.0,
  units = "in",
  device = grDevices::cairo_pdf
)

ggplot2::ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Fig2_observed_microscopy_and_feeding_data_FINAL.tiff"
  ),
  plot = fig2_final,
  width = 8.0,
  height = 7.0,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

print(fig2_final)


