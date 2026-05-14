# build_sensitivity_figures.R
# Cross-reference comparison figures.
# Run from the bundle root.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# Locate bundle root: this script lives in scripts/, so go up one level
HERE <- dirname(normalizePath(sys.frame(1)$ofile %||% rstudioapi::getActiveDocumentContext()$path))
BUNDLE_ROOT <- normalizePath(file.path(HERE, ".."))
RESULTS_DIR <- file.path(BUNDLE_ROOT, "results")
OUTPUT_DIR  <- file.path(BUNDLE_ROOT, "figs_sensitivity")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

RUN_DIRS <- list(
  "Loans (main)" = file.path(RESULTS_DIR, "h8_weekly_seed_symmetric_diagnostic_loans_main"),
  "Cash"         = file.path(RESULTS_DIR, "20260508_142126_h8_sensitivity_ref_cash"),
  "Securities"   = file.path(RESULTS_DIR, "20260511_141952_h8_sensitivity_ref_securities"),
  "Other"        = file.path(RESULTS_DIR, "20260508_142710_h8_sensitivity_ref_other")
)

load_elpd_for_ref <- function(ref_name, run_dir) {
  path <- file.path(run_dir, "tables", "rolling_elpd_cov95.csv")
  df <- readr::read_csv(path, show_col_types = FALSE)
  df <- df %>% arrange(date) %>%
    mutate(reference = ref_name, cum_elpd_diff = cumsum(ELPD_diff))
  df
}

load_divs_for_ref <- function(ref_name, run_dir) {
  path <- file.path(run_dir, "tables", "rolling_diagnostic_summary.csv")
  df <- readr::read_csv(path, show_col_types = FALSE)
  df %>%
    mutate(reference = ref_name,
           spec = ifelse(grepl("Centered", model), "Centered", "Raw")) %>%
    select(reference, spec, total_divergences)
}

elpd_long <- bind_rows(lapply(names(RUN_DIRS), function(r) load_elpd_for_ref(r, RUN_DIRS[[r]])))
div_long  <- bind_rows(lapply(names(RUN_DIRS), function(r) load_divs_for_ref(r, RUN_DIRS[[r]])))
ref_levels <- names(RUN_DIRS)
elpd_long$reference <- factor(elpd_long$reference, levels = ref_levels)
div_long$reference  <- factor(div_long$reference,  levels = ref_levels)

# Figure 1: cumulative ELPD diff
p_cumelpd <- ggplot(elpd_long,
                    aes(x = as.Date(date), y = cum_elpd_diff,
                        color = reference, group = reference)) +
  geom_hline(yintercept = 0, color = "gray60", linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
  scale_color_manual(values = c(
    "Loans (main)" = "#1f77b4", "Cash" = "#d62728",
    "Securities"   = "#2ca02c", "Other" = "#9467bd"
  )) +
  labs(x = NULL,
       y = "Cumulative ELPD difference (Centered \u2212 Raw), nats",
       color = "ALR reference") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave(file.path(OUTPUT_DIR, "sensitivity_cumelpd_by_ref.png"),
       p_cumelpd, width = 8, height = 4.5, dpi = 200, bg = "white")

# Figure 2: total rolling divergences
p_divs <- ggplot(div_long,
                 aes(x = reference, y = total_divergences, fill = spec)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = total_divergences),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c("Centered" = "#1f77b4", "Raw" = "#d62728")) +
  labs(x = "ALR reference",
       y = "Total HMC divergent transitions (rolling fits)",
       fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  expand_limits(y = max(div_long$total_divergences) * 1.12)

ggsave(file.path(OUTPUT_DIR, "sensitivity_divergences_by_ref.png"),
       p_divs, width = 7, height = 4.5, dpi = 200, bg = "white")

cat("Wrote:\n",
    file.path(OUTPUT_DIR, "sensitivity_cumelpd_by_ref.png"), "\n",
    file.path(OUTPUT_DIR, "sensitivity_divergences_by_ref.png"), "\n",
    sep = "")
