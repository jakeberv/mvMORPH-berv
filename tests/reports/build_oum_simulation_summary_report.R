#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[1]))
} else {
  normalizePath("tests/reports/build_oum_simulation_summary_report.R")
}

repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
report_dir <- file.path(repo_root, "tests", "reports")
figure_dir <- file.path(report_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

theme_report <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      legend.position = "bottom"
    )
}

save_plot <- function(plot, filename, width = 9, height = 5.5) {
  path <- file.path(figure_dir, filename)
  ggsave(path, plot = plot, width = width, height = height, dpi = 160, bg = "white")
  path
}

theta_summary_path <- file.path(
  repo_root,
  "archon-pulls", "tg18", "tests", "experimental_oum_theta_recovery_summary.csv"
)
theta_boundary_path <- file.path(
  repo_root,
  "archon-pulls", "tg18", "tests", "experimental_oum_theta_recovery_boundary.csv"
)

if (!file.exists(theta_summary_path) || !file.exists(theta_boundary_path)) {
  stop("Theta recovery CSVs were not found under archon-pulls/tg18/tests.", call. = FALSE)
}

theta_summary <- read.csv(theta_summary_path, stringsAsFactors = FALSE)
theta_boundary <- read.csv(theta_boundary_path, stringsAsFactors = FALSE)

null_df <- theta_summary %>% filter(phase == "focused_null")
recovery_df <- theta_summary %>% filter(phase == "theta_recovery")

historical_covariate_null <- data.frame(
  predictor = rep(c("cont_indep", "cont_confounded", "disc_binary", "disc3_sparse_confounded"), each = 3),
  method = rep(c("LL", "EmpBayes", "LOOCV"), times = 4),
  null_selection_rate = c(
    0.138, 0.200, 0.238,
    0.125, 0.200, 0.250,
    0.150, 0.350, 0.338,
    0.125, 0.275, 0.263
  )
)

historical_bmm_oum_head_to_head <- data.frame(
  truth = rep(c("pure_BMM", "pure_OUM"), each = 3),
  method = rep(c("LL", "EmpBayes", "LOOCV"), times = 2),
  correct_model_rate = c(
    0.961, 0.904, 0.879,
    0.712, 0.813, 0.839
  )
)

historical_four_model_shared <- data.frame(
  truth = c("shared_BM", "shared_OU"),
  correct_model_rate = c(0.466, 0.824)
)

historical_detectability_notes <- list(
  null_false_oum_balanced = 0.047,
  null_false_oum_imbalanced = 0.092,
  frontier_ge50 = 182,
  frontier_total = 2400,
  frontier_ge80 = 42,
  misspec_ge50 = 29,
  misspec_total = 4032,
  misspec_ge80 = 1
)

plot_covariate_null <- ggplot(
  historical_covariate_null,
  aes(x = predictor, y = null_selection_rate, fill = method)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 0.4)) +
  scale_fill_manual(values = c("LL" = "#2F6C8F", "EmpBayes" = "#D97A32", "LOOCV" = "#7B8D42")) +
  labs(
    title = "Near-Null Covariate Selection Inside OUM",
    subtitle = "Dense covariate suite; lower is more conservative",
    x = NULL,
    y = "Full-model selection rate under null",
    fill = "Method"
  ) +
  theme_report() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

plot_bmm_oum_head_to_head <- ggplot(
  historical_bmm_oum_head_to_head,
  aes(x = method, y = correct_model_rate, fill = truth)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = c("pure_BMM" = "#4A7C59", "pure_OUM" = "#8D5A97")) +
  labs(
    title = "BMM vs OUM Head-to-Head Recovery",
    subtitle = "Direct two-model comparison from the BMM/OUM simulation suite",
    x = NULL,
    y = "Correct-model selection rate",
    fill = "Truth"
  ) +
  theme_report()

theta_null_heat <- null_df %>%
  group_by(n_tips, p, alpha) %>%
  summarise(false_oum_rate = mean(select_oum_rate), .groups = "drop")

plot_theta_null_heat <- ggplot(
  theta_null_heat,
  aes(x = factor(p), y = factor(n_tips), fill = false_oum_rate)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f%%", 100 * false_oum_rate)), size = 3.2) +
  facet_wrap(~ alpha, nrow = 1, labeller = label_bquote(alpha == .(alpha))) +
  scale_fill_gradient(low = "#EEF5FB", high = "#B6413B", labels = label_percent(accuracy = 1)) +
  labs(
    title = "False OUM Selection in the Focused Null Study",
    subtitle = "The small-n, high-p corner remains the main failure mode",
    x = "Trait dimension (p)",
    y = "Number of taxa",
    fill = "False OUM"
  ) +
  theme_report()

theta_select_signal <- recovery_df %>%
  group_by(n_tips, alpha, theta_sep) %>%
  summarise(select_oum_rate = mean(select_oum_rate), .groups = "drop")

plot_theta_select_signal <- ggplot(
  theta_select_signal,
  aes(x = theta_sep, y = select_oum_rate, color = factor(alpha))
) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  facet_wrap(~ n_tips, nrow = 1) +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1)) +
  scale_color_manual(values = c("0.5" = "#5C7AEA", "1" = "#2E8B57", "1.5" = "#C75B39")) +
  labs(
    title = "Focused OUM Detectability Improves with Signal, Alpha, and Taxa",
    subtitle = "Average OUM selection rate across predictors, covariate settings, and p",
    x = expression(theta[sep]),
    y = "OUM selected",
    color = expression(alpha)
  ) +
  theme_report()

theta_recovery_signal <- recovery_df %>%
  group_by(alpha, theta_sep) %>%
  summarise(
    theta_cor = mean(theta_cor_oum_mean),
    theta_spread_error = mean(theta_spread_abs_error_oum_mean),
    .groups = "drop"
  )

theta_recovery_long <- rbind(
  data.frame(
    alpha = theta_recovery_signal$alpha,
    theta_sep = theta_recovery_signal$theta_sep,
    metric = "Theta correlation",
    value = theta_recovery_signal$theta_cor
  ),
  data.frame(
    alpha = theta_recovery_signal$alpha,
    theta_sep = theta_recovery_signal$theta_sep,
    metric = "Theta spread abs. error",
    value = theta_recovery_signal$theta_spread_error
  )
)

plot_theta_recovery_signal <- ggplot(
  theta_recovery_long,
  aes(x = theta_sep, y = value, color = factor(alpha))
) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c("0.5" = "#5C7AEA", "1" = "#2E8B57", "1.5" = "#C75B39")) +
  labs(
    title = "Theta Recovery Improves Even Before Selection Becomes Certain",
    subtitle = "Higher alpha and stronger theta separation improve recovered regime optima",
    x = expression(theta[sep]),
    y = NULL,
    color = expression(alpha)
  ) +
  theme_report()

theta_design <- recovery_df %>%
  group_by(n_tips, p) %>%
  summarise(
    select_oum_rate = mean(select_oum_rate),
    theta_cor = mean(theta_cor_oum_mean),
    .groups = "drop"
  )

theta_design_long <- rbind(
  data.frame(
    n_tips = theta_design$n_tips,
    p = theta_design$p,
    metric = "OUM selected",
    value = theta_design$select_oum_rate
  ),
  data.frame(
    n_tips = theta_design$n_tips,
    p = theta_design$p,
    metric = "Theta correlation",
    value = theta_design$theta_cor
  )
)

plot_theta_design <- ggplot(
  theta_design_long,
  aes(x = factor(p), y = factor(n_tips), fill = value)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f%%", 100 * value)), size = 3.2) +
  facet_wrap(~ metric, nrow = 1) +
  scale_fill_gradient(low = "#EEF5FB", high = "#2F6C8F", labels = label_percent(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Selection Strength and Theta Recovery Diverge in the Hardest Corner",
    subtitle = "Averaged across theta separation, alpha, predictors, and covariate settings",
    x = "Trait dimension (p)",
    y = "Number of taxa",
    fill = "Average"
  ) +
  theme_report()

fig_covariate <- save_plot(plot_covariate_null, "covariate_null_selection_rates.png", width = 9, height = 5.5)
fig_head_to_head <- save_plot(plot_bmm_oum_head_to_head, "bmm_oum_head_to_head.png", width = 8.5, height = 5.2)
fig_theta_null <- save_plot(plot_theta_null_heat, "theta_null_heatmap.png", width = 10.5, height = 5.4)
fig_theta_select <- save_plot(plot_theta_select_signal, "theta_selection_vs_signal.png", width = 10.5, height = 5.4)
fig_theta_recovery <- save_plot(plot_theta_recovery_signal, "theta_recovery_vs_signal.png", width = 10.5, height = 5.4)
fig_theta_design <- save_plot(plot_theta_design, "theta_recovery_design_heatmap.png", width = 10.5, height = 5.4)

null_overall <- mean(null_df$select_oum_rate)
null_by_n <- null_df %>% group_by(n_tips) %>% summarise(rate = mean(select_oum_rate), .groups = "drop")
null_by_p <- null_df %>% group_by(p) %>% summarise(rate = mean(select_oum_rate), .groups = "drop")
recovery_by_theta <- recovery_df %>% group_by(theta_sep) %>% summarise(
  select_oum_rate = mean(select_oum_rate),
  theta_cor = mean(theta_cor_oum_mean),
  theta_spread_error = mean(theta_spread_abs_error_oum_mean),
  .groups = "drop"
)
recovery_by_alpha <- recovery_df %>% group_by(alpha) %>% summarise(
  select_oum_rate = mean(select_oum_rate),
  theta_rmse = mean(theta_rmse_oum_mean),
  theta_cor = mean(theta_cor_oum_mean),
  theta_spread_error = mean(theta_spread_abs_error_oum_mean),
  alpha_hat_ou = mean(alpha_hat_ou_mean),
  alpha_hat_oum = mean(alpha_hat_oum_mean),
  .groups = "drop"
)
recovery_by_n <- recovery_df %>% group_by(n_tips) %>% summarise(
  select_oum_rate = mean(select_oum_rate),
  theta_rmse = mean(theta_rmse_oum_mean),
  theta_cor = mean(theta_cor_oum_mean),
  theta_spread_error = mean(theta_spread_abs_error_oum_mean),
  .groups = "drop"
)
recovery_by_p <- recovery_df %>% group_by(p) %>% summarise(
  select_oum_rate = mean(select_oum_rate),
  theta_rmse = mean(theta_rmse_oum_mean),
  theta_cor = mean(theta_cor_oum_mean),
  theta_spread_error = mean(theta_spread_abs_error_oum_mean),
  .groups = "drop"
)

boundary_text <- theta_boundary %>%
  group_by(n_tips, p, alpha) %>%
  summarise(
    empirical50 = mean(empirical_signal50, na.rm = TRUE),
    empirical80 = mean(empirical_signal80, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(alpha, n_tips, p)

boundary_good_p8 <- boundary_text %>%
  filter(n_tips == 240, p == 8, alpha == 1.5)

boundary_good_p24 <- boundary_text %>%
  filter(n_tips == 240, p == 24, alpha == 1.5)

fmt_pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), 100 * x)
fmt_num <- function(x, digits = 3) sprintf(paste0("%.", digits, "f"), x)

md_path <- file.path(report_dir, "oum_simulation_summary_report.md")
html_path <- file.path(report_dir, "oum_simulation_summary_report.html")
css_path <- file.path(report_dir, "report.css")

report_lines <- c(
  "# OUM Simulation Summary Report",
  "",
  paste("Generated on", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "",
  "## Provenance",
  "",
  "- Figures and summaries for the focused theta-recovery study are generated directly from local pulled CSVs in `archon-pulls/tg18/tests/`.",
  "- Earlier Archon campaigns reused the same remote workspace, so some of their raw CSV outputs were overwritten later.",
  "- For those earlier campaigns, the report uses the quantitative results recorded during the session and labels them as reconstructed summaries.",
  "",
  "## Executive Summary",
  "",
  "- The `mvgls(model = \"OUM\")` covariate path is mechanically solid: it fit successfully, preserved covariates through prediction and downstream methods, and behaved sensibly on real data and targeted tests.",
  "- Within OUM, covariate detection is strong once signal is present, but near-null model selection is softer for discrete predictors and for more permissive selection methods (`LOOCV`, `EmpBayes`) than for `LL`.",
  "- Head-to-head `BMM` vs `OUM` comparisons under `GIC` are asymmetric: `BMM` is recovered strongly, while `OUM` is much easier to miss, especially when candidate sets include shared models.",
  "- Direct `OU` vs `OUM` comparison is better calibrated than `BMM` vs `OUM`, but still weak in mixed or high-dimensional small-sample settings.",
  "- The strongest new result is that `theta` recovery is materially better than model selection alone suggests. In favorable settings, regime optima can be estimated well even before `OUM` wins every model comparison.",
  "",
  "## 1. OUM With Covariates: What Worked",
  "",
  "- The implementation now supports standard covariate adjustment inside OUM fits, including painted `simmap` trees and ordinary regression covariates such as `log(mass)`.",
  "- Targeted tests showed that the updated OUM path works with `LL`, penalized fits, `EmpBayes`, `error=TRUE`, `predict()`, `ancestral()`, `manova.gls()`, `GIC()`, and `EIC()`.",
  "- On a real `phyllostomid` example, OUM plus a size covariate fit cleanly and returned both regime-optimum rows and regression-effect rows in the coefficient matrix.",
  "",
  "## 2. Covariate Selection Inside OUM (Reconstructed From Dense Covariate Suite)",
  "",
  paste0(
    "- In the near-null region, false full-model selection ranged from about ",
    fmt_pct(min(historical_covariate_null$null_selection_rate)), " to ",
    fmt_pct(max(historical_covariate_null$null_selection_rate)),
    ", with `LL` consistently more conservative than `EmpBayes` and `LOOCV`."
  ),
  "- Weak, medium, and strong signal settings selected the full covariate model essentially perfectly across methods in that suite.",
  "",
  paste0("![Near-null covariate selection rates](", basename(dirname(fig_covariate)), "/", basename(fig_covariate), ")"),
  "",
  "## 3. BMM vs OUM Head-to-Head (Reconstructed From Direct Two-Model Suite)",
  "",
  "- In the direct `BMM` vs `OUM` comparison, recovery was asymmetric.",
  paste0(
    "- Under pure `BMM` truth, the correct model was selected ",
    paste(sprintf("%s: %s", historical_bmm_oum_head_to_head$method[historical_bmm_oum_head_to_head$truth == "pure_BMM"], fmt_pct(historical_bmm_oum_head_to_head$correct_model_rate[historical_bmm_oum_head_to_head$truth == "pure_BMM"])), collapse = ", "),
    "."
  ),
  paste0(
    "- Under pure `OUM` truth, the correct model was selected ",
    paste(sprintf("%s: %s", historical_bmm_oum_head_to_head$method[historical_bmm_oum_head_to_head$truth == "pure_OUM"], fmt_pct(historical_bmm_oum_head_to_head$correct_model_rate[historical_bmm_oum_head_to_head$truth == "pure_OUM"])), collapse = ", "),
    "."
  ),
  "- In the same study, the shared null was still heavily tilted toward `BMM`, which was the first strong warning that raw `GIC` comparisons were not neutral between these model families.",
  "",
  paste0("![BMM vs OUM head-to-head recovery](", basename(dirname(fig_head_to_head)), "/", basename(fig_head_to_head), ")"),
  "",
  "## 4. Four-Model Calibration (Reconstructed From BM/OU/BMM/OUM Suite)",
  "",
  paste0(
    "- When the candidate set included `BM`, `OU`, `BMM`, and `OUM`, the shared-model truths were recovered asymmetrically: the correct model won ",
    fmt_pct(historical_four_model_shared$correct_model_rate[historical_four_model_shared$truth == "shared_BM"]),
    " under shared BM truth and ",
    fmt_pct(historical_four_model_shared$correct_model_rate[historical_four_model_shared$truth == "shared_OU"]),
    " under shared OU truth."
  ),
  "- In that calibration run, `OUM` was almost never the winner in the OUM frontier or mixed-signal phases. This showed that lack of an `OUM` win is weak evidence against optimum shifts when the candidate set is broad.",
  "",
  "## 5. Direct OU vs OUM Detectability (Reconstructed From Dedicated Detectability Suite)",
  "",
  paste0(
    "- Direct `OU` vs `OUM` comparison was better calibrated than `BMM` vs `OUM`: false `OUM` selection under shared-OU truth averaged ",
    fmt_pct(historical_detectability_notes$null_false_oum_balanced),
    " for balanced painted regimes and ",
    fmt_pct(historical_detectability_notes$null_false_oum_imbalanced),
    " for imbalanced regimes."
  ),
  paste0(
    "- But even there, detectability was limited: only ",
    historical_detectability_notes$frontier_ge50, " of ", historical_detectability_notes$frontier_total,
    " pure-OUM frontier cells crossed 50% `OUM` selection, and only ",
    historical_detectability_notes$frontier_ge80, " crossed 80%."
  ),
  paste0(
    "- Under BM-rate misspecification, detectability collapsed further: only ",
    historical_detectability_notes$misspec_ge50, " of ", historical_detectability_notes$misspec_total,
    " misspecified cells crossed 50% selection, and only ",
    historical_detectability_notes$misspec_ge80, " crossed 80%."
  ),
  "- That suite also showed alpha collapse under misspecification: both `OU` and `OUM` estimated low alpha values regardless of the true simulated alpha.",
  "",
  "## 6. Focused Theta-Recovery Study (Raw CSVs Pulled Locally)",
  "",
  paste0(
    "- This focused study used the regime where OUM looked most promising: 2 imbalanced painted regimes, higher alpha, stronger theta separation, and varying `n` and `p`."
  ),
  paste0(
    "- Overall false `OUM` selection under the null was ", fmt_pct(null_overall),
    ", but the failure mode was very concentrated."
  ),
  paste0(
    "- By number of taxa, false `OUM` selection fell from ",
    fmt_pct(null_by_n$rate[null_by_n$n_tips == 60]), " at `n=60` to ",
    fmt_pct(null_by_n$rate[null_by_n$n_tips == 120]), " at `n=120` and ",
    fmt_pct(null_by_n$rate[null_by_n$n_tips == 240]), " at `n=240`."
  ),
  paste0(
    "- By trait dimension, the main pathological corner was `p=48`, where the average false `OUM` rate was ",
    fmt_pct(null_by_p$rate[null_by_p$p == 48]), "."
  ),
  "",
  paste0("![False OUM selection heatmap](", basename(dirname(fig_theta_null)), "/", basename(fig_theta_null), ")"),
  "",
  paste0(
    "- Under pure simulated OUM, detectability improved strongly with theta separation: `OUM` was selected ",
    fmt_pct(recovery_by_theta$select_oum_rate[recovery_by_theta$theta_sep == 0.08]),
    " at `theta_sep = 0.08`, ",
    fmt_pct(recovery_by_theta$select_oum_rate[recovery_by_theta$theta_sep == 0.18]),
    " at `0.18`, ",
    fmt_pct(recovery_by_theta$select_oum_rate[recovery_by_theta$theta_sep == 0.30]),
    " at `0.30`, and ",
    fmt_pct(recovery_by_theta$select_oum_rate[recovery_by_theta$theta_sep == 0.60]),
    " at `0.60`."
  ),
  paste0(
    "- Higher alpha helped both selection and recovery. Average `OUM` selection rose from ",
    fmt_pct(recovery_by_alpha$select_oum_rate[recovery_by_alpha$alpha == 0.5]),
    " at `alpha = 0.5` to ",
    fmt_pct(recovery_by_alpha$select_oum_rate[recovery_by_alpha$alpha == 1.5]),
    " at `alpha = 1.5`."
  ),
  paste0(
    "- The fitted alpha values tracked the truth well in this focused regime: mean `alpha_hat` for OUM was ",
    paste(sprintf("%s -> %s", recovery_by_alpha$alpha, fmt_num(recovery_by_alpha$alpha_hat_oum, 2)), collapse = "; "),
    "."
  ),
  "",
  paste0("![OUM selection vs signal](", basename(dirname(fig_theta_select)), "/", basename(fig_theta_select), ")"),
  "",
  paste0(
    "- Theta recovery improved even before selection became certain. Average theta correlation rose from ",
    fmt_num(recovery_by_theta$theta_cor[recovery_by_theta$theta_sep == 0.08], 2),
    " at `theta_sep = 0.08` to ",
    fmt_num(recovery_by_theta$theta_cor[recovery_by_theta$theta_sep == 0.60], 2),
    " at `theta_sep = 0.60`."
  ),
  paste0(
    "- At the same time, theta spread absolute error fell from ",
    fmt_num(recovery_by_theta$theta_spread_error[recovery_by_theta$theta_sep == 0.08], 3),
    " to ",
    fmt_num(recovery_by_theta$theta_spread_error[recovery_by_theta$theta_sep == 0.60], 3),
    "."
  ),
  paste0(
    "- More taxa were consistently helpful: theta RMSE was ",
    paste(sprintf("n=%d -> %s", recovery_by_n$n_tips, fmt_num(recovery_by_n$theta_rmse, 3)), collapse = "; "),
    "."
  ),
  paste0(
    "- Trait dimension changed the quality of recovery more than the raw RMSE. Theta correlation averaged ",
    paste(sprintf("p=%d -> %s", recovery_by_p$p, fmt_num(recovery_by_p$theta_cor, 2)), collapse = "; "),
    ", showing that the high-dimensional corner stays hard even when absolute RMSE looks similar."
  ),
  paste0(
    "- In a favorable design (`n=240`, `p=8`, `alpha=1.5`), the empirical 50% and 80% selection thresholds were about ",
    fmt_num(boundary_good_p8$empirical50, 3), " and ", fmt_num(boundary_good_p8$empirical80, 3),
    "; at `p=24` with the same `n` and `alpha`, they rose to about ",
    fmt_num(boundary_good_p24$empirical50, 3), " and ", fmt_num(boundary_good_p24$empirical80, 3), "."
  ),
  "- The `n=60, p=48` corner remained pathological: it often selected `OUM` strongly under the null while yielding poor theta correlation, so those cells should not be interpreted as genuine detectability wins.",
  "",
  paste0("![Theta recovery vs signal](", basename(dirname(fig_theta_recovery)), "/", basename(fig_theta_recovery), ")"),
  "",
  paste0("![Selection versus recovery by design](", basename(dirname(fig_theta_design)), "/", basename(fig_theta_design), ")"),
  "",
  "## Practical Takeaways",
  "",
  "- The new OUM regression path is software-correct enough to use for size-adjusted multivariate analyses with painted regimes.",
  "- For model selection, `LL` is the safest default. `LOOCV` and `EmpBayes` are more permissive near the null.",
  "- Direct `OU` vs `OUM` comparison is more interpretable than `BMM` vs `OUM` when the scientific question is about optimum shifts.",
  "- Failure to select `OUM` is weak evidence against optimum shifts unless the dataset is in a favorable regime.",
  "- The favorable regime looks like: moderate-to-high alpha, sufficiently separated regime optima, at least about 120 taxa, and preferably not extreme small-n/high-p settings.",
  "- If the primary goal is recovering regime optima rather than just model selection, this focused study suggests theta can be estimated usefully even when `OUM` is not yet winning every comparison."
)

writeLines(report_lines, md_path)
writeLines(
  c(
    "body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; color: #1f2933; max-width: 980px; margin: 2rem auto 4rem; padding: 0 1.2rem; line-height: 1.65; background: #fbfbf8; }",
    "h1, h2 { color: #13293d; line-height: 1.25; }",
    "h1 { border-bottom: 3px solid #d7e3ec; padding-bottom: 0.35rem; }",
    "h2 { margin-top: 2.2rem; }",
    "img { display: block; margin: 1.2rem auto 2rem; max-width: 100%; height: auto; box-shadow: 0 8px 24px rgba(19, 41, 61, 0.12); border-radius: 8px; background: white; }",
    "code { background: #eef3f7; padding: 0.1rem 0.3rem; border-radius: 4px; }",
    "ul li { margin: 0.4rem 0; }",
    "p, li { font-size: 1rem; }"
  ),
  css_path
)

pandoc <- Sys.which("pandoc")
if (nzchar(pandoc)) {
  old_wd <- setwd(report_dir)
  on.exit(setwd(old_wd), add = TRUE)
  pandoc_cmd <- paste(
    shQuote(pandoc),
    shQuote(basename(md_path)),
    "--standalone",
    "--from", "gfm",
    "--to", "html5",
    "--css", shQuote(basename(css_path)),
    "--metadata", shQuote("title=OUM Simulation Summary Report"),
    "--output", shQuote(basename(html_path))
  )
  status <- system(pandoc_cmd)
  if (!identical(status, 0L)) {
    stop("pandoc failed while rendering the HTML report.", call. = FALSE)
  }
} else {
  stop("pandoc is required to render the HTML report but was not found on PATH.", call. = FALSE)
}

cat("Wrote report to", md_path, "\n")
cat("Wrote HTML report to", html_path, "\n")
cat("Wrote figures to", figure_dir, "\n")
