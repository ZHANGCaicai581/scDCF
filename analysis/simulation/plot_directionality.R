#!/usr/bin/env Rscript
options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript plot_directionality.R OUTPUT_ROOT [SEED=46]\n")
  quit(status = if (length(args)) 0L else 1L)
}
root <- normalizePath(args[[1]], mustWork = TRUE)
seed <- if (length(args) > 1L) as.integer(args[[2]]) else 46L
out <- file.path(root, "figures")
dir.create(out, showWarnings = FALSE)
status <- read_json(file.path(root, "run_status.json"), simplifyVector = TRUE)
stopifnot(isTRUE(status$complete))
fractions <- c(40L, 60L)
methods <- c("scDRS", "Seurat")
source_rows <- list()

for (fraction in fractions) {
  scenario <- file.path(root,
    sprintf("T2_mixed%ddown_clean_dense_geneDelta_seed%d", fraction, seed))
  meta <- read.csv(file.path(scenario, "data/metadata.csv"), stringsAsFactors = FALSE)
  stopifnot(nrow(meta) == 2400L, !anyDuplicated(meta$Cell_ID))
  for (method in methods) {
    if (method == "scDRS") {
      score <- read.delim(gzfile(file.path(scenario, "scdrs_out/trait1.score.gz")),
                          check.names = FALSE)
      scored <- data.frame(cell_id = score[[1]], score = score$norm_score)
    } else {
      scored <- read.csv(file.path(scenario, "reference_scores/Seurat.csv"))
    }
    stopifnot(nrow(scored) == 2400L, !anyDuplicated(scored$cell_id),
              setequal(scored$cell_id, meta$Cell_ID), all(is.finite(scored$score)))
    aligned <- scored$score[match(meta$Cell_ID, scored$cell_id)]
    keep <- meta$Cell_Type == "CellTypeC"
    part <- data.frame(cell_id = meta$Cell_ID[keep],
      condition = meta$Condition[keep], method = method,
      down_percent = fraction, score = aligned[keep],
      status = ifelse(meta$TruePositive[keep] == "TP", "Perturbed", "Unperturbed"))
    stopifnot(nrow(part) == 800L, sum(part$status == "Perturbed") == 40L,
              sum(part$status == "Unperturbed") == 760L)
    source_rows[[length(source_rows) + 1L]] <- part
  }
}
source_data <- do.call(rbind, source_rows)
write.csv(source_data, file.path(out, "score_distribution_source_data.csv"), row.names = FALSE)

# The fixed display grid aggregates no cells out of the empirical distributions.
cdf_rows <- list()
limits <- list()
for (method in methods) {
  selected <- source_data[source_data$method == method, ]
  ticks <- pretty(range(selected$score), n = 4)
  limits[[method]] <- list(range = range(ticks), ticks = ticks)
  grid <- seq(min(ticks), max(ticks), length.out = 201L)
  for (fraction in fractions) {
    for (group in c("Unperturbed", "Perturbed")) {
      s <- selected$score[selected$down_percent == fraction & selected$status == group]
      cdf_rows[[length(cdf_rows) + 1L]] <- data.frame(method = method,
        down_percent = fraction, status = group, score = grid,
        cumulative_fraction = ecdf(s)(grid), n_cells = length(s))
    }
  }
}
cdf <- do.call(rbind, cdf_rows)
cdf$status <- factor(cdf$status, levels = c("Unperturbed", "Perturbed"))
cdf$setting <- factor(paste0(cdf$down_percent, "%"), levels = c("40%", "60%"))
write.csv(cdf, file.path(out, "ecdf_display_source_data.csv"), row.names = FALSE)

score_colours <- c(Unperturbed = "#397E98", Perturbed = "#C57842")
direction_colours <- c(Down = "#397E98", Unchanged = "#D8DDE0", Up = "#C57842")
theme_set(theme_classic(base_size = 8, base_family = "Helvetica") + theme(
  axis.line = element_line(linewidth = 0.35),
  axis.ticks = element_line(linewidth = 0.3),
  axis.text = element_text(size = 7, colour = "#272727"),
  legend.text = element_text(size = 7), legend.title = element_blank(),
  strip.text = element_text(size = 8, face = "bold"),
  plot.title = element_text(size = 8, face = "bold"),
  plot.tag = element_text(size = 9, face = "bold"),
  plot.margin = margin(5, 5, 5, 5),
  legend.position = "bottom"))

plots <- lapply(methods, function(method) {
  ggplot(cdf[cdf$method == method, ], aes(score, cumulative_fraction,
      colour = status, linetype = setting,
      group = interaction(status, setting))) +
    geom_line(linewidth = 0.65) +
    scale_colour_manual(values = score_colours) +
    scale_linetype_manual(values = c("40%" = "solid", "60%" = "dashed")) +
    scale_x_continuous(limits = limits[[method]]$range,
                       breaks = limits[[method]]$ticks, expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0, 0.01))) +
    labs(x = "Native score", y = "Cumulative fraction", title = method)
})
fig1 <- wrap_plots(plots, nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")

totals <- read.csv(file.path(root, "observed_direction_totals.csv"))
stopifnot(identical(as.integer(totals$down_percent), fractions),
          all(totals$selected_genes == 14000L))
direction_rows <- list()
for (i in seq_len(nrow(totals))) {
  r <- totals[i, ]
  for (stage in c("Intended log-space effects", "Observed count changes")) {
    counts <- if (stage == "Intended log-space effects") {
      c(Down = r$intended_down, Unchanged = 0, Up = r$intended_up)
    } else c(Down = r$observed_down, Unchanged = r$observed_unchanged, Up = r$observed_up)
    stopifnot(sum(counts) == 14000L)
    direction_rows[[length(direction_rows) + 1L]] <- data.frame(
      stage = stage, down_percent = r$down_percent, direction = names(counts),
      count = unname(counts), fraction = unname(counts) / 14000)
  }
}
directions <- do.call(rbind, direction_rows)
directions$direction <- factor(directions$direction, levels = c("Down", "Unchanged", "Up"))
directions$setting <- factor(paste0(directions$down_percent, "% down"),
                             levels = c("60% down", "40% down"))
directions$label <- ifelse(directions$count > 0,
  sprintf("%s\n%.1f%%", format(directions$count, big.mark = ",", trim = TRUE),
          100 * directions$fraction), "")
write.csv(directions, file.path(out, "direction_counts_source_data.csv"), row.names = FALSE)
bar_plots <- lapply(unique(directions$stage), function(stage) {
  d <- directions[directions$stage == stage, ]
  ggplot(d, aes(fraction, setting, fill = direction)) +
    geom_col(width = 0.55, position = position_stack(reverse = TRUE)) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5, reverse = TRUE),
              size = 2.5, lineheight = 0.95, colour = "#1E2930") +
    scale_fill_manual(values = direction_colours, drop = FALSE) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "50", "100"), expand = expansion(mult = 0)) +
    labs(x = "Selected gene-cell effects (%)", y = NULL, title = stage) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
})
fig2 <- wrap_plots(bar_plots, nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")

# One vector export and one preview; no duplicate TikZ/TIFF implementations.
save_figure <- function(plot, name, height_mm) {
  grDevices::pdf(file.path(out, paste0(name, ".pdf")), width = 183 / 25.4,
                height = height_mm / 25.4, family = "Helvetica", useDingbats = FALSE)
  print(plot)
  dev.off()
  ggsave(file.path(out, paste0(name, ".png")), plot, width = 183,
         height = height_mm, units = "mm", dpi = 300)
}
save_figure(fig1, "matched_score_distributions", 84)
save_figure(fig2, "intended_observed_perturbations", 84)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("Figures and their source data written to ", out)
