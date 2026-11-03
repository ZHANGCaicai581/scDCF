#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript generate_scenarios.R OUTPUT_ROOT [SEEDS=46] [DOWN_PCTS=10,20] [PAIR_MODE=nested]\n")
  quit(status = if (length(args)) 0L else 1L)
}
suppressPackageStartupMessages({
  library(splatter)
  library(SingleCellExperiment)
  library(Matrix)
  library(jsonlite)
})

output_root <- normalizePath(args[[1]], mustWork = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
seeds <- if (length(args) >= 2L) as.integer(strsplit(args[[2]], ",")[[1]]) else 46L
down_pcts <- if (length(args) >= 3L) as.integer(strsplit(args[[3]], ",")[[1]]) else c(10L, 20L)
pair_mode <- if (length(args) >= 4L) args[[4]] else "nested"
stopifnot(length(seeds) > 0, !anyNA(seeds), !anyDuplicated(seeds))
stopifnot(length(down_pcts) == 2L, !anyNA(down_pcts),
          !anyDuplicated(down_pcts), all(down_pcts > 0 & down_pcts < 100),
          down_pcts[[1]] < down_pcts[[2]], pair_mode %in% c("nested", "sign_reverse"))
if (pair_mode == "sign_reverse") stopifnot(sum(down_pcts) == 100L)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# These helper definitions are unchanged from the author's original R simulator.
# Their source/hash is recorded in provenance.json; no external Downloads file is needed.
mkdirp <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
clip  <- function(x, lo, hi) pmin(pmax(x, lo), hi)

apply_celltype_baseline <- function(counts, genes, cell_type, gene_set,
                                    delta_A = +log(1.25),
                                    delta_B = -log(1.25),
                                    delta_C = 0,
                                    clamp_min = 0, clamp_max = 12) {
  gix <- match(gene_set, genes)
  gix <- gix[!is.na(gix)]
  if (length(gix) == 0) return(counts)

  shift_cells <- function(idx, delta) {
    if (length(idx) == 0 || delta == 0) return()
    logb <- log1p(counts[, idx, drop = FALSE])
    if (delta > 0) {
      logb[gix, ] <- logb[gix, , drop = FALSE] + delta
    } else {
      logb[gix, ] <- logb[gix, , drop = FALSE] - abs(delta)
    }
    logb <- clip(logb, clamp_min, clamp_max)
    counts[, idx] <<- round(expm1(logb))
  }

  shift_cells(which(cell_type == "CellTypeA"), delta_A)
  shift_cells(which(cell_type == "CellTypeB"), delta_B)
  shift_cells(which(cell_type == "CellTypeC"), delta_C)
  counts
}

build_gene_sets <- function(genes, counts, n_causal = 500,
                            magma_size = 800,
                            magma_causal_frac = 1.0,
                            magma_false_frac  = 0.0,
                            seed = 1) {
  set.seed(seed)

  mean_all <- Matrix::rowMeans(counts)
  pool <- names(sort(mean_all, decreasing = TRUE))[800:min(7000, length(mean_all))]
  pool <- pool[!is.na(pool)]
  if (length(pool) < 50) pool <- genes

  causal_genes <- sample(pool, min(n_causal, length(pool)))

  n_magma <- magma_size
  n_from_causal <- min(round(n_magma * magma_causal_frac), length(causal_genes))
  n_false <- max(0, round(n_magma * magma_false_frac))

  magma_from_causal <- if (n_from_causal > 0) sample(causal_genes, n_from_causal) else character(0)
  false_pool <- setdiff(genes, causal_genes)
  magma_false <- if (n_false > 0) sample(false_pool, min(n_false, length(false_pool))) else character(0)

  magma_genes <- unique(c(magma_from_causal, magma_false))
  if (length(magma_genes) > n_magma) magma_genes <- magma_genes[seq_len(n_magma)]

  magma_df <- data.frame(
    Gene_name = magma_genes,
    causal_flag = magma_genes %in% causal_genes,
    stringsAsFactors = FALSE
  )

  # MAGMA weights (positive ZSTAT, monotone decreasing)
  magma_df$ZSTAT <- sort(abs(rnorm(nrow(magma_df), mean = 5, sd = 2)), decreasing = TRUE)

  list(causal_genes = causal_genes, magma_df = magma_df)
}

provenance_dir <- file.path(output_root, "provenance")
dir.create(provenance_dir, showWarnings = FALSE)
snapshot <- file.path(provenance_dir, "generator.R")
if (!file.exists(snapshot)) stopifnot(file.copy(source_path, snapshot))
stopifnot(unname(tools::md5sum(snapshot)) == unname(tools::md5sum(source_path)))
writeLines(capture.output(sessionInfo()), file.path(provenance_dir, "generation_session.txt"))

for (seed in seeds) {
  completed <- vapply(down_pcts, function(pct) {
    file.exists(file.path(output_root,
      sprintf("T2_mixed%ddown_clean_dense_geneDelta_seed%d", pct, seed),
      "data", "generation_complete.json"))
  }, logical(1))
  # Check configuration even when both cached fractions already exist.
  for (pct in down_pcts[completed]) {
    data_dir <- file.path(output_root,
      sprintf("T2_mixed%ddown_clean_dense_geneDelta_seed%d", pct, seed), "data")
    previous <- read_json(file.path(data_dir, "scenario_manifest.json"), simplifyVector = TRUE)
    previous_mode <- if (is.null(previous$pair_mode)) "nested" else previous$pair_mode
    stopifnot(previous_mode == pair_mode, previous$seed == seed,
              identical(as.integer(previous$paired_down_pcts), down_pcts),
              isTRUE(all.equal(previous$tp$mixed_frac_down, pct / 100)))
    marker <- read_json(file.path(data_dir, "generation_complete.json"), simplifyVector = TRUE)
    stopifnot(isTRUE(marker$validated),
              marker$counts_md5 == unname(tools::md5sum(file.path(data_dir, "counts.mtx"))))
  }
  if (all(completed)) {
    message("Both fractions already generated and validated for seed ", seed)
    next
  }
  message("Generating shared baseline for seed ", seed)
  set.seed(seed)
  params <- newSplatParams()
  params <- setParam(params, "nGenes", 12000L)
  params <- setParam(params, "batchCells", c(A = 800L, B = 800L, C = 800L))
  params <- setParam(params, "batch.facLoc", 0.6)
  params <- setParam(params, "batch.facScale", 0.45)
  sce <- splatSimulate(params, verbose = FALSE)
  counts <- as(assay(sce, "counts"), "dgCMatrix")
  genes <- rownames(counts)
  colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
  batch <- as.character(colData(sce)$Batch)
  types <- unname(setNames(c("CellTypeA", "CellTypeB", "CellTypeC"), sort(unique(batch)))[batch])
  condition <- rep("Healthy", ncol(counts))
  for (ct in unique(types)) {
    ix <- which(types == ct)
    condition[sample(ix, length(ix) / 2)] <- "Disease"
  }
  eligible <- which(types == "CellTypeC" & condition == "Disease")
  tp_idx <- sample(eligible, round(length(eligible) * 0.1))
  meta <- data.frame(Cell_ID = colnames(counts), Cell_Type = types,
                     Condition = condition, TruePositive = "Not_TP")
  meta$TruePositive[tp_idx] <- "TP"
  sets <- build_gene_sets(genes, counts, n_causal = 500, magma_size = 800,
                          magma_causal_frac = 1, magma_false_frac = 0,
                          seed = seed + 100L)
  magma <- sets$magma_df
  counts <- apply_celltype_baseline(counts, genes, types, magma$Gene_name,
                                    delta_A = log(1.25), delta_B = -log(1.25), delta_C = 0)
  stopifnot(nrow(meta) == 2400L, length(tp_idx) == 40L, nrow(magma) == 500L)

  # Sample once, then share selected genes, magnitudes, and ordering across fractions.
  set.seed(seed + 7L)
  k <- round(nrow(magma) * 0.7)
  plan <- lapply(tp_idx, function(cc) {
    gene_idx <- match(sample(magma$Gene_name, k), genes)
    data.frame(cell_index = cc, cell_id = colnames(counts)[cc],
               gene_index = gene_idx, gene_name = genes[gene_idx],
               direction_order = sample.int(k),
               delta_magnitude = clip(rnorm(k, log(1.7), 0.25), log(1.1), log(4)))
  })
  shared_plan <- do.call(rbind, plan)
  non_tp <- setdiff(seq_len(ncol(counts)), tp_idx)

  for (down_pct in down_pcts) {
    scenario <- sprintf("T2_mixed%ddown_clean_dense_geneDelta_seed%d", down_pct, seed)
    scenario_dir <- file.path(output_root, scenario)
    data_dir <- file.path(scenario_dir, "data")
    marker <- file.path(data_dir, "generation_complete.json")
    if (file.exists(marker)) {
      previous <- read_json(file.path(data_dir, "scenario_manifest.json"), simplifyVector = TRUE)
      previous_mode <- if (is.null(previous$pair_mode)) "nested" else previous$pair_mode
      stopifnot(previous_mode == pair_mode)
      message("Already generated: ", scenario)
      next
    }
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    message("Applying ", down_pct, "% down-regulation: ", scenario)
    perturbed <- counts
    effects <- shared_plan
    n_down <- as.integer(k * down_pct / 100)
    effects$direction <- ifelse(effects$direction_order <= n_down, "down", "up")
    if (pair_mode == "sign_reverse" && down_pct == down_pcts[[2]]) {
      first_n_down <- as.integer(k * down_pcts[[1]] / 100)
      effects$direction <- ifelse(effects$direction_order <= first_n_down, "up", "down")
    }
    effects$signed_delta <- ifelse(effects$direction == "down", -1, 1) * effects$delta_magnitude
    audits <- vector("list", length(tp_idx))
    for (i in seq_along(tp_idx)) {
      cc <- tp_idx[[i]]
      e <- effects[effects$cell_index == cc, ]
      before <- as.numeric(counts[e$gene_index, cc])
      after <- round(expm1(clip(log1p(before) + e$signed_delta, 0, 12)))
      perturbed[e$gene_index, cc] <- after
      audits[[i]] <- data.frame(
        cell_id = colnames(counts)[cc], intended_up = sum(e$direction == "up"),
        intended_down = sum(e$direction == "down"), selected_genes = nrow(e),
        observed_up = sum(after > before), observed_down = sum(after < before),
        observed_unchanged = sum(after == before))
    }
    audit <- do.call(rbind, audits)
    stopifnot(all(audit$intended_down == n_down), all(audit$intended_up == k - n_down),
              isTRUE(all.equal(perturbed[, non_tp], counts[, non_tp])),
              all(is.finite(perturbed@x)), all(perturbed@x >= 0),
              all(perturbed@x == round(perturbed@x)))
    Matrix::writeMM(perturbed, file.path(data_dir, "counts.mtx"))
    write.csv(data.frame(Gene_name = genes), file.path(data_dir, "genes.csv"), row.names = FALSE)
    write.csv(meta, file.path(data_dir, "metadata.csv"), row.names = FALSE)
    write.csv(magma, file.path(data_dir, "magma_genes_with_z.csv"), row.names = FALSE)
    write.csv(data.frame(Gene_name = sets$causal_genes), file.path(data_dir, "causal_genes.csv"), row.names = FALSE)
    write.csv(meta[tp_idx, "Cell_ID", drop = FALSE], file.path(data_dir, "tp_cells.csv"), row.names = FALSE)
    write.csv(effects, file.path(data_dir, "perturbation_plan.csv"), row.names = FALSE)
    write.csv(audit, file.path(data_dir, "direction_audit.csv"), row.names = FALSE)
    manifest <- list(
      scenario = scenario, seed = seed, paired_baseline_seed = seed,
      pair_mode = pair_mode, paired_down_pcts = down_pcts,
      n_genes = nrow(perturbed), n_cells = ncol(perturbed), n_TP = length(tp_idx),
      n_cells_per_type = list(A = 800, B = 800, C = 800), prop_disease = 0.5,
      prop_TP_in_diseaseC = 0.1,
      splatter = list(batch.facLoc = 0.6, batch.facScale = 0.45, dropout_enable = FALSE),
      magma = list(requested_size = 800, actual_size = nrow(magma), n_causal = 500,
                   magma_causal_frac = 1, magma_false_frac = 0),
      tp = list(tp_mode = "mixed", mixed_frac_down = down_pct / 100,
                mixed_frac_up = 1 - down_pct / 100, tp_gene_frac_per_cell = 0.7,
                intended_down_per_TP = n_down, intended_up_per_TP = k - n_down,
                fc_target = 1.7, delta_mean = log(1.7), gene_specific_delta = TRUE,
                gene_delta_sd = 0.25, delta_lower = log(1.1), delta_upper = log(4),
                sampling = "Fixed per-cell genes, magnitudes, and ordering across fractions",
                exact_sign_reversal = pair_mode == "sign_reverse"),
      confounds = list(add_celltype_baseline = TRUE, delta_A = log(1.25),
                       delta_B = -log(1.25), delta_C = 0, disease_scale_libsize = 1),
      sparsity = list(thinning_keep = 1),
      generator_file = basename(source_path), generator_md5 = unname(tools::md5sum(source_path)))
    write_json(manifest, file.path(data_dir, "scenario_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
    write_json(list(validated = TRUE, counts_md5 = unname(tools::md5sum(file.path(data_dir, "counts.mtx")))),
               marker, auto_unbox = TRUE, pretty = TRUE)
    message("Validated ", n_down, " down and ", k - n_down, " up genes per TP cell")
  }
  rm(sce, counts, perturbed)
  gc()
}
