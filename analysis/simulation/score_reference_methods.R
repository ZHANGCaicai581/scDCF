#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(AUCell)
  library(VAM)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript score_reference_methods.R SCENARIO_DIR")
scenario_dir <- normalizePath(args[[1]], mustWork = TRUE)
data_dir <- file.path(scenario_dir, "data")
out_dir <- file.path(scenario_dir, "reference_scores")
dir.create(out_dir, showWarnings = FALSE)
meta <- read.csv(file.path(data_dir, "metadata.csv"))
genes <- read.csv(file.path(data_dir, "genes.csv"))$Gene_name
magma <- read.csv(file.path(data_dir, "magma_genes_with_z.csv"))
counts <- as(readMM(file.path(data_dir, "counts.mtx")), "dgCMatrix")
stopifnot(identical(dim(counts), c(length(genes), nrow(meta))))
dimnames(counts) <- list(genes, meta$Cell_ID)
trait <- intersect(magma$Gene_name, genes)
stopifnot(length(trait) == nrow(magma))
obj <- CreateSeuratObject(counts = counts)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- AddModuleScore(obj, features = list(trait), name = "Trait", nbin = 24, ctrl = 100, seed = 1)
save_score <- function(method, score) {
  stopifnot(length(score) == nrow(meta), all(is.finite(score)))
  write.csv(data.frame(cell_id = meta$Cell_ID, score = as.numeric(score)),
            file.path(out_dir, paste0(method, ".csv")), row.names = FALSE)
}
save_score("Seurat", obj$Trait1)
message("Seurat completed")

# AUCell 1.20 passes useNames=NA to newer matrixStats; names do not affect ranks.
trace("colRanks", where = asNamespace("matrixStats"), print = FALSE,
      tracer = quote(if (length(useNames) == 1L && is.na(useNames)) useNames <- TRUE))
set.seed(1)
rankings <- AUCell_buildRankings(as.matrix(counts), plotStats = FALSE, verbose = FALSE,
                                splitByBlocks = FALSE)
auc <- AUCell_calcAUC(list(trait = trait), rankings, aucMaxRank = ceiling(0.05 * nrow(counts)),
                      nCores = 1, normAUC = TRUE, verbose = FALSE)
save_score("AUCell", getAUC(auc)["trait", meta$Cell_ID])
untrace("colRanks", where = asNamespace("matrixStats"))
message("Reference AUCell completed")

# Retain the earlier full-rank statistic under its own, explicit name.
legacy <- vapply(seq_len(ncol(counts)), function(i) {
  ranks <- rank(-as.numeric(counts[, i]), ties.method = "average")
  mean((nrow(counts) - ranks[match(trait, genes)] + 1) / nrow(counts))
}, numeric(1))
save_score("LegacyFullRank", legacy)
norm <- GetAssayData(obj, assay = "RNA", slot = "data")
vam <- VAM::vam(gene.expr = t(as.matrix(norm[trait, , drop = FALSE])), center = TRUE, gamma = FALSE)
save_score("VAM", vam$cdf.value)
save_score("VAM_distance", vam$distance.sq)
message("VAM completed")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
write_json(list(complete = TRUE, trait_genes = length(trait), seurat_seed = 1,
                aucell_seed = 1, aucell_max_rank = ceiling(0.05 * nrow(counts)),
                aucell_implementation = "AUCell package; normalized recovery AUC in top 5% of genes",
                legacy_implementation = "Earlier full-rank mean statistic; not reference AUCell",
                vam = "VAM::vam(center=TRUE, gamma=FALSE), cdf.value",
                matrixStats_compatibility = "Process-local colRanks useNames NA -> TRUE"),
           file.path(out_dir, "complete.json"), auto_unbox = TRUE, pretty = TRUE)
