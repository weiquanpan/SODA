#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(zellkonverter)
})

args <- commandArgs(trailingOnly = TRUE)
output_path <- if (length(args) >= 1) args[[1]] else "outputs/demo_data/pbmc_small_demo.h5ad"

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

data("pbmc_small", package = "SeuratObject")

get_counts <- function(obj) {
  if (exists("LayerData", mode = "function")) {
    tryCatch(
      LayerData(obj, assay = "RNA", layer = "counts"),
      error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
    )
  } else {
    GetAssayData(obj, assay = "RNA", slot = "counts")
  }
}

counts <- get_counts(pbmc_small)
meta <- pbmc_small@meta.data

if ("groups" %in% colnames(meta)) {
  meta$batch <- as.character(meta$groups)
} else {
  meta$batch <- rep(c("batch_a", "batch_b"), length.out = ncol(counts))
}

meta$cell_type <- as.character(Idents(pbmc_small))

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = DataFrame(meta)
)

zellkonverter::writeH5AD(sce, output_path)

cat("Wrote demo dataset to:", output_path, "\n")
cat("Cells:", ncol(counts), "\n")
cat("Genes:", nrow(counts), "\n")
cat("Batch column: batch\n")
cat("Cell-type column: cell_type\n")

