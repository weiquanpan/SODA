#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(Matrix)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    output = "data/pbmcsca_real_10xv2_smartseq2.h5ad",
    nfeatures = 2000L,
    npcs = 30L,
    seed = 42L
  )
  if (length(args) == 0) {
    return(out)
  }
  for (a in args) {
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    key <- kv[1]
    val <- kv[2]
    if (key == "--output") out$output <- val
    if (key == "--nfeatures") out$nfeatures <- as.integer(val)
    if (key == "--npcs") out$npcs <- as.integer(val)
    if (key == "--seed") out$seed <- as.integer(val)
  }
  out
}

load_pbmcsca <- function() {
  installed <- tryCatch(SeuratData::InstalledData(), error = function(e) data.frame())
  if (!("Dataset" %in% colnames(installed)) || !("pbmcsca" %in% installed$Dataset)) {
    SeuratData::InstallData("pbmcsca")
  }
  data("pbmcsca")
  UpdateSeuratObject(pbmcsca)
}

main <- function() {
  args <- parse_args()
  set.seed(args$seed)

  message("[1/5] Loading pbmcsca and subsetting two platforms...")
  pb <- load_pbmcsca()
  pb <- subset(pb, Method %in% c("10x Chromium (v2)", "Smart-seq2"))
  DefaultAssay(pb) <- "RNA"

  message("[2/5] Running Seurat preprocessing for raw baseline...")
  pb <- NormalizeData(pb, verbose = FALSE)
  pb <- FindVariableFeatures(pb, nfeatures = args$nfeatures, verbose = FALSE)
  pb <- ScaleData(pb, verbose = FALSE)
  pb <- RunPCA(pb, npcs = args$npcs, verbose = FALSE)
  pb <- FindNeighbors(pb, dims = seq_len(args$npcs), verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)
  pb <- RunUMAP(pb, dims = seq_len(args$npcs), reduction = "pca", verbose = FALSE)

  message("[3/5] Preparing matrices and metadata...")
  hvg <- VariableFeatures(pb)
  if (length(hvg) < 500) {
    stop("Too few HVGs detected; cannot build robust benchmark dataset.")
  }
  pb_hvg <- pb[hvg, ]

  counts <- GetAssayData(pb_hvg, assay = "RNA", layer = "counts")
  lognorm <- GetAssayData(pb_hvg, assay = "RNA", layer = "data")
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  if (!inherits(lognorm, "dgCMatrix")) lognorm <- as(lognorm, "dgCMatrix")

  counts_cxg <- as(t(counts), "dgCMatrix")
  x_cxg <- as(t(lognorm), "dgCMatrix")
  cell_names <- colnames(pb_hvg)
  gene_names <- rownames(pb_hvg)

  obs <- pb@meta.data
  obs$batch <- as.character(obs$Method)
  obs$cell_type <- as.character(obs$CellType)
  obs$leiden <- as.character(obs$seurat_clusters)
  obs <- obs[cell_names, , drop = FALSE]
  rownames(obs) <- cell_names

  hvg_set <- VariableFeatures(pb_hvg)
  var <- data.frame(
    highly_variable = gene_names %in% hvg_set,
    feature_name = gene_names,
    row.names = gene_names,
    stringsAsFactors = FALSE
  )

  pca <- Embeddings(pb, reduction = "pca")
  pca <- pca[cell_names, seq_len(min(args$npcs, ncol(pca))), drop = FALSE]
  umap <- Embeddings(pb, reduction = "umap")
  umap <- umap[cell_names, 1:2, drop = FALSE]

  out_dir <- dirname(args$output)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  message("[4/5] Exporting temporary files and writing h5ad...")
  tmp_dir <- tempfile("pbmcsca_real_export_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  counts_mtx <- file.path(tmp_dir, "counts_cxg.mtx")
  x_mtx <- file.path(tmp_dir, "x_cxg.mtx")
  obs_csv <- file.path(tmp_dir, "obs.csv")
  var_csv <- file.path(tmp_dir, "var.csv")
  pca_csv <- file.path(tmp_dir, "pca.csv")
  umap_csv <- file.path(tmp_dir, "umap.csv")
  py_script <- file.path(tmp_dir, "write_h5ad.py")

  Matrix::writeMM(counts_cxg, counts_mtx)
  Matrix::writeMM(x_cxg, x_mtx)
  write.csv(obs, obs_csv, row.names = TRUE)
  write.csv(var, var_csv, row.names = TRUE)
  write.csv(pca, pca_csv, row.names = TRUE)
  write.csv(umap, umap_csv, row.names = TRUE)

  py_code <- c(
    "import argparse",
    "import pandas as pd",
    "import scanpy as sc",
    "from scipy import io",
    "",
    "def main():",
    "    p = argparse.ArgumentParser()",
    "    p.add_argument('--counts_mtx', required=True)",
    "    p.add_argument('--x_mtx', required=True)",
    "    p.add_argument('--obs_csv', required=True)",
    "    p.add_argument('--var_csv', required=True)",
    "    p.add_argument('--pca_csv', required=True)",
    "    p.add_argument('--umap_csv', required=True)",
    "    p.add_argument('--output', required=True)",
    "    args = p.parse_args()",
    "",
    "    counts = io.mmread(args.counts_mtx).tocsr().astype('float32')",
    "    x = io.mmread(args.x_mtx).tocsr().astype('float32')",
    "    obs = pd.read_csv(args.obs_csv, index_col=0)",
    "    var = pd.read_csv(args.var_csv, index_col=0)",
    "    pca = pd.read_csv(args.pca_csv, index_col=0)",
    "    umap = pd.read_csv(args.umap_csv, index_col=0)",
    "",
    "    adata = sc.AnnData(X=x, obs=obs, var=var)",
    "    adata.layers['counts'] = counts",
    "    adata.obsm['X_pca'] = pca.loc[adata.obs_names].to_numpy(dtype='float32')",
    "    adata.obsm['X_umap'] = umap.loc[adata.obs_names].to_numpy(dtype='float32')",
    "    adata.write_h5ad(args.output)",
    "",
    "if __name__ == '__main__':",
    "    main()"
  )
  writeLines(py_code, con = py_script, useBytes = TRUE)

  py_bin <- Sys.which("python")
  if (!nzchar(py_bin)) {
    stop("Python executable not found in PATH.")
  }
  py_args <- c(
    py_script,
    "--counts_mtx", counts_mtx,
    "--x_mtx", x_mtx,
    "--obs_csv", obs_csv,
    "--var_csv", var_csv,
    "--pca_csv", pca_csv,
    "--umap_csv", umap_csv,
    "--output", args$output
  )
  status <- system2(py_bin, args = py_args, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Failed to export h5ad.")
  }

  message("[5/5] Done.")
  message(sprintf("Output: %s", args$output))
  message(sprintf("Cells: %d, Genes(HVG): %d", ncol(pb_hvg), nrow(pb_hvg)))
  message(sprintf("Methods: %s", paste(sort(unique(pb$Method)), collapse = ", ")))
}

main()
