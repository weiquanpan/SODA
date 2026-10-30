suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(scDesign3)
  library(SingleCellExperiment)
  library(Matrix)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    output = "data/simulated_scdesign3_benchmark.h5ad",
    seed = 123L,
    nfeatures = 1000L,
    cells_per_group = 300L,
    min_cells_per_group = 120L,
    n_cores = 4L,
    mu_formula = "CellType + Method",
    sigma_formula = "1",
    copula = "gaussian",
    other_covariates = "Method",
    condition_label = "simulated"
  )

  if (length(args) == 0) {
    return(out)
  }

  for (a in args) {
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) {
      next
    }
    key <- kv[1]
    val <- kv[2]
    if (key == "--output") out$output <- val
    if (key == "--seed") out$seed <- as.integer(val)
    if (key == "--nfeatures") out$nfeatures <- as.integer(val)
    if (key == "--cells-per-group") out$cells_per_group <- as.integer(val)
    if (key == "--min-cells-per-group") out$min_cells_per_group <- as.integer(val)
    if (key == "--n-cores") out$n_cores <- as.integer(val)
    if (key == "--mu-formula") out$mu_formula <- val
    if (key == "--sigma-formula") out$sigma_formula <- val
    if (key == "--copula") out$copula <- val
    if (key == "--other-covariates") out$other_covariates <- val
    if (key == "--condition-label") out$condition_label <- val
  }
  out
}

check_and_load_pbmcsca <- function() {
  installed <- tryCatch(
    SeuratData::InstalledData(),
    error = function(e) data.frame()
  )

  if (!("Dataset" %in% colnames(installed)) || !("pbmcsca" %in% installed$Dataset)) {
    SeuratData::InstallData("pbmcsca")
  }
  data("pbmcsca")
  UpdateSeuratObject(pbmcsca)
}

downsample_balanced <- function(obj, methods, min_cells_per_group, cells_per_group, seed) {
  tab <- table(obj$CellType, obj$Method)
  fallback_candidates <- unique(c(
    as.integer(min_cells_per_group),
    120L, 100L, 80L, 60L, 40L, 20L
  ))
  fallback_candidates <- sort(fallback_candidates, decreasing = TRUE)

  used_threshold <- NA_integer_
  keep_ct <- character(0)
  for (thr in fallback_candidates) {
    ct_try <- rownames(tab)[
      rowSums(tab[, methods, drop = FALSE] >= thr) == length(methods)
    ]
    if (length(ct_try) >= 2) {
      used_threshold <- thr
      keep_ct <- ct_try
      break
    }
  }
  if (length(keep_ct) < 2) {
    stop(
      "Too few shared cell types even after fallback thresholds. ",
      "Try changing method pair or lowering --cells-per-group."
    )
  }
  message(sprintf("  - Shared cell-type threshold used: %d", used_threshold))
  message(sprintf("  - Shared cell types kept: %d", length(keep_ct)))

  obj2 <- subset(obj, Method %in% methods & CellType %in% keep_ct)
  set.seed(seed)
  keep_cells <- character(0)

  for (ct in keep_ct) {
    for (m in methods) {
      idx <- which(obj2$CellType == ct & obj2$Method == m)
      cells <- colnames(obj2)[idx]
      take_n <- min(length(cells), cells_per_group)
      if (take_n < 20) {
        next
      }
      keep_cells <- c(keep_cells, sample(cells, size = take_n, replace = FALSE))
    }
  }
  subset(obj2, cells = unique(keep_cells))
}

main <- function() {
  args <- parse_args()
  set.seed(args$seed)
  py_bin <- Sys.which("python")
  if (!nzchar(py_bin)) {
    stop("Python executable not found in PATH. Please ensure `python` is available.")
  }

  message("[1/5] Loading pbmcsca...")
  pb <- check_and_load_pbmcsca()

  methods_use <- c("10x Chromium (v2)", "Smart-seq2")
  message("[2/5] Building balanced training subset...")
  pb_sub <- downsample_balanced(
    obj = pb,
    methods = methods_use,
    min_cells_per_group = args$min_cells_per_group,
    cells_per_group = args$cells_per_group,
    seed = args$seed
  )

  message("[3/5] Selecting HVGs and fitting scDesign3...")
  pb_sub <- FindVariableFeatures(pb_sub, selection.method = "vst", nfeatures = args$nfeatures)
  hvgs <- VariableFeatures(pb_sub)
  pb_sub <- pb_sub[hvgs, ]
  sce <- as.SingleCellExperiment(pb_sub)
  sce$CellType <- as.factor(sce$CellType)
  sce$Method <- as.factor(sce$Method)

  other_cov <- args$other_covariates
  if (is.na(other_cov) || other_cov == "" || toupper(other_cov) == "NONE") {
    other_cov <- NULL
  }

  simu <- scdesign3(
    sce = sce,
    assay_use = "counts",
    celltype = "CellType",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = other_cov,
    mu_formula = args$mu_formula,
    sigma_formula = args$sigma_formula,
    corr_formula = "1",
    family_use = "nb",
    n_cores = args$n_cores,
    copula = args$copula,
    DT = TRUE,
    pseudo_obs = FALSE,
    return_model = FALSE
  )

  message("[4/5] Preparing AnnData-compatible matrices/metadata...")
  sim_counts <- as.matrix(simu$new_count)
  sim_meta <- as.data.frame(simu$new_covariate, stringsAsFactors = FALSE)
  rownames(sim_counts) <- rownames(sce)
  if (!is.null(rownames(sim_meta)) && nrow(sim_meta) == ncol(sim_counts)) {
    colnames(sim_counts) <- rownames(sim_meta)
  }

  if (!("Method" %in% colnames(sim_meta))) {
    method_from_sce <- NULL
    if (!is.null(rownames(sim_meta))) {
      idx <- match(rownames(sim_meta), colnames(sce))
      method_from_sce <- as.character(sce$Method[idx])
    }
    if (is.null(method_from_sce) || all(is.na(method_from_sce))) {
      method_from_sce <- rep("batch0", nrow(sim_meta))
    }
    sim_meta$Method <- method_from_sce
  }

  method_levels <- unique(sim_meta$Method)
  sim_meta$batch <- as.integer(factor(sim_meta$Method, levels = method_levels)) - 1L
  sim_meta$cell_type <- as.character(sim_meta$CellType)
  sim_meta$condition <- args$condition_label
  sim_meta$leiden <- sim_meta$cell_type

  counts_sp <- Matrix(sim_counts, sparse = TRUE)
  lib <- Matrix::colSums(counts_sp)
  norm <- t(t(counts_sp) / pmax(lib, 1)) * 1e4
  lognorm <- log1p(norm)
  gene_names <- rownames(sim_counts)
  cell_names <- colnames(sim_counts)
  if (is.null(cell_names) || length(cell_names) != ncol(sim_counts)) {
    cell_names <- paste0("cell_", seq_len(ncol(sim_counts)))
    colnames(sim_counts) <- cell_names
    colnames(counts_sp) <- cell_names
    colnames(lognorm) <- cell_names
  }

  # AnnData requires cells x genes; convert dgCMatrix to scipy.csr_matrix.
  counts_cxg <- as(t(counts_sp), "dgCMatrix")
  lognorm_cxg <- as(t(lognorm), "dgCMatrix")
  obs_df <- sim_meta
  rownames(obs_df) <- cell_names
  var_df <- data.frame(
    highly_variable = rep(TRUE, length(gene_names)),
    feature_name = gene_names,
    row.names = gene_names,
    stringsAsFactors = FALSE
  )

  out_dir <- dirname(args$output)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  message("[5/5] Writing h5ad...")
  tmp_dir <- tempfile("scdesign3_export_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  counts_mtx <- file.path(tmp_dir, "counts_cxg.mtx")
  obs_csv <- file.path(tmp_dir, "obs.csv")
  var_csv <- file.path(tmp_dir, "var.csv")
  py_script <- file.path(tmp_dir, "write_h5ad.py")

  Matrix::writeMM(counts_cxg, counts_mtx)
  write.csv(obs_df, obs_csv, row.names = TRUE)
  write.csv(var_df, var_csv, row.names = TRUE)

  py_code <- c(
    "import argparse",
    "import pandas as pd",
    "import scanpy as sc",
    "from scipy import sparse, io",
    "",
    "def main():",
    "    p = argparse.ArgumentParser()",
    "    p.add_argument('--counts_mtx', required=True)",
    "    p.add_argument('--obs_csv', required=True)",
    "    p.add_argument('--var_csv', required=True)",
    "    p.add_argument('--output', required=True)",
    "    args = p.parse_args()",
    "",
    "    counts = io.mmread(args.counts_mtx).tocsr().astype('float64')",
    "    obs = pd.read_csv(args.obs_csv, index_col=0)",
    "    var = pd.read_csv(args.var_csv, index_col=0)",
    "",
    "    adata = sc.AnnData(X=counts.copy(), obs=obs, var=var)",
    "    adata.layers['counts'] = counts.copy()",
    "    sc.pp.normalize_total(adata, target_sum=1e4)",
    "    sc.pp.log1p(adata)",
    "",
    "    top_hvg = min(2000, adata.n_vars)",
    "    sc.pp.highly_variable_genes(",
    "        adata,",
    "        n_top_genes=top_hvg,",
    "        flavor='cell_ranger',",
    "        layer='counts',",
    "    )",
    "    if int(adata.var['highly_variable'].sum()) == 0:",
    "        adata.var['highly_variable'] = True",
    "",
    "    sc.pp.scale(adata, max_value=10.0)",
    "    sc.tl.pca(adata, n_comps=min(50, adata.n_vars), random_state=42)",
    "    sc.pp.neighbors(adata, n_pcs=min(30, adata.obsm['X_pca'].shape[1]), random_state=42)",
    "    sc.tl.umap(adata, random_state=42)",
    "    sc.tl.leiden(adata, resolution=0.5, random_state=42)",
    "    adata.write_h5ad(args.output)",
    "",
    "if __name__ == '__main__':",
    "    main()"
  )
  writeLines(py_code, con = py_script, useBytes = TRUE)

  py_args <- c(
    py_script,
    "--counts_mtx", counts_mtx,
    "--obs_csv", obs_csv,
    "--var_csv", var_csv,
    "--output", args$output
  )
  status <- system2(py_bin, args = py_args, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Failed to convert matrices to h5ad via Python helper script.")
  }

  message("Done.")
  message(sprintf("  output: %s", args$output))
  message(sprintf("  shape: %d cells x %d genes", nrow(obs_df), nrow(var_df)))
  message(sprintf("  mu_formula: %s", args$mu_formula))
  message(sprintf("  copula: %s", args$copula))
  message(sprintf("  other_covariates: %s", ifelse(is.null(other_cov), "NULL", args$other_covariates)))
  message(sprintf("  batches: %s", paste(sort(unique(sim_meta$batch)), collapse = ", ")))
  message(sprintf("  cell types: %d", length(unique(sim_meta$cell_type))))
}

main()