#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(ggplot2)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    output_dir = "results/rebuttal/pbmcsca_full_multimethod",
    nfeatures = 2000L,
    dims = 30L,
    seed = 42L,
    pt_size = 0.2
  )
  if (length(args) == 0) {
    return(out)
  }
  for (a in args) {
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    key <- kv[1]
    val <- kv[2]
    if (key == "--output-dir") out$output_dir <- val
    if (key == "--nfeatures") out$nfeatures <- as.integer(val)
    if (key == "--dims") out$dims <- as.integer(val)
    if (key == "--seed") out$seed <- as.integer(val)
    if (key == "--pt-size") out$pt_size <- as.numeric(val)
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

save_plot <- function(p, path, width = 10, height = 8, dpi = 300) {
  ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi)
}

run_raw_umap <- function(obj, nfeatures, dims_use) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = max(dims_use), verbose = FALSE)
  obj <- RunUMAP(
    obj,
    dims = dims_use,
    reduction = "pca",
    reduction.name = "umap_raw",
    reduction.key = "UMAPRAW_",
    verbose = FALSE
  )
  obj
}

build_split_list <- function(obj, nfeatures) {
  obj_list <- SplitObject(obj, split.by = "Method")
  obj_list <- lapply(
    obj_list,
    function(x) {
      DefaultAssay(x) <- "RNA"
      x <- NormalizeData(x, verbose = FALSE)
      x <- FindVariableFeatures(x, nfeatures = nfeatures, verbose = FALSE)
      x
    }
  )
  obj_list
}

run_cca_umap <- function(obj_list, dims_use) {
  anchors <- FindIntegrationAnchors(
    object.list = obj_list,
    dims = dims_use
  )
  int_obj <- IntegrateData(anchorset = anchors, dims = dims_use)
  DefaultAssay(int_obj) <- "integrated"
  int_obj <- ScaleData(int_obj, verbose = FALSE)
  int_obj <- RunPCA(int_obj, npcs = max(dims_use), verbose = FALSE)
  int_obj <- RunUMAP(
    int_obj,
    dims = dims_use,
    reduction = "pca",
    reduction.name = "umap_cca",
    reduction.key = "UMAPCCA_",
    verbose = FALSE
  )
  int_obj
}

run_rpca_umap <- function(obj_list, dims_use, nfeatures) {
  features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures)
  obj_list <- lapply(
    obj_list,
    function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, npcs = max(dims_use), verbose = FALSE)
      x
    }
  )
  anchors <- FindIntegrationAnchors(
    object.list = obj_list,
    anchor.features = features,
    reduction = "rpca",
    dims = dims_use
  )
  int_obj <- IntegrateData(anchorset = anchors, dims = dims_use)
  DefaultAssay(int_obj) <- "integrated"
  int_obj <- ScaleData(int_obj, verbose = FALSE)
  int_obj <- RunPCA(int_obj, npcs = max(dims_use), verbose = FALSE)
  int_obj <- RunUMAP(
    int_obj,
    dims = dims_use,
    reduction = "pca",
    reduction.name = "umap_rpca",
    reduction.key = "UMAPRPCA_",
    verbose = FALSE
  )
  int_obj
}

run_harmony_umap <- function(raw_obj, dims_use) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    return(NULL)
  }
  obj <- raw_obj
  obj <- harmony::RunHarmony(
    object = obj,
    group.by.vars = "Method",
    reduction = "pca",
    dims.use = dims_use,
    verbose = FALSE
  )
  obj <- RunUMAP(
    obj,
    reduction = "harmony",
    dims = dims_use,
    reduction.name = "umap_harmony",
    reduction.key = "UMAPHARM_",
    verbose = FALSE
  )
  obj
}

main <- function() {
  args <- parse_args()
  set.seed(args$seed)
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }
  dims_use <- seq_len(max(2L, args$dims))

  message("[1/6] Loading full pbmcsca...")
  pb <- check_and_load_pbmcsca()
  message(sprintf("  - Cells: %d", ncol(pb)))
  message(sprintf("  - Genes: %d", nrow(pb)))

  batch_tab <- as.data.frame(table(pb$Method), stringsAsFactors = FALSE)
  colnames(batch_tab) <- c("batch_definition_method", "n_cells")
  write.csv(
    batch_tab,
    file.path(args$output_dir, "pbmcsca_batch_definition_counts.csv"),
    row.names = FALSE
  )

  message("[2/6] Uncorrected (raw) UMAP...")
  pb_raw <- run_raw_umap(pb, nfeatures = args$nfeatures, dims_use = dims_use)

  p_raw_batch <- DimPlot(
    pb_raw,
    reduction = "umap_raw",
    group.by = "Method",
    pt.size = args$pt_size
  ) + ggtitle("Original Uncorrected (group.by = Method)")
  p_raw_cell <- DimPlot(
    pb_raw,
    reduction = "umap_raw",
    group.by = "CellType",
    pt.size = args$pt_size,
    label = TRUE,
    repel = TRUE
  ) + ggtitle("Original Uncorrected (group.by = CellType)")

  save_plot(p_raw_batch, file.path(args$output_dir, "umap_raw_by_method.png"))
  save_plot(p_raw_cell, file.path(args$output_dir, "umap_raw_by_celltype.png"), width = 11, height = 8)

  message("[3/6] Preparing split object list...")
  pb_list <- build_split_list(pb, nfeatures = args$nfeatures)

  message("[4/6] Seurat CCA integration...")
  pb_cca <- run_cca_umap(pb_list, dims_use = dims_use)
  p_cca_batch <- DimPlot(
    pb_cca,
    reduction = "umap_cca",
    group.by = "Method",
    pt.size = args$pt_size
  ) + ggtitle("Seurat CCA (group.by = Method)")
  save_plot(p_cca_batch, file.path(args$output_dir, "umap_cca_by_method.png"))

  message("[5/6] Seurat RPCA integration...")
  pb_rpca <- run_rpca_umap(pb_list, dims_use = dims_use, nfeatures = args$nfeatures)
  p_rpca_batch <- DimPlot(
    pb_rpca,
    reduction = "umap_rpca",
    group.by = "Method",
    pt.size = args$pt_size
  ) + ggtitle("Seurat RPCA (group.by = Method)")
  save_plot(p_rpca_batch, file.path(args$output_dir, "umap_rpca_by_method.png"))

  message("[6/6] Harmony integration (optional)...")
  pb_harmony <- run_harmony_umap(pb_raw, dims_use = dims_use)
  p_harmony_batch <- NULL
  if (!is.null(pb_harmony)) {
    p_harmony_batch <- DimPlot(
      pb_harmony,
      reduction = "umap_harmony",
      group.by = "Method",
      pt.size = args$pt_size
    ) + ggtitle("Harmony (group.by = Method)")
    save_plot(p_harmony_batch, file.path(args$output_dir, "umap_harmony_by_method.png"))
  } else {
    message("  - Package 'harmony' not installed; skip Harmony panel.")
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    panel_list <- list(p_raw_batch, p_cca_batch, p_rpca_batch)
    if (!is.null(p_harmony_batch)) {
      panel_list <- c(panel_list, list(p_harmony_batch))
    }
    panel <- patchwork::wrap_plots(panel_list, ncol = 2)
    save_plot(panel, file.path(args$output_dir, "umap_multimethod_panel_by_method.png"), width = 16, height = 12)
  } else {
    message("  - package 'patchwork' not installed; skip combined panel.")
  }

  message("Done.")
  message(sprintf("Output dir: %s", args$output_dir))
}

main()
