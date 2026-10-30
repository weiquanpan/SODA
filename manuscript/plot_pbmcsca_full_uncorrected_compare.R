#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(ggplot2)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    output_dir = "results/rebuttal/pbmcsca_full",
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
    if (length(kv) != 2) {
      next
    }
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

save_dimplot <- function(plot_obj, out_path, width = 10, height = 8, dpi = 300) {
  ggsave(
    filename = out_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
}

main <- function() {
  args <- parse_args()
  set.seed(args$seed)
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }

  message("[1/4] Loading full pbmcsca...")
  pbmcsca <- check_and_load_pbmcsca()
  message(sprintf("  - Cells: %d", ncol(pbmcsca)))
  message(sprintf("  - Genes: %d", nrow(pbmcsca)))
  message(sprintf("  - Methods: %s", paste(sort(unique(pbmcsca$Method)), collapse = ", ")))

  dims_use <- seq_len(max(2L, args$dims))

  message("[2/4] Building UNCORRECTED (raw) PCA/UMAP...")
  DefaultAssay(pbmcsca) <- "RNA"
  pbmcsca <- NormalizeData(pbmcsca, verbose = FALSE)
  pbmcsca <- FindVariableFeatures(pbmcsca, nfeatures = args$nfeatures, verbose = FALSE)
  pbmcsca <- ScaleData(pbmcsca, verbose = FALSE)
  pbmcsca <- RunPCA(pbmcsca, npcs = max(dims_use), verbose = FALSE)
  pbmcsca <- RunUMAP(
    pbmcsca,
    dims = dims_use,
    reduction = "pca",
    reduction.name = "umap_raw",
    reduction.key = "UMAPRAW_",
    verbose = FALSE
  )

  p_raw_batch <- DimPlot(
    pbmcsca,
    reduction = "umap_raw",
    group.by = "Method",
    pt.size = args$pt_size
  ) + ggtitle("Original Uncorrected (PBMC Cross-platform)")

  p_raw_celltype <- DimPlot(
    pbmcsca,
    reduction = "umap_raw",
    group.by = "CellType",
    pt.size = args$pt_size,
    label = TRUE,
    repel = TRUE
  ) + ggtitle("Original Uncorrected (colored by CellType)")

  save_dimplot(
    p_raw_batch,
    file.path(args$output_dir, "pbmcsca_full_umap_raw_by_method.png"),
    width = 10,
    height = 8
  )
  save_dimplot(
    p_raw_celltype,
    file.path(args$output_dir, "pbmcsca_full_umap_raw_by_celltype.png"),
    width = 11,
    height = 8
  )

  message("[3/4] Building Seurat CCA integrated PCA/UMAP for comparison...")
  pb_list <- SplitObject(pbmcsca, split.by = "Method")
  pb_list <- lapply(
    pb_list,
    function(x) {
      x <- NormalizeData(x, verbose = FALSE)
      x <- FindVariableFeatures(x, nfeatures = args$nfeatures, verbose = FALSE)
      x
    }
  )

  anchors <- FindIntegrationAnchors(
    object.list = pb_list,
    dims = dims_use
  )
  pb_int <- IntegrateData(anchorset = anchors, dims = dims_use)
  DefaultAssay(pb_int) <- "integrated"
  pb_int <- ScaleData(pb_int, verbose = FALSE)
  pb_int <- RunPCA(pb_int, npcs = max(dims_use), verbose = FALSE)
  pb_int <- RunUMAP(
    pb_int,
    dims = dims_use,
    reduction = "pca",
    reduction.name = "umap_cca",
    reduction.key = "UMAPCCA_",
    verbose = FALSE
  )

  p_cca_batch <- DimPlot(
    pb_int,
    reduction = "umap_cca",
    group.by = "Method",
    pt.size = args$pt_size
  ) + ggtitle("Seurat CCA Integrated")

  p_cca_celltype <- DimPlot(
    pb_int,
    reduction = "umap_cca",
    group.by = "CellType",
    pt.size = args$pt_size,
    label = TRUE,
    repel = TRUE
  ) + ggtitle("Seurat CCA Integrated (colored by CellType)")

  save_dimplot(
    p_cca_batch,
    file.path(args$output_dir, "pbmcsca_full_umap_cca_by_method.png"),
    width = 10,
    height = 8
  )
  save_dimplot(
    p_cca_celltype,
    file.path(args$output_dir, "pbmcsca_full_umap_cca_by_celltype.png"),
    width = 11,
    height = 8
  )

  batch_tab <- as.data.frame(table(pbmcsca$Method), stringsAsFactors = FALSE)
  colnames(batch_tab) <- c("method", "n_cells")
  write.csv(
    batch_tab,
    file.path(args$output_dir, "pbmcsca_full_method_counts.csv"),
    row.names = FALSE
  )

  message("[4/4] Trying side-by-side panel export...")
  if (requireNamespace("patchwork", quietly = TRUE)) {
    panel <- p_raw_batch + p_cca_batch + patchwork::plot_layout(ncol = 2)
    save_dimplot(
      panel,
      file.path(args$output_dir, "pbmcsca_full_raw_vs_cca_method_panel.png"),
      width = 16,
      height = 7
    )
  } else {
    message("  - patchwork not installed; skipped panel plot.")
  }

  message("Done.")
  message(sprintf("Output dir: %s", args$output_dir))
}

main()
