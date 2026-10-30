# -*- coding: utf-8 -*-
"""SMAI 批次校正方法（通过 rpy2 调用 R SMAI）。"""

import numpy as np
import scanpy as sc
from scipy import sparse

from soda_benchmark.methods.base import BaseBatchMethod, register_method


def _check_rpy2_seurat_smai():
    """检查 rpy2、R Seurat 和 SMAI 是否可用。"""
    try:
        import rpy2.robjects as ro  # noqa: F401
    except ImportError:
        return False, "rpy2 not installed. Install with: pip install rpy2"

    try:
        ro.r('suppressPackageStartupMessages(library(Seurat))')
    except Exception:
        return False, (
            "R package 'Seurat' is not available. "
            "Install in R with: install.packages('Seurat')"
        )

    try:
        ro.r('suppressPackageStartupMessages(library(SMAI))')
    except Exception:
        return False, (
            "R package 'SMAI' is not available. "
            "Install in R with: devtools::install_github('rongstat/SMAI')"
        )

    return True, ""


def _extract_expression_matrix(adata):
    """
    从 AnnData 提取表达矩阵，优先使用原始计数。

    Returns
    -------
    X : np.ndarray, shape (n_cells, n_genes)
    gene_names : list[str]
    """
    X = None
    gene_names = list(adata.var_names)

    if adata.raw is not None:
        try:
            raw_adata = adata.raw.to_adata()
            hvg = [g for g in gene_names if g in raw_adata.var_names]
            if len(hvg) > 0:
                X = raw_adata[:, hvg].X
                gene_names = hvg
        except Exception:
            pass

    if X is None and "counts" in adata.layers:
        X = adata.layers["counts"]

    if X is None:
        X = adata.X

    if sparse.issparse(X):
        X = np.asarray(X.todense(), dtype=np.float64)
    else:
        X = np.asarray(X, dtype=np.float64)

    return X, gene_names


def _r_matrix_to_numpy(ro, r_varname):
    """从 R 全局环境取出矩阵并转为 numpy ndarray。"""
    r_obj = ro.globalenv[r_varname]
    arr = np.array(r_obj, dtype=np.float64)
    if arr.ndim == 2:
        return arr
    dims = ro.r("dim")(r_obj)
    nrow, ncol = int(dims[0]), int(dims[1])
    return arr.reshape((nrow, ncol), order="F")


def _select_common_features(X, gene_names, nfeatures):
    """
    Select a fixed feature set shared across all batches for multi-batch SMAI.
    """
    n_genes = X.shape[1]
    if n_genes <= nfeatures:
        return np.arange(n_genes, dtype=int), list(gene_names)

    # Use variance on log1p(counts) to pick informative genes deterministically.
    X_log = np.log1p(np.asarray(X, dtype=np.float64))
    gene_var = X_log.var(axis=0)
    top_idx = np.argsort(gene_var)[::-1][:nfeatures]
    top_idx = np.sort(top_idx)
    top_genes = [gene_names[i] for i in top_idx]
    return top_idx, top_genes


@register_method("SMAI")
class SMAIMethod(BaseBatchMethod):
    """
    SMAI integration method.

    当前实现按 SMAI 教程流程支持 2 个 batch 的整合。
    """

    def __init__(
        self,
        integrate_key="batch",
        nfeatures=2000,
        t_iter=3,
        knn=30,
        r_max=200,
        prop_align=0.3,
        cutoff=1.001,
        cutoff_t=1.5,
        **kwargs,
    ):
        self.integrate_key = integrate_key
        self.nfeatures = nfeatures
        self.t_iter = t_iter
        self.knn = knn
        self.r_max = r_max
        self.prop_align = prop_align
        self.cutoff = cutoff
        self.cutoff_t = cutoff_t

    @classmethod
    def check_dependencies(cls):
        return _check_rpy2_seurat_smai()

    def _run_smai_pair(self, ro, X_ref, X_tar, gene_names):
        """
        Run SMAI alignment for one (reference, target) batch pair.
        """
        ro.globalenv["x1"] = ro.r.matrix(
            ro.FloatVector(np.ascontiguousarray(X_ref).flatten(order="F")),
            nrow=X_ref.shape[0],
            ncol=X_ref.shape[1],
        )
        ro.globalenv["x2"] = ro.r.matrix(
            ro.FloatVector(np.ascontiguousarray(X_tar).flatten(order="F")),
            nrow=X_tar.shape[0],
            ncol=X_tar.shape[1],
        )
        ro.globalenv["gene_names"] = ro.StrVector(gene_names)
        ro.globalenv["t_iter_val"] = ro.IntVector([int(self.t_iter)])
        ro.globalenv["knn_val"] = ro.IntVector([int(self.knn)])
        ro.globalenv["rmax_val"] = ro.IntVector([int(self.r_max)])
        ro.globalenv["prop_align_val"] = ro.FloatVector([float(self.prop_align)])
        ro.globalenv["cutoff_val"] = ro.FloatVector([float(self.cutoff)])
        ro.globalenv["cutoff_t_val"] = ro.FloatVector([float(self.cutoff_t)])

        ro.r(
            """
            suppressPackageStartupMessages({
              library(Seurat)
              library(SMAI)
            })

            rownames(x1) <- gene_names
            rownames(x2) <- gene_names

            get_norm_data <- function(mat) {
              obj <- CreateSeuratObject(
                counts = mat,
                project = "smai_align",
                min.cells = 0,
                min.features = 0
              )
              obj <- NormalizeData(obj, verbose = FALSE)
              out <- tryCatch(
                as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")),
                error = function(e) {
                  as.matrix(GetAssayData(obj, assay = "RNA", slot = "data"))
                }
              )
              # Robust row alignment: keep target order and fill missing genes with 0.
              out_aligned <- matrix(
                0,
                nrow = length(gene_names),
                ncol = ncol(out),
                dimnames = list(gene_names, colnames(out))
              )
              present <- intersect(gene_names, rownames(out))
              if (length(present) > 0) {
                out_aligned[present, ] <- out[present, , drop = FALSE]
              }
              out <- out_aligned
              return(out)
            }

            data_x1_smai <- get_norm_data(x1)
            data_x2_smai <- get_norm_data(x2)

            sel1 <- seq_len(ncol(data_x1_smai))
            sel2 <- seq_len(ncol(data_x2_smai))

            out_smai <- align(
              data_x1_smai,
              data_x2_smai,
              sel1 = sel1,
              sel2 = sel2,
              t = t_iter_val[1],
              knn = knn_val[1],
              r.max = rmax_val[1],
              dir.map = "auto",
              denoise = "scree",
              prop.align = prop_align_val[1],
              cutoff = cutoff_val[1],
              cutoff.t = cutoff_t_val[1]
            )

            smai_data1 <- out_smai$data1.integrate
            smai_data2 <- out_smai$data2.integrate
            """
        )

        ref_int = _r_matrix_to_numpy(ro, "smai_data1")
        tar_int = _r_matrix_to_numpy(ro, "smai_data2")
        ro.r("rm(list = ls()); invisible(gc())")
        return ref_int, tar_int

    def apply(self, adata):
        print("\nApplying SMAI...")
        try:
            import rpy2.robjects as ro

            if self.integrate_key not in adata.obs:
                raise KeyError(
                    f"'{self.integrate_key}' not found in adata.obs."
                )

            batch_series = adata.obs[self.integrate_key].astype(str)
            batches = list(batch_series.unique())
            if len(batches) < 2:
                raise ValueError(
                    f"SMAI requires at least 2 batches, got {len(batches)}."
                )

            X, gene_names = _extract_expression_matrix(adata)
            feat_idx, feature_names = _select_common_features(
                X=X,
                gene_names=gene_names,
                nfeatures=int(self.nfeatures),
            )
            X = X[:, feat_idx]

            # Deterministic reference: largest batch; then lexicographic tie-break.
            batch_sizes = batch_series.value_counts()
            ref_batch = sorted(
                batch_sizes.index.tolist(),
                key=lambda b: (-int(batch_sizes[b]), str(b))
            )[0]
            ref_idx = np.where(batch_series.values == ref_batch)[0]

            embedding = np.zeros((adata.n_obs, len(feature_names)), dtype=np.float64)
            ref_accum = []

            for batch in batches:
                if batch == ref_batch:
                    continue

                tar_idx = np.where(batch_series.values == batch)[0]
                if len(tar_idx) == 0:
                    continue

                X_ref = X[ref_idx, :].T  # genes x cells
                X_tar = X[tar_idx, :].T

                print(f"  - SMAI aligning reference '{ref_batch}' vs '{batch}'...")
                ref_int, tar_int = self._run_smai_pair(
                    ro=ro,
                    X_ref=X_ref,
                    X_tar=X_tar,
                    gene_names=feature_names,
                )

                if ref_int.shape[1] != len(ref_idx) or tar_int.shape[1] != len(tar_idx):
                    raise RuntimeError(
                        f"SMAI output cell count mismatch for pair ({ref_batch}, {batch})."
                    )

                ref_accum.append(ref_int.T)
                embedding[tar_idx, :] = tar_int.T

            if len(ref_accum) == 0:
                raise RuntimeError("SMAI failed to align any non-reference batch.")

            # Average reference embedding across all pairwise alignments.
            embedding[ref_idx, :] = np.mean(np.stack(ref_accum, axis=0), axis=0)
            adata.obsm["X_smai"] = embedding

            # 尝试在维度一致时保存为表达层，便于后续分析。
            if embedding.shape[1] == adata.n_vars:
                adata.layers["smai_integrated"] = embedding.copy()

            sc.pp.neighbors(adata, use_rep="X_smai", random_state=42)
            sc.tl.leiden(adata, resolution=0.5, random_state=42)
            sc.tl.umap(adata, random_state=42)

            print("  ✓ SMAI completed")
            return adata
        except Exception as e:
            print(f"  ✗ SMAI failed: {e}")
            return None
