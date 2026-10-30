# -*- coding: utf-8 -*-
"""
LIGER (iNMF) 批次校正方法（通过 rpy2 调用 R rliger）。

LIGER 使用 integrative Non-negative Matrix Factorization (iNMF) 对多批次数据
进行整合，再通过 quantile normalization 对齐各批次的细胞因子负载。

参考: https://welch-lab.github.io/liger/articles/liger_with_seurat.html

依赖
----
Python: rpy2 (>= 3.4)
R:      Seurat (>= 5.0.0), rliger (>= 2.0.0)
"""

import numpy as np
import scanpy as sc
from scipy import sparse

from soda_benchmark.methods.base import BaseBatchMethod, register_method


def _check_rpy2_seurat_liger():
    """检查 rpy2、R Seurat 和 rliger 是否同时可用。"""
    try:
        import rpy2.robjects as ro  # noqa: F401
    except ImportError:
        return False, "rpy2 not installed. Install with: pip install rpy2"

    try:
        ro.r('suppressPackageStartupMessages(library(Seurat))')
    except Exception:
        return False, (
            "R package 'Seurat' (>= 5.0) not available. "
            "Install in R with: install.packages('Seurat')"
        )

    try:
        ro.r('suppressPackageStartupMessages(library(rliger))')
    except Exception:
        return False, (
            "R package 'rliger' (>= 2.0) not available. "
            "Install in R with: install.packages('rliger')"
        )

    return True, ""


def _r_matrix_to_numpy(ro, r_varname):
    """从 R 全局环境取出矩阵并转为 numpy ndarray。"""
    r_obj = ro.globalenv[r_varname]
    arr = np.array(r_obj, dtype=np.float64)
    if arr.ndim == 2:
        return arr
    dims = ro.r('dim')(r_obj)
    nrow, ncol = int(dims[0]), int(dims[1])
    return arr.reshape((nrow, ncol), order='F')


def _extract_expression_matrix(adata):
    """从 AnnData 提取表达矩阵，优先使用原始计数。"""
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

    if X is None and 'counts' in adata.layers:
        X = adata.layers['counts']

    if X is None:
        X = adata.X

    if sparse.issparse(X):
        X = np.asarray(X.todense(), dtype=np.float64)
    else:
        X = np.asarray(X, dtype=np.float64)

    return X, gene_names


@register_method("LIGER")
class LigerMethod(BaseBatchMethod):
    """LIGER (iNMF + quantile normalization) 整合方法。"""

    def __init__(self, k=20, **kwargs):
        self.k = k

    @classmethod
    def check_dependencies(cls):
        return _check_rpy2_seurat_liger()

    def apply(self, adata):
        print("\nApplying LIGER (iNMF) Integration...")
        try:
            return self._run(adata)
        except Exception as e:
            print(f"  ✗ LIGER failed: {e}")
            return None

    def _run(self, adata):
        import rpy2.robjects as ro

        X, gene_names = _extract_expression_matrix(adata)

        X_t = np.ascontiguousarray(X.T, dtype=np.float64)
        n_genes, n_cells = X_t.shape

        cell_names = list(adata.obs_names)
        batch_labels = list(adata.obs['batch'].astype(str))

        r_vec = ro.FloatVector(X_t.flatten(order='F'))
        ro.globalenv['expr_mat'] = ro.r.matrix(r_vec, nrow=n_genes, ncol=n_cells)
        ro.globalenv['gene_names_vec'] = ro.StrVector(gene_names)
        ro.globalenv['cell_names_vec'] = ro.StrVector(cell_names)
        ro.globalenv['batch_vec'] = ro.StrVector(batch_labels)
        ro.globalenv['k_val'] = ro.IntVector([self.k])

        ro.r('''
            suppressPackageStartupMessages({
                library(Seurat)
                library(rliger)
            })
            options(future.globals.maxSize = 8 * 1024^3)

            rownames(expr_mat) <- gene_names_vec
            colnames(expr_mat) <- cell_names_vec

            obj <- CreateSeuratObject(counts = expr_mat)
            obj$batch <- batch_vec
            obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)

            # LIGER 专用预处理：normalize → selectGenes → scaleNotCenter
            obj <- normalize(obj)
            obj <- selectGenes(obj)
            obj <- scaleNotCenter(obj)

            # iNMF 整合 + quantile normalization
            obj <- runINMF(obj, k = k_val[1])
            obj <- quantileNorm(obj)

            liger_emb <- Embeddings(obj, reduction = "inmfNorm")
        ''')

        embeddings = _r_matrix_to_numpy(ro, 'liger_emb')
        ro.r('rm(list = ls()); invisible(gc())')

        adata.obsm['X_liger'] = embeddings
        sc.pp.neighbors(adata, use_rep='X_liger', random_state=42)
        sc.tl.leiden(adata, resolution=0.5, random_state=42)
        sc.tl.umap(adata, random_state=42)

        print("  ✓ LIGER completed")
        return adata
