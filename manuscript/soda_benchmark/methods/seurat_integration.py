# -*- coding: utf-8 -*-
"""
Seurat CCA / RPCA 批次校正方法（通过 rpy2 调用 R）。

依赖
----
Python: rpy2 (>= 3.4)
R:      Seurat (>= 5.0.0)
"""

import numpy as np
import scanpy as sc
from scipy import sparse

from soda_benchmark.methods.base import BaseBatchMethod, register_method


# ------------------------------------------------------------------ #
# 依赖检查
# ------------------------------------------------------------------ #

def _check_rpy2_and_seurat():
    """检查 rpy2 和 R Seurat 包是否同时可用。"""
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

    return True, ""


# ------------------------------------------------------------------ #
# 表达矩阵提取
# ------------------------------------------------------------------ #

def _extract_expression_matrix(adata):
    """
    从 AnnData 提取表达矩阵，优先使用原始计数。

    Returns
    -------
    X : np.ndarray, shape (n_cells, n_genes), dtype float64
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

    if X is None and 'counts' in adata.layers:
        X = adata.layers['counts']

    if X is None:
        X = adata.X

    if sparse.issparse(X):
        X = np.asarray(X.todense(), dtype=np.float64)
    else:
        X = np.asarray(X, dtype=np.float64)

    return X, gene_names


# ------------------------------------------------------------------ #
# 核心整合流程
# ------------------------------------------------------------------ #

def _r_matrix_to_numpy(ro, r_varname):
    """从 R 全局环境取出矩阵并转为 numpy ndarray (cells × dims)。"""
    r_obj = ro.globalenv[r_varname]
    arr = np.array(r_obj, dtype=np.float64)
    if arr.ndim == 2:
        return arr
    dims = ro.r('dim')(r_obj)
    nrow, ncol = int(dims[0]), int(dims[1])
    return arr.reshape((nrow, ncol), order='F')


def _run_integration(adata, r_method_name, r_reduction_name, obsm_key):
    """
    通用 Seurat v5 整合：AnnData → R Seurat → IntegrateLayers → 嵌入回写 Python。

    Parameters
    ----------
    adata : AnnData
        输入数据（调用方已做 copy）。
    r_method_name : str
        R 端方法对象名，"CCAIntegration" 或 "RPCAIntegration"。
    r_reduction_name : str
        整合降维结果名，"integrated.cca" 或 "integrated.rpca"。
    obsm_key : str
        嵌入写入 adata.obsm 的键名。

    Returns
    -------
    AnnData
        含整合嵌入、leiden 聚类和 UMAP 的 AnnData。
    """
    import rpy2.robjects as ro

    X, gene_names = _extract_expression_matrix(adata)

    # R 矩阵格式: genes × cells
    X_t = np.ascontiguousarray(X.T, dtype=np.float64)
    n_genes, n_cells = X_t.shape

    cell_names = list(adata.obs_names)
    batch_labels = list(adata.obs['batch'].astype(str))

    # ---- 传输数据到 R（显式构造 R 对象，无需 converter） ----
    r_vec = ro.FloatVector(X_t.flatten(order='F'))
    ro.globalenv['expr_mat'] = ro.r.matrix(r_vec, nrow=n_genes, ncol=n_cells)
    ro.globalenv['gene_names_vec'] = ro.StrVector(gene_names)
    ro.globalenv['cell_names_vec'] = ro.StrVector(cell_names)
    ro.globalenv['batch_vec'] = ro.StrVector(batch_labels)

    # ---- Seurat 预处理 + 整合 ----
    ro.r(f'''
        suppressPackageStartupMessages(library(Seurat))
        options(future.globals.maxSize = 8 * 1024^3)

        rownames(expr_mat) <- gene_names_vec
        colnames(expr_mat) <- cell_names_vec

        obj <- CreateSeuratObject(counts = expr_mat)
        obj$batch <- batch_vec
        obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)

        obj <- NormalizeData(obj, verbose = FALSE)
        obj <- FindVariableFeatures(obj, verbose = FALSE)
        obj <- ScaleData(obj, verbose = FALSE)
        obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

        obj <- IntegrateLayers(
            object = obj,
            method = {r_method_name},
            orig.reduction = "pca",
            new.reduction = "{r_reduction_name}",
            verbose = FALSE
        )

        result_embeddings <- Embeddings(obj, reduction = "{r_reduction_name}")
    ''')

    embeddings = _r_matrix_to_numpy(ro, 'result_embeddings')

    # 清理 R 内存
    ro.r('rm(list = ls()); invisible(gc())')

    # ---- Python 后处理 ----
    adata.obsm[obsm_key] = embeddings
    sc.pp.neighbors(adata, use_rep=obsm_key, random_state=42)
    sc.tl.leiden(adata, resolution=0.5, random_state=42)
    sc.tl.umap(adata, random_state=42)

    return adata


# ------------------------------------------------------------------ #
# Seurat CCA
# ------------------------------------------------------------------ #

@register_method("Seurat_CCA")
class SeuratCCAMethod(BaseBatchMethod):
    """Seurat CCA (Canonical Correlation Analysis) 整合。"""

    def __init__(self, **kwargs):
        pass

    @classmethod
    def check_dependencies(cls):
        return _check_rpy2_and_seurat()

    def apply(self, adata):
        print("\nApplying Seurat CCA Integration...")
        try:
            result = _run_integration(
                adata,
                r_method_name="CCAIntegration",
                r_reduction_name="integrated.cca",
                obsm_key="X_seurat_cca",
            )
            print("  ✓ Seurat CCA completed")
            return result
        except Exception as e:
            print(f"  ✗ Seurat CCA failed: {e}")
            return None


# ------------------------------------------------------------------ #
# Seurat RPCA
# ------------------------------------------------------------------ #

@register_method("Seurat_RPCA")
class SeuratRPCAMethod(BaseBatchMethod):
    """Seurat RPCA (Reciprocal PCA) 整合。"""

    def __init__(self, **kwargs):
        pass

    @classmethod
    def check_dependencies(cls):
        return _check_rpy2_and_seurat()

    def apply(self, adata):
        print("\nApplying Seurat RPCA Integration...")
        try:
            result = _run_integration(
                adata,
                r_method_name="RPCAIntegration",
                r_reduction_name="integrated.rpca",
                obsm_key="X_seurat_rpca",
            )
            print("  ✓ Seurat RPCA completed")
            return result
        except Exception as e:
            print(f"  ✗ Seurat RPCA failed: {e}")
            return None
