# -*- coding: utf-8 -*-
"""Harmony 批次校正方法。"""

import numpy as np
import pandas as pd
import scanpy as sc

from soda_benchmark.methods.base import BaseBatchMethod, register_method

try:
    import harmonypy as _hm
    _HAS_HARMONYPY = True
except ImportError:
    _HAS_HARMONYPY = False


@register_method("Harmony")
class HarmonyMethod(BaseBatchMethod):

    def __init__(self, **kwargs):
        pass

    @classmethod
    def check_dependencies(cls):
        # 优先检查 harmonypy 是否可用；其次检查 scanpy 包装器
        if _HAS_HARMONYPY:
            return True, ""
        if hasattr(sc.external.pp, 'harmony_integrate'):
            return True, ""
        return False, "harmony not available. Install with: pip install harmonypy"

    def apply(self, adata):
        print("\nApplying Harmony...")
        try:
            # 确保 X_pca 存在且维度正确 (n_obs, n_comps)
            need_pca = True
            if 'X_pca' in adata.obsm:
                xpca = adata.obsm['X_pca']
                if xpca.ndim == 2 and xpca.shape[0] == adata.n_obs:
                    need_pca = False
            if need_pca:
                sc.tl.pca(adata, n_comps=50, random_state=42)

            # 确保 X_pca 为 float64
            adata.obsm['X_pca'] = np.asarray(
                adata.obsm['X_pca'], dtype=np.float64)

            if not isinstance(adata.obs['batch'].dtype, pd.CategoricalDtype):
                adata.obs['batch'] = (
                    adata.obs['batch'].astype(str).astype('category')
                )

            # 直接调用 harmonypy，绕过 scanpy 包装器
            # scanpy 的 harmony_integrate 对 Z_corr 做了 .T 转置，
            # 但新版 harmonypy 已将输出改为 (n_cells, n_components)，
            # 转置后反而变成 (n_components, n_cells) 导致 obsm 赋值失败。
            if _HAS_HARMONYPY:
                ho = _hm.run_harmony(
                    adata.obsm['X_pca'],
                    adata.obs,
                    'batch',
                    max_iter_harmony=10,
                )
                Z_corr = ho.Z_corr
                # PyTorch CUDA tensor → CPU numpy
                if hasattr(Z_corr, 'cpu'):
                    Z_corr = Z_corr.cpu().detach().numpy()
                Z_corr = np.asarray(Z_corr, dtype=np.float64)
                # 兼容新旧版 harmonypy：
                #   新版 Z_corr shape = (n_cells, n_components)
                #   旧版 Z_corr shape = (n_components, n_cells)
                if Z_corr.shape[0] == adata.n_obs:
                    adata.obsm['X_harmony'] = Z_corr
                else:
                    adata.obsm['X_harmony'] = Z_corr.T
            else:
                sc.external.pp.harmony_integrate(
                    adata,
                    key='batch',
                    basis='X_pca',
                    adjusted_basis='X_harmony',
                    max_iter_harmony=10,
                )

            sc.pp.neighbors(adata, use_rep='X_harmony', random_state=42)
            sc.tl.leiden(adata, resolution=0.5, random_state=42)
            sc.tl.umap(adata, random_state=42)
            print("  ✓ Harmony completed")
            return adata
        except Exception as e:
            print(f"  ✗ Harmony failed: {e}")
            return None
