# -*- coding: utf-8 -*-
"""邻域保持指标（Neighborhood Preservation）。

支持两种比较口径：
1) expression：在 `adata.X` 空间比较邻域（默认）
2) representation：在方法对应 embedding 空间比较邻域（优先使用 `obsm`）
"""

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

from soda_benchmark.metrics.base import BaseMetric, MetricResult, register_metric


@register_metric("neighborhood_preservation")
class NeighborhoodPreservation(BaseMetric):
    """衡量校正前后 k 近邻排名变化的中位数。"""

    name = "Neighborhood Preservation"
    higher_is_better = False  # 越小越好

    _METHOD_REP_KEYS = {
        "Original": ["X_pca"],
        "SODA": ["X_soda"],
        "scVI": ["X_scVI"],
        "Harmony": ["X_harmony"],
        "LIGER": ["X_liger"],
        "Seurat_CCA": ["X_seurat_cca"],
        "Seurat_RPCA": ["X_seurat_rpca"],
        "CellANOVA": ["X_cellanova"],
        "SMAI": ["X_smai"],
    }

    def __init__(self, k: int = 30, neighborhood_mode: str = "expression", **kwargs):
        super().__init__(**kwargs)
        self.k = k
        mode = str(neighborhood_mode).lower().strip()
        if mode not in {"expression", "representation"}:
            raise ValueError(
                "neighborhood_mode must be 'expression' or 'representation'."
            )
        self.neighborhood_mode = mode

    def _to_dense(self, X):
        if sparse.issparse(X):
            return X.toarray()
        return np.asarray(X)

    def _get_matrix(self, adata, method_name):
        if self.neighborhood_mode == "expression":
            return self._to_dense(adata.X), "X"

        # representation mode: 优先使用已知方法对应 embedding
        candidates = list(self._METHOD_REP_KEYS.get(method_name, []))
        generic_key = f"X_{method_name.lower()}"
        if generic_key not in candidates:
            candidates.append(generic_key)

        for key in candidates:
            if key in adata.obsm and adata.obsm[key].shape[0] == adata.n_obs:
                return np.asarray(adata.obsm[key]), key

        # 未找到 embedding 时，回退到表达空间，避免方法被跳过
        return self._to_dense(adata.X), "X(fallback)"

    def calculate(self, adata_original, adata_dict) -> MetricResult:
        print(f"\nCalculating neighborhood preservation... (mode={self.neighborhood_mode})")

        X_orig, orig_source = self._get_matrix(adata_original, "Original")
        n_cells = X_orig.shape[0]
        n_neighbors = max(2, min(self.k + 1, n_cells))
        orig_nn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
        orig_nn.fit(X_orig)
        top_k_orig = orig_nn.kneighbors(X_orig, return_distance=False)[:, 1:]  # (n, k_eff)
        k_eff = top_k_orig.shape[1]
        orig_ranks = np.arange(1, k_eff + 1)[None, :]  # (1, k_eff)

        per_cell_results = {}
        summary_rows = {}
        matrix_source = {"Original": orig_source}

        for method_name, adata in adata_dict.items():
            X_corrected, source_key = self._get_matrix(adata, method_name)
            matrix_source[method_name] = source_key
            corr_nn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
            corr_nn.fit(X_corrected)
            top_k_corr = corr_nn.kneighbors(X_corrected, return_distance=False)[:, 1:]

            # 内存友好近似：仅在 corrected 的 top-k 邻域内给出排名；
            # 若原邻居不在 corrected top-k，赋予 k+1 惩罚排名。
            new_ranks = np.full((n_cells, k_eff), k_eff + 1, dtype=np.int32)
            for i in range(n_cells):
                corr_rank_map = {idx: rank for rank, idx in enumerate(top_k_corr[i], start=1)}
                for j, neigh in enumerate(top_k_orig[i]):
                    new_ranks[i, j] = corr_rank_map.get(int(neigh), k_eff + 1)

            rank_changes = np.median(np.abs(new_ranks - orig_ranks), axis=1)

            per_cell_results[method_name] = rank_changes.tolist()
            median_value = float(np.median(rank_changes))
            ci_lower, ci_upper = self._bootstrap_ci(
                rank_changes, reducer=np.median)
            summary_rows[method_name] = {
                'score': median_value,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
            }

        scores_df = pd.DataFrame.from_dict(summary_rows, orient='index')
        per_cell_df = pd.DataFrame(per_cell_results)

        print("  ✓ Neighborhood preservation calculated")
        return MetricResult(
            name=self.name,
            key="neighborhood_preservation",
            higher_is_better=self.higher_is_better,
            scores=scores_df,
            extra={
                'per_cell': per_cell_df,
                'mode': self.neighborhood_mode,
                'matrix_source': matrix_source,
            },
        )

    def grade(self, value: float) -> str:
        if value < 10:
            return "Excellent"
        if value < 30:
            return "Good"
        if value < 50:
            return "Fair"
        return "Poor"
