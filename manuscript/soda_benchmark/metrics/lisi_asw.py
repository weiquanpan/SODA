# -*- coding: utf-8 -*-
"""LISI 与 ASW(batch) 指标实现。

参考 scib-metrics 的常见定义：
- iLISI: 批次混合度，越大越好（此处归一化到 [0, 1]）
- cLISI: 细胞类型保真，越大越好（由原始 cLISI 反向归一化到 [0, 1]）
- ASW_batch: 按细胞类型分组后的 batch silhouette，越大越好
"""

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import silhouette_samples
from sklearn.neighbors import NearestNeighbors

from soda_benchmark.metrics.base import BaseMetric, MetricResult, register_metric


class _EmbeddingMetricMixin:
    """提供 embedding 选择与标签解析的公共工具。"""

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
    _LABEL_CANDIDATES = ["cell_type", "celltype", "CellType", "leiden", "label"]

    def __init__(self, neighborhood_mode="representation", label_key=None,
                 lisi_n_neighbors=30, **kwargs):
        super().__init__(**kwargs)
        mode = str(neighborhood_mode).lower().strip()
        if mode not in {"expression", "representation"}:
            raise ValueError("neighborhood_mode must be 'expression' or 'representation'.")
        self.neighborhood_mode = mode
        self.label_key = label_key
        self.lisi_n_neighbors = int(lisi_n_neighbors)

    @staticmethod
    def _to_dense(X):
        if sparse.issparse(X):
            return X.toarray()
        return np.asarray(X)

    def _get_matrix(self, adata, method_name):
        if self.neighborhood_mode == "expression":
            return self._to_dense(adata.X), "X"

        candidates = list(self._METHOD_REP_KEYS.get(method_name, []))
        generic_key = f"X_{method_name.lower()}"
        if generic_key not in candidates:
            candidates.append(generic_key)

        for key in candidates:
            if key in adata.obsm and adata.obsm[key].shape[0] == adata.n_obs:
                return np.asarray(adata.obsm[key]), key

        return self._to_dense(adata.X), "X(fallback)"

    def _resolve_label_key(self, adata):
        if self.label_key:
            if self.label_key not in adata.obs.columns:
                raise KeyError(
                    f"label_key '{self.label_key}' not found in adata.obs. "
                    f"Available: {list(adata.obs.columns)}"
                )
            return self.label_key
        for key in self._LABEL_CANDIDATES:
            if key in adata.obs.columns:
                return key
        raise KeyError(
            "No label column found for cLISI/ASW_batch. "
            f"Tried: {self._LABEL_CANDIDATES}"
        )

    def _knn_indices(self, X):
        n_cells = X.shape[0]
        # 至少为 2（含自身 + 1 个邻居），最大不超过样本数
        n_neighbors = max(2, min(self.lisi_n_neighbors + 1, n_cells))
        nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
        nn.fit(X)
        idx = nn.kneighbors(X, return_distance=False)
        # 去掉自身（第一列）
        return idx[:, 1:]

    @staticmethod
    def _compute_lisi_values(knn_idx, label_codes, n_categories):
        n_cells, k = knn_idx.shape
        values = np.zeros(n_cells, dtype=float)
        for i in range(n_cells):
            neigh_codes = label_codes[knn_idx[i]]
            counts = np.bincount(neigh_codes, minlength=n_categories).astype(float)
            probs = counts / k
            simpson = np.sum(probs * probs)
            values[i] = 1.0 / simpson if simpson > 0 else 1.0
        return values


@register_metric("ilisi")
class ILISIMetric(_EmbeddingMetricMixin, BaseMetric):
    """iLISI（批次混合）指标。"""

    name = "iLISI"
    higher_is_better = True

    def calculate(self, adata_original, adata_dict) -> MetricResult:
        print(f"\nCalculating iLISI... (mode={self.neighborhood_mode})")
        rows = {}
        per_cell = {}
        matrix_source = {}

        for method_name, adata in adata_dict.items():
            X, source_key = self._get_matrix(adata, method_name)
            matrix_source[method_name] = source_key
            knn_idx = self._knn_indices(X)

            batch_codes, batch_uniques = pd.factorize(adata.obs["batch"].astype(str), sort=False)
            n_batches = len(batch_uniques)
            raw_ilisi = self._compute_lisi_values(knn_idx, batch_codes, n_batches)
            if n_batches <= 1:
                ilisi_norm = np.ones_like(raw_ilisi)
            else:
                ilisi_norm = (raw_ilisi - 1.0) / (n_batches - 1.0)
                ilisi_norm = np.clip(ilisi_norm, 0.0, 1.0)

            per_cell[method_name] = ilisi_norm.tolist()
            score = float(np.mean(ilisi_norm))
            ci_lower, ci_upper = self._bootstrap_ci(ilisi_norm, reducer=np.mean)
            rows[method_name] = {
                "score": score,
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
            }

        scores_df = pd.DataFrame.from_dict(rows, orient="index")
        print("  ✓ iLISI calculated")
        return MetricResult(
            name=self.name,
            key="ilisi",
            higher_is_better=self.higher_is_better,
            scores=scores_df,
            extra={
                "per_cell": pd.DataFrame(per_cell),
                "mode": self.neighborhood_mode,
                "matrix_source": matrix_source,
            },
        )

    def grade(self, value: float) -> str:
        if value > 0.8:
            return "Excellent"
        if value > 0.6:
            return "Good"
        if value > 0.4:
            return "Fair"
        return "Poor"


@register_metric("clisi")
class CLISIMetric(_EmbeddingMetricMixin, BaseMetric):
    """cLISI（细胞类型保真）指标。"""

    name = "cLISI"
    higher_is_better = True

    def calculate(self, adata_original, adata_dict) -> MetricResult:
        label_key = self._resolve_label_key(adata_original)
        print(
            f"\nCalculating cLISI... (mode={self.neighborhood_mode}, label_key={label_key})"
        )
        rows = {}
        per_cell = {}
        matrix_source = {}

        for method_name, adata in adata_dict.items():
            X, source_key = self._get_matrix(adata, method_name)
            matrix_source[method_name] = source_key
            knn_idx = self._knn_indices(X)

            label_codes, label_uniques = pd.factorize(
                adata.obs[label_key].astype(str), sort=False
            )
            n_labels = len(label_uniques)
            raw_clisi = self._compute_lisi_values(knn_idx, label_codes, n_labels)
            if n_labels <= 1:
                clisi_norm = np.ones_like(raw_clisi)
            else:
                # 原始 cLISI 越小越好，这里转成 [0,1] 的“越大越好”
                clisi_norm = (n_labels - raw_clisi) / (n_labels - 1.0)
                clisi_norm = np.clip(clisi_norm, 0.0, 1.0)

            per_cell[method_name] = clisi_norm.tolist()
            score = float(np.mean(clisi_norm))
            ci_lower, ci_upper = self._bootstrap_ci(clisi_norm, reducer=np.mean)
            rows[method_name] = {
                "score": score,
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
            }

        scores_df = pd.DataFrame.from_dict(rows, orient="index")
        print("  ✓ cLISI calculated")
        return MetricResult(
            name=self.name,
            key="clisi",
            higher_is_better=self.higher_is_better,
            scores=scores_df,
            extra={
                "per_cell": pd.DataFrame(per_cell),
                "mode": self.neighborhood_mode,
                "label_key": label_key,
                "matrix_source": matrix_source,
            },
        )

    def grade(self, value: float) -> str:
        if value > 0.8:
            return "Excellent"
        if value > 0.6:
            return "Good"
        if value > 0.4:
            return "Fair"
        return "Poor"


@register_metric("asw_batch")
class ASWBatchMetric(_EmbeddingMetricMixin, BaseMetric):
    """ASW_batch：按细胞类型分组后衡量批次混合。"""

    name = "ASW batch"
    higher_is_better = True

    def calculate(self, adata_original, adata_dict) -> MetricResult:
        label_key = self._resolve_label_key(adata_original)
        print(
            f"\nCalculating ASW_batch... (mode={self.neighborhood_mode}, label_key={label_key})"
        )
        rows = {}
        per_cell = {}
        matrix_source = {}

        for method_name, adata in adata_dict.items():
            X, source_key = self._get_matrix(adata, method_name)
            matrix_source[method_name] = source_key
            labels = adata.obs[label_key].astype(str).values
            batches = adata.obs["batch"].astype(str).values

            cell_scores = np.full(adata.n_obs, np.nan, dtype=float)
            for lab in np.unique(labels):
                idx = np.where(labels == lab)[0]
                if len(idx) < 3:
                    continue
                sub_batches = batches[idx]
                if len(np.unique(sub_batches)) < 2:
                    continue
                sub_codes, _ = pd.factorize(sub_batches, sort=False)
                try:
                    sil = silhouette_samples(X[idx], sub_codes, metric="euclidean")
                except Exception:
                    continue
                # 与 scIB 常用定义一致：1 - |silhouette|
                cell_scores[idx] = 1.0 - np.abs(sil)

            valid = cell_scores[~np.isnan(cell_scores)]
            if valid.size == 0:
                score = 0.0
                ci_lower, ci_upper = 0.0, 0.0
                per_cell[method_name] = [None] * adata.n_obs
            else:
                score = float(np.mean(valid))
                ci_lower, ci_upper = self._bootstrap_ci(valid, reducer=np.mean)
                per_cell[method_name] = cell_scores.tolist()

            rows[method_name] = {
                "score": score,
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
            }

        scores_df = pd.DataFrame.from_dict(rows, orient="index")
        print("  ✓ ASW_batch calculated")
        return MetricResult(
            name=self.name,
            key="asw_batch",
            higher_is_better=self.higher_is_better,
            scores=scores_df,
            extra={
                "per_cell": pd.DataFrame(per_cell),
                "mode": self.neighborhood_mode,
                "label_key": label_key,
                "matrix_source": matrix_source,
            },
        )

    def grade(self, value: float) -> str:
        if value > 0.8:
            return "Excellent"
        if value > 0.6:
            return "Good"
        if value > 0.4:
            return "Fair"
        return "Poor"
