# -*- coding: utf-8 -*-
"""调整兰德指数（Adjusted Rand Index）。"""

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from soda_benchmark.metrics.base import BaseMetric, MetricResult, register_metric


_ARI_REF_CANDIDATES = ("cell_type", "CellType", "celltype", "cell_ontology_class")

@register_metric("ari")
class ARIMetric(BaseMetric):
    """衡量校正后聚类与参考标签的一致性（ARI）。优先使用真实细胞类型作为参考。"""

    name = "Adjusted Rand Index"
    higher_is_better = True  # 越大越好

    def __init__(self, label_key=None, **kwargs):
        super().__init__(**kwargs)
        self.label_key = label_key

    def _resolve_ari_reference(self, adata):
        """解析 ARI 参考标签：优先 cell_type，否则 leiden。"""
        if self.label_key and self.label_key in adata.obs.columns:
            return adata.obs[self.label_key].astype(str).values
        for key in _ARI_REF_CANDIDATES:
            if key in adata.obs.columns:
                return adata.obs[key].astype(str).values
        return adata.obs['leiden'].astype(str).values

    def calculate(self, adata_original, adata_dict) -> MetricResult:
        print("\nCalculating ARI...")

        ref_labels = self._resolve_ari_reference(adata_original)
        rows = {}

        for method_name, adata in adata_dict.items():
            pred_labels = adata.obs['leiden'].astype(str).values
            score = adjusted_rand_score(ref_labels, pred_labels)
            if method_name == 'Original':
                ci_lower, ci_upper = float(score), float(score)
            else:
                ci_lower, ci_upper = self._bootstrap_ari(ref_labels, pred_labels)
            rows[method_name] = {
                'score': score,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
            }

        scores_df = pd.DataFrame.from_dict(rows, orient='index')
        print("  ✓ ARI calculated")
        return MetricResult(
            name=self.name,
            key="ari",
            higher_is_better=self.higher_is_better,
            scores=scores_df,
        )

    def grade(self, value: float) -> str:
        if value > 0.8:
            return "Excellent"
        if value > 0.5:
            return "Good"
        if value > 0.3:
            return "Fair"
        return "Poor"

    # ------------------------------------------------------------------ #
    # 内部计算
    # ------------------------------------------------------------------ #

    def _bootstrap_ari(self, orig, corrected):
        stats = []
        n = len(orig)
        for _ in range(self.n_bootstrap):
            idx = self.rng.integers(0, n, size=n)
            stats.append(adjusted_rand_score(orig[idx], corrected[idx]))
        return self._ci_bounds(np.array(stats))

