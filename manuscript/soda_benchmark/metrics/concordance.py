# -*- coding: utf-8 -*-
"""聚类一致性指标（Cluster Concordance）。"""

import numpy as np
import pandas as pd

from soda_benchmark.metrics.base import BaseMetric, MetricResult, register_metric


@register_metric("cluster_concordance")
class ClusterConcordance(BaseMetric):
    """衡量校正前后聚类标签的匹配程度。"""

    name = "Cluster Concordance"
    higher_is_better = False  # 越小越好

    def calculate(self, adata_original, adata_dict) -> MetricResult:
        print("\nCalculating cluster concordance...")

        orig_clusters = np.asarray(adata_original.obs['leiden'].values)
        rows = {}

        for method_name, adata in adata_dict.items():
            corrected_clusters = np.asarray(adata.obs['leiden'].values)
            score = self._cluster_concordance_value(orig_clusters, corrected_clusters)
            if method_name == 'Original':
                # Baseline against itself is deterministic.
                ci_lower, ci_upper = float(score), float(score)
            else:
                ci_lower, ci_upper = self._bootstrap_cluster_concordance(
                    orig_clusters, corrected_clusters)
            rows[method_name] = {
                'score': score,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
            }

        scores_df = pd.DataFrame.from_dict(rows, orient='index')
        print("  ✓ Cluster concordance calculated")
        return MetricResult(
            name=self.name,
            key="cluster_concordance",
            higher_is_better=self.higher_is_better,
            scores=scores_df,
        )

    def grade(self, value: float) -> str:
        if value < 0.1:
            return "Excellent"
        if value < 0.3:
            return "Good"
        if value < 0.5:
            return "Fair"
        return "Poor"

    # ------------------------------------------------------------------ #
    # 内部计算
    # ------------------------------------------------------------------ #

    def _reindex_labels(self, labels):
        """
        将任意标签重映射为连续整数 0..K-1，避免非连续标签索引越界。
        """
        codes, _ = pd.factorize(labels, sort=False)
        return np.asarray(codes, dtype=np.int64)

    def _cluster_concordance_value(self, orig_clusters, corrected_clusters):
        orig_idx = self._reindex_labels(orig_clusters)
        corr_idx = self._reindex_labels(corrected_clusters)
        if len(orig_idx) != len(corr_idx):
            raise ValueError("orig_clusters and corrected_clusters must have same length.")

        n_orig = len(np.unique(orig_idx))
        n_corr = len(np.unique(corr_idx))

        confusion = np.zeros((n_orig, n_corr))
        for i in range(len(orig_idx)):
            confusion[orig_idx[i], corr_idx[i]] += 1

        cluster_changes = []
        for i in range(n_orig):
            if confusion[i].sum() > 0:
                cluster_changes.append(
                    (confusion[i].sum() - confusion[i].max()) / confusion[i].sum()
                )

        cluster_changes_reverse = []
        for j in range(n_corr):
            if confusion[:, j].sum() > 0:
                cluster_changes_reverse.append(
                    (confusion[:, j].sum() - confusion[:, j].max()) / confusion[:, j].sum()
                )

        forward = np.mean(cluster_changes) if cluster_changes else 0.0
        backward = np.mean(cluster_changes_reverse) if cluster_changes_reverse else 0.0
        return (forward + backward) / 2

    def _bootstrap_cluster_concordance(self, orig, corrected):
        stats = []
        n = len(orig)
        for _ in range(self.n_bootstrap):
            idx = self.rng.integers(0, n, size=n)
            stats.append(self._cluster_concordance_value(orig[idx], corrected[idx]))
        return self._ci_bounds(np.array(stats))

