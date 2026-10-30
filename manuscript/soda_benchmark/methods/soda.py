# -*- coding: utf-8 -*-
"""
SODA: Sinkhorn Optimal Transport with Divergence-Based Adaptation

包含 SODABatchCorrector 核心算法与对应的方法注册。
"""

import numpy as np
import scanpy as sc
from sklearn.cluster import KMeans

from soda_benchmark.methods.base import BaseBatchMethod, register_method

try:
    import ot
    HAS_OT = True
except ImportError:
    HAS_OT = False
    print("Warning: POT library not installed. SODA will be skipped.")
    print("Install with: pip install POT")


# ------------------------------------------------------------------ #
# SODA 核心算法
# ------------------------------------------------------------------ #

class SODABatchCorrector:
    """
    SODA: Sinkhorn Optimal Transport with Divergence-Based Adaptation

    Harmony式迭代聚类 + Sinkhorn 散度驱动的 OT 对齐。
    """

    def __init__(self, n_clusters=50, max_iter=10, epsilon=0.05,
                 divergence_threshold=0.05, alpha=0.7):
        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.epsilon = epsilon
        self.divergence_threshold = divergence_threshold
        if not 0.0 <= alpha <= 1.0:
            raise ValueError("alpha must be in [0, 1]")
        self.alpha = alpha
        self.correction_history = []

    def compute_sinkhorn_divergence(self, X_i, X_j):
        """计算 Sinkhorn 散度，用于衡量批次差异。"""
        try:
            M_ij = ot.dist(X_i, X_j, metric='sqeuclidean')
            a = np.ones(len(X_i)) / len(X_i)
            b = np.ones(len(X_j)) / len(X_j)
            W2_ij = ot.sinkhorn2(a, b, M_ij, self.epsilon)

            M_ii = ot.dist(X_i, X_i, metric='sqeuclidean')
            M_jj = ot.dist(X_j, X_j, metric='sqeuclidean')
            W2_ii = ot.sinkhorn2(a, a, M_ii, self.epsilon)
            W2_jj = ot.sinkhorn2(b, b, M_jj, self.epsilon)

            result = W2_ij - 0.5 * W2_ii - 0.5 * W2_jj
            return max(0.0, float(result))  # 避免数值精度导致的负值
        except Exception:
            return 0.0

    def detect_global_batch_effect(self, adata, batch_key='batch', use_rep='X_pca'):
        """通过跨批次 Sinkhorn 散度估计整体批次效应。"""
        X = adata.obsm[use_rep] if use_rep in adata.obsm else (
            adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X)

        divergences = []
        batches = adata.obs[batch_key].unique()
        if len(batches) < 2:
            return 0.0

        for i, batch_i in enumerate(batches):
            for j, batch_j in enumerate(batches):
                if i >= j:
                    continue
                mask_i = adata.obs[batch_key] == batch_i
                mask_j = adata.obs[batch_key] == batch_j
                X_i, X_j = X[mask_i], X[mask_j]

                if len(X_i) < 5 or len(X_j) < 5:
                    continue

                sinkhorn_div = self.compute_sinkhorn_divergence(X_i, X_j)
                if sinkhorn_div > -1e-10:  # 包含数值零
                    divergences.append(max(1e-12, sinkhorn_div))

        return np.mean(divergences) if divergences else 0.0

    def compute_cluster_divergence(self, X_cluster, batches):
        """在聚类内部计算跨批次 Sinkhorn 散度，用于自适应校正。"""
        unique_batches = batches.unique()
        if len(unique_batches) < 2:
            return 0.0

        divergences = []
        for i, batch_i in enumerate(unique_batches):
            for j, batch_j in enumerate(unique_batches):
                if i >= j:
                    continue
                X_i = X_cluster[batches == batch_i]
                X_j = X_cluster[batches == batch_j]

                if len(X_i) < 5 or len(X_j) < 5:
                    continue

                sinkhorn_div = self.compute_sinkhorn_divergence(X_i, X_j)
                if sinkhorn_div > -1e-10:  # 包含数值零
                    divergences.append(max(1e-12, sinkhorn_div))

        raw = np.mean(divergences) if divergences else 0.0
        # 当有两批次但散度为 0（数值精度）时，使用阈值触发校正
        if len(unique_batches) >= 2 and raw < 1e-10:
            return max(0.01, self.divergence_threshold)
        return raw

    def correct(self, adata, batch_key='batch', use_rep='X_pca'):
        """执行 SODA 校正流程。"""
        if use_rep not in adata.obsm:
            raise ValueError(f"{use_rep} not found in adata.obsm")

        X = np.asarray(adata.obsm[use_rep].copy(), dtype=np.float64)
        if np.any(np.isnan(X)):
            X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

        for iteration in range(self.max_iter):
            if np.any(np.isnan(X)) or np.any(np.isinf(X)):
                X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
            kmeans = KMeans(n_clusters=self.n_clusters, random_state=42, n_init=10)
            clusters = kmeans.fit_predict(X)

            for cluster_id in range(self.n_clusters):
                mask_cluster = clusters == cluster_id
                if mask_cluster.sum() < 10:
                    continue

                cluster_divergence = self.compute_cluster_divergence(
                    X[mask_cluster], adata.obs.loc[mask_cluster, batch_key])

                if cluster_divergence < self.divergence_threshold:
                    continue

                X_cluster = X[mask_cluster]
                batches_cluster = adata.obs.loc[mask_cluster, batch_key]
                ref_batch = batches_cluster.value_counts().idxmax()
                mask_ref = (batches_cluster == ref_batch).values
                X_ref = X_cluster[mask_ref]

                for batch in batches_cluster.unique():
                    if batch == ref_batch:
                        continue
                    mask_batch = (batches_cluster == batch).values
                    X_batch = X_cluster[mask_batch]
                    if len(X_batch) == 0 or len(X_ref) == 0:
                        continue
                    try:
                        M = ot.dist(X_batch, X_ref, metric='sqeuclidean')
                        a = np.ones(len(X_batch)) / len(X_batch)
                        b = np.ones(len(X_ref)) / len(X_ref)
                        T = ot.sinkhorn(a, b, M, reg=self.epsilon)
                        T_sum = np.maximum(T.sum(axis=1, keepdims=True), 1e-12)
                        X_corrected_batch = (T / T_sum) @ X_ref

                        batch_indices = np.where(mask_cluster)[0][mask_batch]
                        corrected = (
                            self.alpha * X_corrected_batch
                            + (1.0 - self.alpha) * X_batch
                        )
                        bad = np.isnan(corrected) | np.isinf(corrected)
                        corrected[bad] = X_batch[bad]
                        X[batch_indices] = corrected
                    except Exception:
                        continue

            if iteration >= 3:
                adata_temp = adata.copy()
                adata_temp.obsm[use_rep] = X
                global_divergence = self.detect_global_batch_effect(
                    adata_temp, batch_key=batch_key, use_rep=use_rep)
                if global_divergence < self.divergence_threshold:
                    break

        return X


# ------------------------------------------------------------------ #
# 注册 SODA 方法
# ------------------------------------------------------------------ #

@register_method("SODA")
class SODAMethod(BaseBatchMethod):

    def __init__(self, n_clusters=50, max_iter=10, epsilon=0.05,
                 divergence_threshold=0.05, tau=None, alpha=0.7, **kwargs):
        # tau 是论文中的符号，这里作为 divergence_threshold 的别名暴露。
        if tau is not None:
            divergence_threshold = tau
        self.soda_params = dict(
            n_clusters=n_clusters,
            max_iter=max_iter,
            epsilon=epsilon,
            divergence_threshold=divergence_threshold,
            alpha=alpha,
        )

    @classmethod
    def check_dependencies(cls):
        if HAS_OT:
            return True, ""
        return False, "POT library not installed. Install with: pip install POT"

    def apply(self, adata):
        print("\nApplying SODA...")
        try:
            if 'X_pca' not in adata.obsm:
                sc.tl.pca(adata, n_comps=50, random_state=42)

            corrector = SODABatchCorrector(**self.soda_params)
            X_corrected = corrector.correct(adata, batch_key='batch', use_rep='X_pca')
            adata.obsm['X_soda'] = X_corrected

            sc.pp.neighbors(adata, use_rep='X_soda', random_state=42)
            sc.tl.leiden(adata, resolution=0.5, random_state=42)
            sc.tl.umap(adata, random_state=42)
            print("  ✓ SODA completed")
            return adata
        except Exception as e:
            print(f"  ✗ SODA failed: {e}")
            return None

