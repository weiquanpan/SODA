"""Core SODA algorithm implementation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.cluster import KMeans

try:
    import ot
except ImportError:  # pragma: no cover - handled at runtime
    ot = None


@dataclass(slots=True)
class SODAResult:
    """Result bundle returned by the SODA integrator."""

    embedding: np.ndarray
    n_iterations: int
    correction_history: list[float]


class SODA:
    """Sinkhorn Optimal Transport with Divergence-Based Adaptation."""

    def __init__(
        self,
        n_clusters: int = 50,
        max_iter: int = 10,
        epsilon: float = 0.05,
        tau: float = 0.05,
        alpha: float = 0.7,
        random_state: int = 42,
    ) -> None:
        if n_clusters < 2:
            raise ValueError("n_clusters must be at least 2.")
        if max_iter < 1:
            raise ValueError("max_iter must be at least 1.")
        if epsilon <= 0:
            raise ValueError("epsilon must be positive.")
        if tau < 0:
            raise ValueError("tau must be non-negative.")
        if not 0.0 <= alpha <= 1.0:
            raise ValueError("alpha must be in [0, 1].")

        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.epsilon = epsilon
        self.tau = tau
        self.alpha = alpha
        self.random_state = random_state

    def fit_transform(self, embedding: np.ndarray, batches: np.ndarray) -> np.ndarray:
        """Return the corrected embedding only."""
        return self.integrate(embedding, batches).embedding

    def integrate(self, embedding: np.ndarray, batches: np.ndarray) -> SODAResult:
        """Run SODA on a precomputed embedding and batch labels."""
        self._require_pot()

        x = self._validate_embedding(embedding)
        batch_labels = self._validate_batches(batches, n_obs=x.shape[0])

        if x.shape[0] < 2 or np.unique(batch_labels).size < 2:
            return SODAResult(
                embedding=x.copy(),
                n_iterations=0,
                correction_history=[],
            )

        corrected = x.copy()
        correction_history: list[float] = []
        n_clusters = min(self.n_clusters, corrected.shape[0])

        for iteration in range(1, self.max_iter + 1):
            corrected = self._sanitize(corrected)
            kmeans = KMeans(
                n_clusters=n_clusters,
                random_state=self.random_state,
                n_init=10,
            )
            clusters = kmeans.fit_predict(corrected)

            for cluster_id in range(n_clusters):
                cluster_mask = clusters == cluster_id
                if int(cluster_mask.sum()) < 10:
                    continue

                cluster_embedding = corrected[cluster_mask]
                cluster_batches = batch_labels[cluster_mask]
                cluster_divergence = self._cluster_divergence(
                    cluster_embedding,
                    cluster_batches,
                )
                if cluster_divergence < self.tau:
                    continue

                corrected[cluster_mask] = self._correct_cluster(
                    cluster_embedding,
                    cluster_batches,
                )

            global_divergence = self._detect_global_batch_effect(corrected, batch_labels)
            correction_history.append(global_divergence)
            if iteration >= 3 and global_divergence < self.tau:
                return SODAResult(
                    embedding=corrected,
                    n_iterations=iteration,
                    correction_history=correction_history,
                )

        return SODAResult(
            embedding=corrected,
            n_iterations=self.max_iter,
            correction_history=correction_history,
        )

    def _correct_cluster(
        self,
        cluster_embedding: np.ndarray,
        cluster_batches: np.ndarray,
    ) -> np.ndarray:
        corrected_cluster = cluster_embedding.copy()
        labels, counts = np.unique(cluster_batches, return_counts=True)
        ref_batch = labels[np.argmax(counts)]
        ref_mask = cluster_batches == ref_batch
        ref_points = cluster_embedding[ref_mask]

        if ref_points.shape[0] == 0:
            return corrected_cluster

        for batch in labels:
            if batch == ref_batch:
                continue
            batch_mask = cluster_batches == batch
            batch_points = cluster_embedding[batch_mask]
            if batch_points.shape[0] == 0:
                continue

            try:
                cost = ot.dist(batch_points, ref_points, metric="sqeuclidean")
                source_weights = np.full(batch_points.shape[0], 1.0 / batch_points.shape[0])
                target_weights = np.full(ref_points.shape[0], 1.0 / ref_points.shape[0])
                transport = ot.sinkhorn(
                    source_weights,
                    target_weights,
                    cost,
                    reg=self.epsilon,
                )
                row_sums = np.maximum(transport.sum(axis=1, keepdims=True), 1e-12)
                aligned = (transport / row_sums) @ ref_points
                blended = self.alpha * aligned + (1.0 - self.alpha) * batch_points
                corrected_cluster[batch_mask] = self._replace_invalid(
                    blended,
                    fallback=batch_points,
                )
            except Exception:
                continue

        return corrected_cluster

    def _cluster_divergence(
        self,
        cluster_embedding: np.ndarray,
        cluster_batches: np.ndarray,
    ) -> float:
        labels = np.unique(cluster_batches)
        if labels.size < 2:
            return 0.0

        divergences: list[float] = []
        for left_index, left_batch in enumerate(labels):
            left_points = cluster_embedding[cluster_batches == left_batch]
            if left_points.shape[0] < 5:
                continue
            for right_batch in labels[left_index + 1 :]:
                right_points = cluster_embedding[cluster_batches == right_batch]
                if right_points.shape[0] < 5:
                    continue
                divergence = self._sinkhorn_divergence(left_points, right_points)
                divergences.append(max(1e-12, divergence))

        if not divergences:
            return 0.0

        raw_divergence = float(np.mean(divergences))
        if raw_divergence < 1e-10:
            return max(0.01, self.tau)
        return raw_divergence

    def _detect_global_batch_effect(
        self,
        embedding: np.ndarray,
        batches: np.ndarray,
    ) -> float:
        labels = np.unique(batches)
        if labels.size < 2:
            return 0.0

        divergences: list[float] = []
        for left_index, left_batch in enumerate(labels):
            left_points = embedding[batches == left_batch]
            if left_points.shape[0] < 5:
                continue
            for right_batch in labels[left_index + 1 :]:
                right_points = embedding[batches == right_batch]
                if right_points.shape[0] < 5:
                    continue
                divergence = self._sinkhorn_divergence(left_points, right_points)
                divergences.append(max(1e-12, divergence))

        if not divergences:
            return 0.0
        return float(np.mean(divergences))

    def _sinkhorn_divergence(
        self,
        left_points: np.ndarray,
        right_points: np.ndarray,
    ) -> float:
        try:
            cross_cost = ot.dist(left_points, right_points, metric="sqeuclidean")
            left_weights = np.full(left_points.shape[0], 1.0 / left_points.shape[0])
            right_weights = np.full(right_points.shape[0], 1.0 / right_points.shape[0])
            cross_value = ot.sinkhorn2(
                left_weights,
                right_weights,
                cross_cost,
                self.epsilon,
            )

            left_cost = ot.dist(left_points, left_points, metric="sqeuclidean")
            right_cost = ot.dist(right_points, right_points, metric="sqeuclidean")
            left_self = ot.sinkhorn2(left_weights, left_weights, left_cost, self.epsilon)
            right_self = ot.sinkhorn2(
                right_weights,
                right_weights,
                right_cost,
                self.epsilon,
            )
            divergence = float(cross_value - 0.5 * left_self - 0.5 * right_self)
            return max(0.0, divergence)
        except Exception:
            return 0.0

    @staticmethod
    def _replace_invalid(values: np.ndarray, fallback: np.ndarray) -> np.ndarray:
        invalid = np.isnan(values) | np.isinf(values)
        if not invalid.any():
            return values
        repaired = values.copy()
        repaired[invalid] = fallback[invalid]
        return repaired

    @staticmethod
    def _sanitize(values: np.ndarray) -> np.ndarray:
        return np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)

    @staticmethod
    def _validate_embedding(embedding: np.ndarray) -> np.ndarray:
        x = np.asarray(embedding, dtype=np.float64)
        if x.ndim != 2:
            raise ValueError("embedding must be a 2D array.")
        return SODA._sanitize(x)

    @staticmethod
    def _validate_batches(batches: np.ndarray, n_obs: int) -> np.ndarray:
        labels = np.asarray(batches).astype(str)
        if labels.ndim != 1:
            raise ValueError("batches must be a 1D array.")
        if labels.shape[0] != n_obs:
            raise ValueError("embedding and batches must contain the same number of cells.")
        return labels

    @staticmethod
    def _require_pot() -> None:
        if ot is None:
            raise ImportError(
                "POT is required for SODA. Install it with `pip install POT`."
            )

