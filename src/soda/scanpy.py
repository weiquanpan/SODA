"""Scanpy-facing helpers for applying SODA to AnnData objects."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from soda.model import SODA

if TYPE_CHECKING:
    from anndata import AnnData


def integrate_adata(
    adata: "AnnData",
    integrator: SODA,
    batch_key: str = "batch",
    use_rep: str = "X_pca",
    output_key: str = "X_soda",
    copy: bool = True,
    n_pcs: int = 50,
    compute_neighbors: bool = True,
    neighbors_key: str = "soda_neighbors",
    compute_leiden: bool = True,
    leiden_key: str = "leiden_soda",
    leiden_resolution: float = 0.5,
    compute_umap: bool = True,
    umap_key: str = "X_umap_soda",
) -> "AnnData":
    """Apply SODA to an AnnData object and store results in `obsm`."""
    import scanpy as sc

    if batch_key not in adata.obs.columns:
        raise KeyError(f"Batch key '{batch_key}' was not found in adata.obs.")

    target = adata.copy() if copy else adata

    if use_rep not in target.obsm:
        sc.tl.pca(target, n_comps=n_pcs, random_state=integrator.random_state)
        use_rep = "X_pca"

    embedding = np.asarray(target.obsm[use_rep], dtype=np.float64)
    batches = target.obs[batch_key].astype(str).to_numpy()
    result = integrator.integrate(embedding, batches)
    target.obsm[output_key] = result.embedding

    target.uns.setdefault("soda", {})
    target.uns["soda"][output_key] = {
        "batch_key": batch_key,
        "input_representation": use_rep,
        "n_clusters": integrator.n_clusters,
        "max_iter": integrator.max_iter,
        "epsilon": integrator.epsilon,
        "tau": integrator.tau,
        "alpha": integrator.alpha,
        "random_state": integrator.random_state,
        "n_iterations": result.n_iterations,
        "correction_history": result.correction_history,
    }

    if compute_neighbors:
        sc.pp.neighbors(
            target,
            use_rep=output_key,
            key_added=neighbors_key,
            random_state=integrator.random_state,
        )

    if compute_leiden:
        if not compute_neighbors:
            raise ValueError("compute_neighbors must be True when compute_leiden is True.")
        sc.tl.leiden(
            target,
            resolution=leiden_resolution,
            neighbors_key=neighbors_key,
            key_added=leiden_key,
            random_state=integrator.random_state,
        )

    if compute_umap:
        if not compute_neighbors:
            raise ValueError("compute_neighbors must be True when compute_umap is True.")
        sc.tl.umap(
            target,
            neighbors_key=neighbors_key,
            random_state=integrator.random_state,
        )
        if umap_key != "X_umap":
            target.obsm[umap_key] = target.obsm["X_umap"].copy()

    return target

