"""
Generate benchmark-ready simulated scRNA-seq data for this project.


"""

import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def _sample_cell_types(rng, n_cell_types, n_cells, composition):
    probs = np.asarray(composition, dtype=np.float64)
    probs = probs / probs.sum()
    return rng.choice(n_cell_types, size=n_cells, replace=True, p=probs)


def _build_base_means(
    rng,
    n_cell_types,
    n_genes,
    n_marker_genes_per_ct,
    marker_mean_range,
    background_mean_range,
):
    max_marker_per_ct = max(1, n_genes // n_cell_types)
    marker_per_ct = min(n_marker_genes_per_ct, max_marker_per_ct)
    n_marker_total = n_cell_types * marker_per_ct

    marker_assignments = np.full(n_genes, -1, dtype=int)
    marker_assignments[:n_marker_total] = np.repeat(
        np.arange(n_cell_types), marker_per_ct
    )

    mu_base = np.zeros((n_cell_types, n_genes), dtype=np.float64)
    for ct in range(n_cell_types):
        mu_base[ct, :] = rng.uniform(
            background_mean_range[0], background_mean_range[1], size=n_genes
        )
        marker_idx = np.where(marker_assignments == ct)[0]
        if len(marker_idx) > 0:
            mu_base[ct, marker_idx] = rng.uniform(
                marker_mean_range[0], marker_mean_range[1], size=len(marker_idx)
            )
    return mu_base, marker_assignments


def simulate_cellanova_data(
    n_cell_types: int = 6,
    n_control_batches: int = 5,
    n_treatment_batches: int = 2,
    n_genes: int = 5000,
    n_marker_genes_per_ct: int = 60,
    cells_per_batch: int = 3000,
    theta: float = 1.0,
    marker_mean_range: tuple = (2.0, 4.0),
    background_mean_range: tuple = (0.2, 1.2),
    treatment_effect_additive: float = 1.25,
    batch_effect_log_mean: float = 0.0,
    batch_effect_log_sd: float = 0.35,
    seed: int = 42,
):
    """
    Simulate counts with biological + batch + treatment structure.

    Returns raw-count AnnData (before preprocessing).
    """
    rng = np.random.default_rng(seed)
    n_batches = n_control_batches + n_treatment_batches
    if n_batches < 2:
        raise ValueError("Need at least 2 batches for integration benchmark.")

    mu_base, marker_assignments = _build_base_means(
        rng=rng,
        n_cell_types=n_cell_types,
        n_genes=n_genes,
        n_marker_genes_per_ct=n_marker_genes_per_ct,
        marker_mean_range=marker_mean_range,
        background_mean_range=background_mean_range,
    )

    # CT6 marker genes carry treatment signal in treatment batches.
    target_ct = n_cell_types - 1
    target_marker_idx = np.where(marker_assignments == target_ct)[0]
    mu_treatment_target = mu_base[target_ct].copy()
    mu_treatment_target[target_marker_idx] += treatment_effect_additive

    all_counts = []
    all_obs = []

    # Slight composition shift between control/treatment to make task realistic.
    base_comp = np.ones(n_cell_types, dtype=np.float64)
    comp_control = base_comp.copy()
    comp_treatment = base_comp.copy()
    comp_treatment[target_ct] *= 1.35

    for batch_id in range(n_batches):
        is_treatment = batch_id >= n_control_batches
        condition = "treatment" if is_treatment else "control"
        batch_label = batch_id  # int labels are compatible with benchmark prints.

        comp = comp_treatment if is_treatment else comp_control
        ct_labels = _sample_cell_types(
            rng=rng,
            n_cell_types=n_cell_types,
            n_cells=cells_per_batch,
            composition=comp,
        )

        # Additive random shift on log-scale -> multiplicative effect on mean.
        log_shift = rng.normal(
            loc=batch_effect_log_mean, scale=batch_effect_log_sd, size=n_genes
        )
        batch_multiplier = np.exp(log_shift)

        counts = np.zeros((cells_per_batch, n_genes), dtype=np.int32)
        for i, ct in enumerate(ct_labels):
            if is_treatment and ct == target_ct:
                mu = mu_treatment_target.copy()
            else:
                mu = mu_base[ct].copy()

            # Cell-specific library-size factor.
            lib_factor = np.exp(rng.normal(loc=0.0, scale=0.3))
            mu = np.maximum(mu * batch_multiplier * lib_factor, 1e-8)

            p = theta / (theta + mu)
            counts[i, :] = rng.negative_binomial(theta, p, size=n_genes)

        all_counts.append(counts)
        all_obs.append(
            pd.DataFrame(
                {
                    "batch": batch_label,
                    "condition": condition,
                    "cell_type": [f"CT{x + 1}" for x in ct_labels],
                }
            )
        )

    X_counts = np.vstack(all_counts)
    obs = pd.concat(all_obs, axis=0, ignore_index=True)
    obs.index = [f"cell_{i}" for i in range(obs.shape[0])]

    is_marker = marker_assignments >= 0
    marker_for = np.array(["background"] * n_genes, dtype=object)
    marker_for[is_marker] = [f"CT{x + 1}" for x in marker_assignments[is_marker]]

    var = pd.DataFrame(
        {
            "is_marker": is_marker,
            "marker_for": marker_for,
            "is_target_ct_deg": marker_for == f"CT{target_ct + 1}",
        },
        index=[f"gene_{i}" for i in range(n_genes)],
    )

    adata = ad.AnnData(X=sparse.csr_matrix(X_counts), obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    return adata


def preprocess_for_benchmark(adata, n_top_hvg=2000):
    """
    Convert raw-count AnnData into benchmark-ready preprocessed AnnData.
    """
    adata = adata.copy()
    adata.raw = adata.copy()  # keep all original genes/counts

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    top_hvg = min(n_top_hvg, adata.n_vars)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=top_hvg,
        flavor="cell_ranger",
        layer="counts",
    )
    if "highly_variable" not in adata.var:
        adata.var["highly_variable"] = True
    if int(adata.var["highly_variable"].sum()) == 0:
        adata.var["highly_variable"] = True

    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10.0)
    sc.tl.pca(adata, n_comps=min(50, adata.n_vars), random_state=42)
    sc.pp.neighbors(adata, n_pcs=min(30, adata.obsm["X_pca"].shape[1]), random_state=42)
    sc.tl.umap(adata, random_state=42)
    sc.tl.leiden(adata, resolution=0.5, random_state=42)
    return adata


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate benchmark-ready simulated scRNA-seq data."
    )
    parser.add_argument(
        "--output",
        default="data/simulated_cellanova_benchmark.h5ad",
        help="Output h5ad path (default: data/simulated_cellanova_benchmark.h5ad).",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--cells-per-batch", type=int, default=3000)
    parser.add_argument("--genes", type=int, default=5000)
    parser.add_argument(
        "--control-batches",
        type=int,
        default=1,
        help="Control batch count (default 1).",
    )
    parser.add_argument(
        "--treatment-batches",
        type=int,
        default=1,
        help="Treatment batch count (default 1).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    adata_counts = simulate_cellanova_data(
        n_control_batches=args.control_batches,
        n_treatment_batches=args.treatment_batches,
        n_genes=args.genes,
        cells_per_batch=args.cells_per_batch,
        seed=args.seed,
    )
    adata_ready = preprocess_for_benchmark(adata_counts)

    output_path = args.output
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    adata_ready.write_h5ad(output_path)

    print("\nGenerated benchmark-ready dataset")
    print(f"  shape: {adata_ready.n_obs} cells x {adata_ready.n_vars} HVGs")
    print(f"  output: {output_path}")
    print(f"  batch counts: {adata_ready.obs['batch'].value_counts().to_dict()}")
    print(f"  condition counts: {adata_ready.obs['condition'].value_counts().to_dict()}")
    print(f"  cell types: {adata_ready.obs['cell_type'].value_counts().sort_index().to_dict()}")
    print(f"  clusters (leiden): {adata_ready.obs['leiden'].nunique()}")


if __name__ == "__main__":
    main()
