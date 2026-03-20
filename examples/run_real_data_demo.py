"""Run a practical SODA demo on an h5ad dataset."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc

from soda import SODA, integrate_adata


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run SODA on an AnnData dataset and export example outputs.",
    )
    parser.add_argument("--input", required=True, help="Input .h5ad file.")
    parser.add_argument(
        "--output-dir",
        default="outputs/real_data_demo",
        help="Directory for corrected data and figures.",
    )
    parser.add_argument(
        "--batch-key",
        default="batch",
        help="Batch label column in adata.obs.",
    )
    parser.add_argument(
        "--celltype-key",
        default=None,
        help="Optional biological label column for additional UMAP plots.",
    )
    parser.add_argument(
        "--use-rep",
        default="X_pca",
        help="Embedding key in adata.obsm used as SODA input.",
    )
    parser.add_argument("--n-clusters", type=int, default=50)
    parser.add_argument("--max-iter", type=int, default=10)
    parser.add_argument("--epsilon", type=float, default=0.05)
    parser.add_argument("--tau", type=float, default=0.05)
    parser.add_argument("--alpha", type=float, default=0.7)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--n-pcs", type=int, default=50)
    return parser


def ensure_original_umap(adata, use_rep: str, n_pcs: int, random_state: int) -> None:
    if "X_umap" in adata.obsm:
        return
    if use_rep not in adata.obsm:
        sc.tl.pca(adata, n_comps=n_pcs, random_state=random_state)
        use_rep = "X_pca"
    sc.pp.neighbors(adata, use_rep=use_rep, random_state=random_state)
    sc.tl.umap(adata, random_state=random_state)


def save_umap(adata, color: str, output_path: Path) -> None:
    sc.pl.umap(adata, color=color, show=False)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def write_summary(
    output_path: Path,
    input_path: Path,
    adata_original,
    adata_soda,
    batch_key: str,
    celltype_key: str | None,
) -> None:
    lines = [
        "SODA demo summary",
        f"input_file: {input_path}",
        f"n_cells: {adata_original.n_obs}",
        f"n_genes: {adata_original.n_vars}",
        f"batch_key: {batch_key}",
        f"n_batches: {adata_original.obs[batch_key].astype(str).nunique()}",
        f"output_embedding: X_soda",
        f"n_iterations: {adata_soda.uns['soda']['X_soda']['n_iterations']}",
    ]
    if celltype_key and celltype_key in adata_original.obs.columns:
        lines.append(
            f"n_celltypes: {adata_original.obs[celltype_key].astype(str).nunique()}"
        )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = build_parser().parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(input_path)
    if args.batch_key not in adata.obs.columns:
        raise KeyError(f"Batch key '{args.batch_key}' was not found in adata.obs.")

    ensure_original_umap(
        adata=adata,
        use_rep=args.use_rep,
        n_pcs=args.n_pcs,
        random_state=args.random_state,
    )

    integrator = SODA(
        n_clusters=args.n_clusters,
        max_iter=args.max_iter,
        epsilon=args.epsilon,
        tau=args.tau,
        alpha=args.alpha,
        random_state=args.random_state,
    )

    adata_soda = integrate_adata(
        adata=adata,
        integrator=integrator,
        batch_key=args.batch_key,
        use_rep=args.use_rep,
        output_key="X_soda",
        copy=True,
        n_pcs=args.n_pcs,
        compute_neighbors=True,
        compute_leiden=True,
        compute_umap=True,
        umap_key="X_umap_soda",
    )

    adata_soda.write_h5ad(output_dir / "integrated.h5ad")

    save_umap(adata, args.batch_key, output_dir / "umap_original_batch.png")
    save_umap(adata_soda, args.batch_key, output_dir / "umap_soda_batch.png")

    if args.celltype_key and args.celltype_key in adata.obs.columns:
        save_umap(adata, args.celltype_key, output_dir / "umap_original_celltype.png")
        save_umap(
            adata_soda,
            args.celltype_key,
            output_dir / "umap_soda_celltype.png",
        )

    write_summary(
        output_path=output_dir / "summary.txt",
        input_path=input_path,
        adata_original=adata,
        adata_soda=adata_soda,
        batch_key=args.batch_key,
        celltype_key=args.celltype_key,
    )

    print(f"Demo complete. Outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()
