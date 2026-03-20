"""Command-line entry point for running SODA on h5ad files."""

from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc

from soda.model import SODA
from soda.scanpy import integrate_adata


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Apply SODA batch integration to an h5ad dataset.",
    )
    parser.add_argument("--input", required=True, help="Input .h5ad file.")
    parser.add_argument(
        "--output",
        default=None,
        help="Output .h5ad file. Defaults to <input>_soda.h5ad.",
    )
    parser.add_argument(
        "--batch-key",
        default="batch",
        help="Column in adata.obs that stores batch labels.",
    )
    parser.add_argument(
        "--use-rep",
        default="X_pca",
        help="Embedding in adata.obsm used as SODA input.",
    )
    parser.add_argument(
        "--output-key",
        default="X_soda",
        help="Key used to store the corrected embedding in adata.obsm.",
    )
    parser.add_argument("--n-clusters", type=int, default=50)
    parser.add_argument("--max-iter", type=int, default=10)
    parser.add_argument("--epsilon", type=float, default=0.05)
    parser.add_argument("--tau", type=float, default=0.05)
    parser.add_argument("--alpha", type=float, default=0.7)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=50,
        help="Number of PCs to compute if the requested representation is missing.",
    )
    parser.add_argument(
        "--no-neighbors",
        action="store_true",
        help="Skip neighbor graph construction on the corrected embedding.",
    )
    parser.add_argument(
        "--no-leiden",
        action="store_true",
        help="Skip Leiden clustering.",
    )
    parser.add_argument(
        "--no-umap",
        action="store_true",
        help="Skip UMAP generation.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else input_path.with_name(
        f"{input_path.stem}_soda.h5ad"
    )

    integrator = SODA(
        n_clusters=args.n_clusters,
        max_iter=args.max_iter,
        epsilon=args.epsilon,
        tau=args.tau,
        alpha=args.alpha,
        random_state=args.random_state,
    )

    adata = sc.read_h5ad(input_path)
    corrected = integrate_adata(
        adata=adata,
        integrator=integrator,
        batch_key=args.batch_key,
        use_rep=args.use_rep,
        output_key=args.output_key,
        copy=True,
        n_pcs=args.n_pcs,
        compute_neighbors=not args.no_neighbors,
        compute_leiden=not args.no_leiden,
        compute_umap=not args.no_umap,
    )
    corrected.write_h5ad(output_path)

    print(f"SODA integration finished. Output saved to: {output_path}")

