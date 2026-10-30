#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import math
import os
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scanpy.external as sce


def with_v1_prefix(name: str) -> str:
    return name if name.startswith("v1_") else f"v1_{name}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Quick reviewer-style Figure1 plotting for pbmc4k."
    )
    parser.add_argument("--data", default="data/pbmc4k.h5ad", help="Input h5ad path.")
    parser.add_argument("--output-root", default="results/rebuttal/data", help="Output root.")
    parser.add_argument("--batch-key", default="batch", help="Batch key in adata.obs.")
    parser.add_argument("--celltype-key", default="leiden", help="Color key for type-like panel.")
    return parser.parse_args()


def ensure_umap_from_pca(adata) -> None:
    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata, n_comps=50, random_state=42)
    sc.pp.neighbors(adata, n_pcs=min(30, adata.obsm["X_pca"].shape[1]), random_state=42)
    sc.tl.umap(adata, random_state=42)


def ensure_leiden(adata) -> None:
    if "leiden" in adata.obs.columns:
        return
    if "X_umap" not in adata.obsm:
        ensure_umap_from_pca(adata)
    sc.tl.leiden(adata, random_state=42)


def draw_method_grid_batch(adata_dict: dict, out_png: str, out_pdf: str, batch_key: str) -> None:
    methods = list(adata_dict.keys())
    n = len(methods)
    ncols = 4 if n >= 4 else n
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.2 * nrows))
    axes = np.array(axes).reshape(-1)
    for i, method in enumerate(methods):
        ax = axes[i]
        adata = adata_dict[method]
        cats = sorted(adata.obs[batch_key].astype(str).unique())
        for j, cat in enumerate(cats):
            mask = adata.obs[batch_key].astype(str) == cat
            ax.scatter(
                adata.obsm["X_umap"][mask, 0],
                adata.obsm["X_umap"][mask, 1],
                s=6,
                alpha=0.75,
                c=f"C{j}",
                edgecolors="none",
                label=cat,
            )
        ax.set_title(method, fontsize=11, fontweight="bold")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.grid(alpha=0.2)
        ax.set_aspect("equal", "box")
    for k in range(n, len(axes)):
        axes[k].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=min(4, len(labels)), frameon=False)
    plt.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def draw_method_grid_celltype(
    adata_dict: dict, out_png: str, out_pdf: str, celltype_key: str
) -> None:
    methods = list(adata_dict.keys())
    n = len(methods)
    ncols = 4 if n >= 4 else n
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.2 * nrows))
    axes = np.array(axes).reshape(-1)
    cmap = plt.get_cmap("tab20")
    for i, method in enumerate(methods):
        ax = axes[i]
        adata = adata_dict[method]
        ct_col = adata.obs[celltype_key].astype(str)
        cats = sorted(ct_col.unique())
        for j, ct in enumerate(cats):
            mask = ct_col == ct
            xy = adata.obsm["X_umap"][mask, :2]
            ax.scatter(xy[:, 0], xy[:, 1], s=6, alpha=0.7, c=[cmap(j % 20)], edgecolors="none")
        ax.set_title(method, fontsize=11, fontweight="bold")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.grid(alpha=0.2)
        ax.set_aspect("equal", "box")
    for k in range(n, len(axes)):
        axes[k].axis("off")
    plt.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def draw_four_panel_main(adata_dict: dict, out_png: str, out_pdf: str, batch_key: str) -> None:
    panel_order = [
        ("Original", "Uncorrected (Raw)"),
        ("ComBat", "ComBat"),
        ("Harmony", "Harmony"),
        ("BBKNN", "BBKNN"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    axes = axes.flatten()
    for ax, (method, title) in zip(axes, panel_order):
        if method not in adata_dict:
            ax.axis("off")
            ax.set_title(f"{title}\n(missing)")
            continue
        adata = adata_dict[method]
        cats = sorted(adata.obs[batch_key].astype(str).unique())
        for i, cat in enumerate(cats):
            mask = adata.obs[batch_key].astype(str) == cat
            ax.scatter(
                adata.obsm["X_umap"][mask, 0],
                adata.obsm["X_umap"][mask, 1],
                s=6,
                alpha=0.75,
                c=f"C{i}",
                edgecolors="none",
                label=cat,
            )
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.grid(alpha=0.2)
        ax.set_aspect("equal", "box")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    plt.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    adata = sc.read_h5ad(args.data)
    if args.batch_key not in adata.obs.columns:
        raise KeyError(f"batch key '{args.batch_key}' not found.")
    adata.obs["batch"] = adata.obs[args.batch_key].astype(str)

    data_name = os.path.splitext(os.path.basename(args.data))[0]
    run_folder = with_v1_prefix(f"{data_name}_reviewer_quick_fig1_{datetime.now().strftime('%Y%m%d-%H%M%S')}")
    out_dir = os.path.join(args.output_root, run_folder)
    os.makedirs(out_dir, exist_ok=True)

    adata_dict = {}

    ad0 = adata.copy()
    ensure_umap_from_pca(ad0)
    ensure_leiden(ad0)
    adata_dict["Original"] = ad0

    ad1 = adata.copy()
    sc.pp.combat(ad1, key="batch")
    ensure_umap_from_pca(ad1)
    ensure_leiden(ad1)
    adata_dict["ComBat"] = ad1

    try:
        ad2 = adata.copy()
        if "X_pca" not in ad2.obsm:
            sc.tl.pca(ad2, n_comps=50, random_state=42)
        sce.pp.harmony_integrate(ad2, "batch")
        sc.pp.neighbors(ad2, use_rep="X_pca_harmony", random_state=42)
        sc.tl.umap(ad2, random_state=42)
        ensure_leiden(ad2)
        adata_dict["Harmony"] = ad2
    except Exception as e:
        print(f"Harmony skipped due to error: {e}")

    try:
        ad3 = adata.copy()
        if "X_pca" not in ad3.obsm:
            sc.tl.pca(ad3, n_comps=50, random_state=42)
        sce.pp.bbknn(ad3, batch_key="batch")
        sc.tl.umap(ad3, random_state=42)
        ensure_leiden(ad3)
        adata_dict["BBKNN"] = ad3
    except Exception as e:
        print(f"BBKNN skipped due to error: {e}")

    if args.celltype_key not in adata_dict["Original"].obs.columns:
        args.celltype_key = "leiden"

    draw_four_panel_main(
        adata_dict=adata_dict,
        out_png=os.path.join(out_dir, "figure1_main_2x2_batch.png"),
        out_pdf=os.path.join(out_dir, "figure1_main_2x2_batch.pdf"),
        batch_key="batch",
    )
    draw_method_grid_batch(
        adata_dict=adata_dict,
        out_png=os.path.join(out_dir, "figure1_full_methods_grid_batch.png"),
        out_pdf=os.path.join(out_dir, "figure1_full_methods_grid_batch.pdf"),
        batch_key="batch",
    )
    draw_method_grid_celltype(
        adata_dict=adata_dict,
        out_png=os.path.join(out_dir, "figure1_full_methods_grid_celltype.png"),
        out_pdf=os.path.join(out_dir, "figure1_full_methods_grid_celltype.pdf"),
        celltype_key=args.celltype_key,
    )

    with open(os.path.join(out_dir, "methods_ran.txt"), "w", encoding="utf-8") as f:
        for m in adata_dict.keys():
            f.write(f"{m}\n")

    print("Done.")
    print(f"Output directory: {out_dir}")
    print(f"Methods ran: {', '.join(adata_dict.keys())}")


if __name__ == "__main__":
    main()

