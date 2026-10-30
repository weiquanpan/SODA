#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Run reviewer Figure 1 pipeline on real pbmcsca subset with ALL methods
that are available in run_batch.py (i.e., benchmark default available methods).
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

from soda_benchmark.core.benchmark import SingleCellBatchBenchmark


def with_v1_prefix(name: str) -> str:
    return name if name.startswith("v1_") else f"v1_{name}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run reviewer Figure1 on real pbmcsca with all run_batch methods."
    )
    parser.add_argument(
        "--data",
        default="data/pbmcsca_real_10xv2_smartseq2.h5ad",
        help="Input h5ad path.",
    )
    parser.add_argument(
        "--output-root",
        default="results/rebuttal/data",
        help="Output root directory.",
    )
    parser.add_argument(
        "--batch-key",
        default="batch",
        help="Batch key in adata.obs (default: batch).",
    )
    parser.add_argument(
        "--celltype-key",
        default="CellType",
        help=(
            "Cell type key in adata.obs (default: CellType). "
            "If missing, script auto-falls back to common keys and then leiden."
        ),
    )
    parser.add_argument(
        "--skip-celltype-plot",
        action="store_true",
        help="Skip full-method celltype grid figure.",
    )
    parser.add_argument(
        "--batch-label-map",
        default="10x Chromium (v2):10x Chromium (v2),Smart-seq2:Smart-seq2",
        help="Batch label mapping string (key:label,key2:label2).",
    )
    parser.add_argument(
        "--neighborhood-mode",
        choices=["expression", "representation"],
        default="expression",
    )
    parser.add_argument("--n-bootstrap", type=int, default=100)
    parser.add_argument(
        "--rebuild-data",
        action="store_true",
        help="Force regenerate real pbmcsca h5ad.",
    )
    parser.add_argument(
        "--soda-tau",
        type=float,
        default=0.05,
        help="SODA tau (same as run_batch defaults).",
    )
    parser.add_argument(
        "--soda-alpha",
        type=float,
        default=0.7,
        help="SODA alpha (same as run_batch defaults).",
    )
    parser.add_argument(
        "--soda-epsilon",
        type=float,
        default=0.05,
        help="SODA epsilon (same as run_batch defaults).",
    )
    parser.add_argument(
        "--exclude-smai",
        action="store_true",
        help="Exclude SMAI method to avoid long/stuck runs on some environments.",
    )
    parser.add_argument(
        "--selected-methods",
        default="",
        help=(
            "Comma-separated method names to run, e.g. "
            "'SODA,Seurat_CCA,LIGER,Harmony,ComBat,BBKNN'. "
            "Empty means all available methods."
        ),
    )
    return parser.parse_args()


def parse_label_map(raw: str) -> dict[str, str]:
    if not raw:
        return {}
    mapping = {}
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError("Invalid batch-label-map format: key:label,key2:label2")
        key, val = item.split(":", 1)
        mapping[key.strip()] = val.strip()
    return mapping


def resolve_celltype_key(adata_dict: dict, preferred_key: str) -> str | None:
    """Pick a cell-type-like key available in all methods."""
    candidate_keys = [preferred_key, "CellType", "cell_type", "celltype", "labels", "leiden"]
    for key in candidate_keys:
        if not key:
            continue
        if all(key in adata.obs.columns for adata in adata_dict.values()):
            return key
    return None


def ensure_real_pbmcsca_h5ad(path: str, rebuild: bool = False) -> None:
    if os.path.exists(path) and not rebuild:
        return
    cmd = [
        "Rscript",
        "export_real_pbmcsca_subset_h5ad.R",
        f"--output={path}",
        "--nfeatures=2000",
        "--npcs=30",
        "--seed=42",
    ]
    print("Generating real pbmcsca subset h5ad...")
    subprocess.run(cmd, check=True)


def ensure_umap(adata) -> None:
    if "X_umap" in adata.obsm:
        return
    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata, n_comps=50, random_state=42)
    sc.pp.neighbors(adata, n_pcs=min(30, adata.obsm["X_pca"].shape[1]), random_state=42)
    sc.tl.umap(adata, random_state=42)


def draw_method_grid_batch(
    adata_dict: dict,
    out_png: str,
    out_pdf: str,
    batch_key: str,
    batch_label_map: dict[str, str],
) -> None:
    method_names = list(adata_dict.keys())
    n = len(method_names)
    ncols = 4 if n >= 4 else n
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.2 * nrows))
    axes = np.array(axes).reshape(-1)

    colors = {
        "10x Chromium (v2)": "#0072B2",
        "Smart-seq2": "#E69F00",
    }
    markers = ["o", "^", "s", "D", "P", "X", "v", "<", ">"]

    for idx, method in enumerate(method_names):
        ax = axes[idx]
        adata = adata_dict[method]
        ensure_umap(adata)
        if batch_key not in adata.obs.columns:
            raise KeyError(f"Batch key '{batch_key}' not found in method '{method}'.")
        batch_col = adata.obs[batch_key].astype(str)
        cats = sorted(batch_col.unique())
        for j, cat in enumerate(cats):
            mask = batch_col == cat
            label = batch_label_map.get(cat, cat)
            ax.scatter(
                adata.obsm["X_umap"][mask, 0],
                adata.obsm["X_umap"][mask, 1],
                s=6,
                alpha=0.75,
                c=colors.get(cat, f"C{j}"),
                marker=markers[j % len(markers)],
                edgecolors="none",
                label=label,
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
    plt.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def draw_method_grid_celltype(
    adata_dict: dict,
    out_png: str,
    out_pdf: str,
    celltype_key: str,
) -> None:
    method_names = list(adata_dict.keys())
    n = len(method_names)
    ncols = 4 if n >= 4 else n
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.2 * nrows))
    axes = np.array(axes).reshape(-1)

    for idx, method in enumerate(method_names):
        ax = axes[idx]
        adata = adata_dict[method]
        ensure_umap(adata)
        if celltype_key not in adata.obs.columns:
            raise KeyError(f"Celltype key '{celltype_key}' not found in method '{method}'.")
        ct_col = adata.obs[celltype_key].astype(str)
        cats = sorted(ct_col.unique())
        cmap = plt.get_cmap("tab20")
        texts = []
        for j, ct in enumerate(cats):
            mask = ct_col == ct
            xy = adata.obsm["X_umap"][mask, :2]
            ax.scatter(
                xy[:, 0],
                xy[:, 1],
                s=6,
                alpha=0.7,
                c=[cmap(j % 20)],
                edgecolors="none",
            )
            x0 = float(np.median(xy[:, 0]))
            y0 = float(np.median(xy[:, 1]))
            txt = ax.text(
                x0,
                y0,
                ct,
                fontsize=7,
                ha="center",
                va="center",
                weight="bold",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.55, pad=0.2),
            )
            texts.append(txt)

        try:
            from adjustText import adjust_text

            adjust_text(
                texts,
                ax=ax,
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.4, alpha=0.6),
                force_text=(0.15, 0.3),
            )
        except Exception:
            pass

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


def draw_four_panel_main(
    adata_dict: dict,
    out_png: str,
    out_pdf: str,
    batch_key: str,
    batch_label_map: dict[str, str],
) -> None:
    panel_order = [
        ("Original", "Uncorrected (Raw)"),
        ("Seurat_CCA", "Seurat v5 (CCA)"),
        ("LIGER", "LIGER (iNMF)"),
        ("SODA", "SODA (Ours)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    axes = axes.flatten()
    colors = {
        "10x Chromium (v2)": "#0072B2",
        "Smart-seq2": "#E69F00",
    }

    for ax, (method, title) in zip(axes, panel_order):
        if method not in adata_dict:
            ax.axis("off")
            ax.set_title(f"{title}\n(missing)")
            continue
        adata = adata_dict[method]
        ensure_umap(adata)
        batch_col = adata.obs[batch_key].astype(str)
        cats = sorted(batch_col.unique())
        for i, cat in enumerate(cats):
            mask = batch_col == cat
            ax.scatter(
                adata.obsm["X_umap"][mask, 0],
                adata.obsm["X_umap"][mask, 1],
                s=6,
                alpha=0.75,
                c=colors.get(cat, f"C{i}"),
                edgecolors="none",
                label=batch_label_map.get(cat, cat),
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


def write_batch_summary(adata, out_csv: str, batch_key: str) -> None:
    vc = adata.obs[batch_key].astype(str).value_counts()
    with open(out_csv, "w", encoding="utf-8") as f:
        f.write("dataset,batch_key,n_batches,batch_label,n_cells,batch_definition_note\n")
        for label, n_cells in vc.items():
            f.write(
                f"pbmcsca_real_10xv2_smartseq2,{batch_key},{len(vc)},"
                f"\"{label}\",{int(n_cells)},\"batch defined by sequencing technology (Method)\"\n"
            )


def main() -> None:
    args = parse_args()
    label_map = parse_label_map(args.batch_label_map)
    ensure_real_pbmcsca_h5ad(args.data, rebuild=args.rebuild_data)

    dataset_name = os.path.splitext(os.path.basename(args.data))[0]
    run_folder = with_v1_prefix(
        f"{dataset_name}_reviewer_full_fig1_{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    )
    out_dir = os.path.join(args.output_root, run_folder)
    os.makedirs(out_dir, exist_ok=True)

    metric_params = {
        "n_bootstrap": args.n_bootstrap,
        "ci_level": 95,
        "random_state": 42,
        "neighborhood_mode": args.neighborhood_mode,
    }
    method_params = {
        "SODA": {
            "tau": args.soda_tau,
            "alpha": args.soda_alpha,
            "epsilon": args.soda_epsilon,
        },
    }

    selected_methods = [m.strip() for m in args.selected_methods.split(",") if m.strip()]
    selected_methods = selected_methods or None

    # selected_methods=None means run all available methods (same behavior as run_batch default).
    bench = SingleCellBatchBenchmark(
        h5ad_path=args.data,
        method_params=method_params,
        metric_params=metric_params,
        selected_methods=selected_methods,
        batch_key=args.batch_key,
        batch_label_map=label_map,
        recompute_original_pca=True,
    )
    if args.exclude_smai and "SMAI" in bench.available_methods:
        del bench.available_methods["SMAI"]
        print("Excluding method: SMAI")
    bench.run_full_benchmark(output_dir=out_dir, generate_plots=False)

    # Main manuscript style 2x2 panel.
    draw_four_panel_main(
        adata_dict=bench.adata_dict,
        out_png=os.path.join(out_dir, "figure1_main_2x2_batch.png"),
        out_pdf=os.path.join(out_dir, "figure1_main_2x2_batch.pdf"),
        batch_key="batch",
        batch_label_map=label_map,
    )

    # Full-method reviewer panels.
    draw_method_grid_batch(
        adata_dict=bench.adata_dict,
        out_png=os.path.join(out_dir, "figure1_full_methods_grid_batch.png"),
        out_pdf=os.path.join(out_dir, "figure1_full_methods_grid_batch.pdf"),
        batch_key="batch",
        batch_label_map=label_map,
    )
    celltype_key = resolve_celltype_key(bench.adata_dict, args.celltype_key)
    if args.skip_celltype_plot:
        print("Skipping celltype plot because --skip-celltype-plot is set.")
    elif celltype_key is None:
        print(
            "Warning: No shared celltype-like key found across methods; "
            "skipping celltype grid plot."
        )
    else:
        if celltype_key != args.celltype_key:
            print(
                f"Requested celltype key '{args.celltype_key}' not available across all "
                f"methods, using fallback key '{celltype_key}'."
            )
        draw_method_grid_celltype(
            adata_dict=bench.adata_dict,
            out_png=os.path.join(out_dir, "figure1_full_methods_grid_celltype.png"),
            out_pdf=os.path.join(out_dir, "figure1_full_methods_grid_celltype.pdf"),
            celltype_key=celltype_key,
        )
    write_batch_summary(
        adata=bench.adata_dict["Original"],
        out_csv=os.path.join(out_dir, "dataset_batch_summary.csv"),
        batch_key="batch",
    )

    with open(os.path.join(out_dir, "methods_ran.txt"), "w", encoding="utf-8") as f:
        for name in bench.adata_dict.keys():
            f.write(f"{name}\n")

    print("Done.")
    print(f"Output directory: {out_dir}")
    print(f"Methods ran: {', '.join(bench.adata_dict.keys())}")


if __name__ == "__main__":
    main()
