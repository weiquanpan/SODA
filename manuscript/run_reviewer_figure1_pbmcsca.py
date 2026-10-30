#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generate reviewer-focused Figure 1 on real pbmcsca subset:
Uncorrected vs Seurat CCA vs LIGER vs SODA.
"""

from __future__ import annotations

import argparse
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
        description="Run reviewer Figure 1 on real pbmcsca subset."
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
        "--batch-label-map",
        default="10x Chromium (v2):10x Chromium (v2),Smart-seq2:Smart-seq2",
        help="Batch label mapping string (key:label,key2:label2).",
    )
    parser.add_argument(
        "--neighborhood-mode",
        choices=["expression", "representation"],
        default="expression",
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=300,
        help="Bootstrap rounds for metrics.",
    )
    parser.add_argument(
        "--rebuild-data",
        action="store_true",
        help="Force regenerate real pbmcsca h5ad before benchmark.",
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


def ensure_umap(adata, batch_key: str) -> None:
    if "X_umap" in adata.obsm:
        return
    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata, n_comps=50, random_state=42)
    if batch_key not in adata.obs.columns and "batch" in adata.obs.columns:
        batch_key = "batch"
    sc.pp.neighbors(adata, n_pcs=min(30, adata.obsm["X_pca"].shape[1]), random_state=42)
    sc.tl.umap(adata, random_state=42)


def draw_four_panel(adata_dict: dict, out_png: str, out_pdf: str, batch_key: str, batch_label_map: dict) -> None:
    panel_order = [
        ("Original", "Uncorrected (Raw)"),
        ("Seurat_CCA", "Seurat v5 (CCA)"),
        ("LIGER", "LIGER (iNMF)"),
        ("SODA", "SODA (Ours)"),
    ]
    colors = {
        "10x Chromium (v2)": "#0072B2",
        "Smart-seq2": "#E69F00",
    }

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    axes = axes.flatten()

    for ax, (method, title) in zip(axes, panel_order):
        if method not in adata_dict:
            ax.axis("off")
            ax.set_title(f"{title}\n(missing)")
            continue

        adata = adata_dict[method]
        ensure_umap(adata, batch_key=batch_key)
        if batch_key not in adata.obs.columns:
            raise KeyError(f"Batch key '{batch_key}' not found in method '{method}'.")

        batch_col = adata.obs[batch_key].astype(str)
        categories = sorted(batch_col.unique())
        for i, cat in enumerate(categories):
            mask = batch_col == cat
            label = batch_label_map.get(cat, cat)
            color = colors.get(cat, f"C{i}")
            ax.scatter(
                adata.obsm["X_umap"][mask, 0],
                adata.obsm["X_umap"][mask, 1],
                s=6,
                alpha=0.75,
                c=color,
                label=label,
                edgecolors="none",
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


def draw_four_panel_celltype(adata_dict: dict, out_png: str, out_pdf: str, celltype_key: str = "CellType") -> None:
    panel_order = [
        ("Original", "Uncorrected (Raw)"),
        ("Seurat_CCA", "Seurat v5 (CCA)"),
        ("LIGER", "LIGER (iNMF)"),
        ("SODA", "SODA (Ours)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    axes = axes.flatten()

    for ax, (method, title) in zip(axes, panel_order):
        if method not in adata_dict:
            ax.axis("off")
            ax.set_title(f"{title}\n(missing)")
            continue

        adata = adata_dict[method]
        ensure_umap(adata, batch_key="batch")
        if celltype_key not in adata.obs.columns:
            raise KeyError(f"Cell type key '{celltype_key}' not found in method '{method}'.")

        ct_col = adata.obs[celltype_key].astype(str)
        categories = sorted(ct_col.unique())
        palette = plt.get_cmap("tab20")
        label_positions = []
        texts = []
        for i, ct in enumerate(categories):
            mask = ct_col == ct
            xy = adata.obsm["X_umap"][mask, :2]
            color = palette(i % 20)
            ax.scatter(
                xy[:, 0],
                xy[:, 1],
                s=6,
                alpha=0.7,
                c=[color],
                edgecolors="none",
            )
            x_median = float(np.median(xy[:, 0]))
            y_median = float(np.median(xy[:, 1]))
            label_positions.append((x_median, y_median))
            txt = ax.text(
                x_median,
                y_median,
                ct,
                fontsize=8,
                weight="bold",
                ha="center",
                va="center",
                color="black",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.55, pad=0.3),
            )
            texts.append(txt)

        # Optional repel if adjustText is installed.
        try:
            from adjustText import adjust_text

            adjust_text(
                texts,
                ax=ax,
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.6),
                only_move={"points": "y", "text": "xy"},
                force_text=(0.2, 0.4),
            )
        except Exception:
            # Fallback without dependency; keep median labels.
            pass

        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.grid(alpha=0.2)
        ax.set_aspect("equal", "box")

    plt.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_batch_summary(adata, out_csv: str, batch_key: str) -> None:
    batch_counts = adata.obs[batch_key].astype(str).value_counts()
    with open(out_csv, "w", encoding="utf-8") as f:
        f.write("dataset,batch_key,n_batches,batch_label,n_cells,batch_definition_note\n")
        for label, n_cells in batch_counts.items():
            f.write(
                f"pbmcsca_real_10xv2_smartseq2,{batch_key},{len(batch_counts)},"
                f"\"{label}\",{int(n_cells)},\"batch defined by sequencing technology (Method)\"\n"
            )


def main() -> None:
    args = parse_args()
    label_map = parse_label_map(args.batch_label_map)
    ensure_real_pbmcsca_h5ad(args.data, rebuild=args.rebuild_data)

    dataset_name = os.path.splitext(os.path.basename(args.data))[0]
    run_folder = with_v1_prefix(f"{dataset_name}_reviewer_fig1_{datetime.now().strftime('%Y%m%d-%H%M%S')}")
    out_dir = os.path.join(args.output_root, run_folder)
    os.makedirs(out_dir, exist_ok=True)

    metric_params = {
        "n_bootstrap": args.n_bootstrap,
        "ci_level": 95,
        "random_state": 42,
        "neighborhood_mode": args.neighborhood_mode,
    }
    method_params = {
        "SODA": {"tau": 0.05, "alpha": 0.7, "epsilon": 0.05},
    }
    selected_methods = ["Seurat_CCA", "LIGER", "SODA"]

    bench = SingleCellBatchBenchmark(
        h5ad_path=args.data,
        method_params=method_params,
        metric_params=metric_params,
        selected_methods=selected_methods,
        batch_key=args.batch_key,
        batch_label_map=label_map,
        recompute_original_pca=True,
    )
    bench.run_full_benchmark(output_dir=out_dir, generate_plots=False)

    out_png = os.path.join(out_dir, "figure1_reviewer_real_pbmcsca_2x2.png")
    out_pdf = os.path.join(out_dir, "figure1_reviewer_real_pbmcsca_2x2.pdf")
    draw_four_panel(
        adata_dict=bench.adata_dict,
        out_png=out_png,
        out_pdf=out_pdf,
        batch_key="batch",
        batch_label_map=label_map,
    )
    out_bio_png = os.path.join(out_dir, "figure1_reviewer_real_pbmcsca_2x2_celltype.png")
    out_bio_pdf = os.path.join(out_dir, "figure1_reviewer_real_pbmcsca_2x2_celltype.pdf")
    draw_four_panel_celltype(
        adata_dict=bench.adata_dict,
        out_png=out_bio_png,
        out_pdf=out_bio_pdf,
        celltype_key="CellType",
    )
    write_batch_summary(
        adata=bench.adata_dict["Original"],
        out_csv=os.path.join(out_dir, "dataset_batch_summary.csv"),
        batch_key="batch",
    )

    print("Done.")
    print(f"Output directory: {out_dir}")
    print(f"Figure (PNG): {out_png}")
    print(f"Figure (PDF): {out_pdf}")
    print(f"Figure CellType (PNG): {out_bio_png}")
    print(f"Figure CellType (PDF): {out_bio_pdf}")


if __name__ == "__main__":
    main()
