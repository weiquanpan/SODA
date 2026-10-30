# -*- coding: utf-8 -*-
"""
主入口脚本：对单个数据集执行完整基准测试流程。

用法
----
    python run_benchmark.py
"""

import os
import sys
import argparse
from datetime import datetime

from soda_benchmark import SingleCellBatchBenchmark


def with_v1_prefix(name):
    """确保输出目录名统一带 v1_ 前缀。"""
    return name if name.startswith("v1_") else f"v1_{name}"


def parse_label_map(raw):
    """解析 --batch-label-map，格式: key1:label1,key2:label2"""
    if not raw:
        return {}
    mapping = {}
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(
                "Invalid --batch-label-map format. Use key:label,key2:label2"
            )
        k, v = item.split(":", 1)
        mapping[k.strip()] = v.strip()
    return mapping


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Run single-dataset batch correction benchmark."
    )
    parser.add_argument(
        "--data",
        default="data/pbmc4k.h5ad",
        help="Input h5ad file path (default: data/pbmc4k.h5ad)",
    )
    parser.add_argument(
        "--neighborhood-mode",
        choices=["expression", "representation"],
        default="expression",
        help=(
            "Neighborhood metric mode: "
            "'expression' compares adata.X; "
            "'representation' compares method embeddings in obsm."
        ),
    )
    parser.add_argument(
        "--soda-tau",
        type=float,
        default=0.05,
        help="SODA tau (divergence threshold) parameter (default: 0.05)",
    )
    parser.add_argument(
        "--soda-alpha",
        type=float,
        default=0.7,
        help="SODA soft update alpha parameter (default: 0.7)",
    )
    parser.add_argument(
        "--soda-epsilon",
        type=float,
        default=0.05,
        help="SODA entropy regularization epsilon parameter (default: 0.05)",
    )
    parser.add_argument(
        "--output-dir",
        default="results",
        help="Output root directory (default: results)",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Disable plot generation (plots are enabled by default).",
    )
    parser.add_argument(
        "--batch-key",
        default="batch",
        help="Batch column name in adata.obs (default: batch).",
    )
    parser.add_argument(
        "--batch-label-map",
        default="",
        help="Optional legend label mapping: key1:label1,key2:label2",
    )
    parser.add_argument(
        "--recompute-original-pca",
        action="store_true",
        help="Recompute uncorrected PCA/UMAP for original baseline.",
    )
    parser.add_argument(
        "--label-key",
        default=None,
        help=(
            "Optional label column for cLISI/ASW_batch "
            "(e.g., cell_type or leiden). Defaults to auto-detect."
        ),
    )
    parser.add_argument(
        "--lisi-n-neighbors",
        type=int,
        default=30,
        help="Number of neighbors for iLISI/cLISI (default: 30).",
    )
    args = parser.parse_args()

    h5ad_path = args.data
    # h5ad_path = "data/simulated_data.h5ad"

    if not os.path.exists(h5ad_path):
        print(f"Error: Data file not found: {h5ad_path}")
        print("\nPlease ensure you have:")
        print("1. Downloaded the required single-cell dataset")
        print("2. Run preprocessing to generate the H5AD file")
        sys.exit(1)

    metric_params = {
        "n_bootstrap": 1000,
        "ci_level": 95,
        "random_state": 42,
        "neighborhood_mode": args.neighborhood_mode,
        "label_key": args.label_key,
        "lisi_n_neighbors": args.lisi_n_neighbors,
    }
    method_params = {
        "SODA": {
            "tau": args.soda_tau,
            "alpha": args.soda_alpha,
            "epsilon": args.soda_epsilon,
        }
    }

    benchmark = SingleCellBatchBenchmark(
        h5ad_path=h5ad_path,
        method_params=method_params,
        metric_params=metric_params,
        batch_key=args.batch_key,
        batch_label_map=parse_label_map(args.batch_label_map),
        recompute_original_pca=args.recompute_original_pca,
    )
    dataset_name = os.path.splitext(os.path.basename(h5ad_path))[0]
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    run_folder = with_v1_prefix(f"{dataset_name}_{timestamp}")
    output_dir = os.path.join(args.output_dir, run_folder)

    benchmark.run_full_benchmark(
        output_dir=output_dir,
        generate_plots=not args.no_plots,
    )


if __name__ == '__main__':
    main()

