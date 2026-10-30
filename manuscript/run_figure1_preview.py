# -*- coding: utf-8 -*-
"""
Figure 1 预览脚本（先跑一个数据集验证 UMAP 展示效果）。

默认：
  - 数据：data/heart.h5ad
  - 输出：results/rebuttal/data
  - 强制重算 Original 的 uncorrected PCA/UMAP
  - 生成 Figure 1 相关图：umap_comparison.png / umap_clusters.png
"""

import argparse
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Run one-dataset benchmark preview for Figure 1."
    )
    parser.add_argument(
        "--data",
        default="data/heart.h5ad",
        help="Input h5ad file path (default: data/heart.h5ad)",
    )
    parser.add_argument(
        "--output-root",
        default="results/rebuttal/data",
        help="Output root directory (default: results/rebuttal/data)",
    )
    parser.add_argument(
        "--batch-key",
        default="batch",
        help="Batch column name in adata.obs (default: batch).",
    )
    parser.add_argument(
        "--batch-label-map",
        default="0:Batch 0,1:Batch 1,0.0:Batch 0,1.0:Batch 1",
        help="Legend label mapping: key1:label1,key2:label2",
    )
    parser.add_argument(
        "--neighborhood-mode",
        choices=["expression", "representation"],
        default="expression",
        help="Neighborhood metric mode.",
    )
    args = parser.parse_args()

    cmd = [
        sys.executable,
        "run_benchmark.py",
        "--data",
        args.data,
        "--output-dir",
        args.output_root,
        "--batch-key",
        args.batch_key,
        "--batch-label-map",
        args.batch_label_map,
        "--recompute-original-pca",
        "--neighborhood-mode",
        args.neighborhood_mode,
    ]

    print("Running command:")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
