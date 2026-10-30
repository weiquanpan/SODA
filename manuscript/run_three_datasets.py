# -*- coding: utf-8 -*-
"""
一键运行固定三数据集：
  - data/pbmc4k.h5ad
  - data/heart.h5ad
  - data/neuro.h5ad

默认执行普通批量基准；可选开启 tau 敏感性扫描。
"""

import argparse
import subprocess
import sys


DATASETS = ["pbmc4k.h5ad", "heart.h5ad", "neuro.h5ad"]


def main():
    parser = argparse.ArgumentParser(
        description="Run benchmark on pbmc4k/heart/neuro datasets."
    )
    parser.add_argument(
        "--neighborhood-mode",
        choices=["expression", "representation"],
        default="expression",
        help="Neighborhood metric mode.",
    )
    parser.add_argument(
        "--soda-tau",
        type=float,
        default=0.05,
        help="Default SODA tau for normal benchmark mode.",
    )
    parser.add_argument(
        "--soda-alpha",
        type=float,
        default=0.7,
        help="Default SODA alpha for normal benchmark mode.",
    )
    parser.add_argument(
        "--soda-epsilon",
        type=float,
        default=0.05,
        help="Default SODA epsilon for normal benchmark mode.",
    )
    parser.add_argument(
        "--tau-sensitivity",
        action="store_true",
        help="Enable reviewer-style tau sweep: 0.01,0.02,0.05,0.1,0.2",
    )
    parser.add_argument(
        "--tau-values",
        default="0.01,0.02,0.05,0.1,0.2",
        help="Comma-separated tau values for sensitivity sweep.",
    )
    args = parser.parse_args()

    cmd = [
        sys.executable,
        "run_batch.py",
        "--neighborhood-mode",
        args.neighborhood_mode,
        "--soda-tau",
        str(args.soda_tau),
        "--soda-alpha",
        str(args.soda_alpha),
        "--soda-epsilon",
        str(args.soda_epsilon),
        "--data-files",
        *DATASETS,
    ]

    if args.tau_sensitivity:
        cmd.extend(
            [
                "--sensitivity",
                "--sensitivity-param",
                "tau",
                "--tau-values",
                args.tau_values,
            ]
        )

    print("Running command:")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()

