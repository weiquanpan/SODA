# -*- coding: utf-8 -*-
"""
Reviewer one-click sensitivity runner.

一次命令按顺序跑完 SODA 三组敏感性分析：
  1) tau
  2) alpha
  3) epsilon

固定数据集：
  - data/pbmc4k.h5ad
  - data/heart.h5ad
  - data/neuro.h5ad
"""

import argparse
import subprocess
import sys


DATASETS = ["pbmc4k.h5ad", "heart.h5ad", "neuro.h5ad"]


def run_one_sweep(param_name, values, neighborhood_mode, base_tau, base_alpha, base_epsilon):
    cmd = [
        sys.executable,
        "run_batch.py",
        "--sensitivity",
        "--sensitivity-param",
        param_name,
        "--neighborhood-mode",
        neighborhood_mode,
        "--soda-tau",
        str(base_tau),
        "--soda-alpha",
        str(base_alpha),
        "--soda-epsilon",
        str(base_epsilon),
        "--data-files",
        *DATASETS,
    ]

    if param_name == "tau":
        cmd.extend(["--tau-values", values])
    elif param_name == "alpha":
        cmd.extend(["--alpha-values", values])
    elif param_name == "epsilon":
        cmd.extend(["--epsilon-values", values])
    else:
        raise ValueError(f"Unsupported sensitivity parameter: {param_name}")

    print("\n" + "=" * 120)
    print(f"RUNNING {param_name.upper()} SENSITIVITY")
    print("=" * 120)
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description="Run tau/alpha/epsilon sensitivity on pbmc4k/heart/neuro."
    )
    parser.add_argument(
        "--neighborhood-mode",
        choices=["expression", "representation"],
        default="representation",
        help="Neighborhood metric mode.",
    )
    parser.add_argument(
        "--soda-tau",
        type=float,
        default=0.05,
        help="Base tau value used in non-tau sweeps.",
    )
    parser.add_argument(
        "--soda-alpha",
        type=float,
        default=0.7,
        help="Base alpha value used in non-alpha sweeps.",
    )
    parser.add_argument(
        "--soda-epsilon",
        type=float,
        default=0.05,
        help="Base epsilon value used in non-epsilon sweeps.",
    )
    parser.add_argument(
        "--tau-values",
        default="0.01,0.02,0.05,0.1,0.2",
        help="Comma-separated tau values.",
    )
    parser.add_argument(
        "--alpha-values",
        default="0.3,0.5,0.7,0.9",
        help="Comma-separated alpha values.",
    )
    parser.add_argument(
        "--epsilon-values",
        default="0.01,0.02,0.05,0.1,0.2",
        help="Comma-separated epsilon values.",
    )
    args = parser.parse_args()

    run_one_sweep(
        param_name="tau",
        values=args.tau_values,
        neighborhood_mode=args.neighborhood_mode,
        base_tau=args.soda_tau,
        base_alpha=args.soda_alpha,
        base_epsilon=args.soda_epsilon,
    )
    run_one_sweep(
        param_name="alpha",
        values=args.alpha_values,
        neighborhood_mode=args.neighborhood_mode,
        base_tau=args.soda_tau,
        base_alpha=args.soda_alpha,
        base_epsilon=args.soda_epsilon,
    )
    run_one_sweep(
        param_name="epsilon",
        values=args.epsilon_values,
        neighborhood_mode=args.neighborhood_mode,
        base_tau=args.soda_tau,
        base_alpha=args.soda_alpha,
        base_epsilon=args.soda_epsilon,
    )

    print("\nAll reviewer sensitivity sweeps completed successfully.")


if __name__ == "__main__":
    main()

