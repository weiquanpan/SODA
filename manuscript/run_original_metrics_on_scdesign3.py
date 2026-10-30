#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run original project metrics on newly generated scDesign3 datasets.

This script DOES NOT re-implement metrics. It reuses the existing benchmark
framework (registered methods + registered metrics) and writes consolidated
result tables for rebuttal use.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

from soda_benchmark.core.benchmark import SingleCellBatchBenchmark


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run original benchmark metrics on scDesign3 h5ad files."
    )
    parser.add_argument(
        "--data-dir",
        default="results/rebuttal/data",
        help="Directory containing scDesign3 h5ad files.",
    )
    parser.add_argument(
        "--pattern",
        default="scdesign3_*.h5ad",
        help="Glob pattern for selecting input h5ad files.",
    )
    parser.add_argument(
        "--output-dir",
        default="results/rebuttal_recalc",
        help="Directory to store benchmark outputs and summary tables.",
    )
    parser.add_argument("--n-bootstrap", type=int, default=200)
    parser.add_argument("--scvi-max-epochs", type=int, default=50)
    parser.add_argument(
        "--neighborhood-mode",
        choices=["expression", "representation"],
        default="representation",
    )
    return parser.parse_args()


def parse_scenario_seed(filename: str) -> Tuple[str, int]:
    # Expected pattern: scdesign3_<scenario>_seed<seed>.h5ad
    m = re.match(r"scdesign3_(.+)_seed(\d+)\.h5ad$", filename)
    if not m:
        return "unknown", -1
    return m.group(1), int(m.group(2))


def compute_overall_ranking(adata_dict, metrics) -> Dict[str, float]:
    scores = {}
    for method in adata_dict.keys():
        if method == "Original":
            continue
        total = 0.0
        n_counted = 0
        for mr in metrics.values():
            if method not in mr.scores.index:
                continue
            method_scores = mr.scores["score"]
            score_range = method_scores.max() - method_scores.min()
            if score_range <= 0:
                continue
            normalized = (method_scores[method] - method_scores.min()) / score_range
            if not mr.higher_is_better:
                normalized = 1 - normalized
            total += float(normalized)
            n_counted += 1
        if n_counted > 0:
            scores[method] = total / n_counted
    return scores


def format_mean_std(mean_val: float, std_val: float) -> str:
    if pd.isna(std_val):
        return f"{mean_val:.4f}"
    return f"{mean_val:.4f} ± {std_val:.4f}"


def main() -> None:
    args = parse_args()
    data_dir = Path(args.data_dir)
    out_root = Path(args.output_dir)
    run_dir = out_root / "benchmarks"
    out_root.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(data_dir.glob(args.pattern))
    if not files:
        raise FileNotFoundError(f"No files matched: {data_dir / args.pattern}")

    metric_rows: List[dict] = []
    rank_rows: List[dict] = []

    for f in files:
        scenario, seed = parse_scenario_seed(f.name)
        this_run_dir = run_dir / f"{scenario}_seed{seed}"
        this_run_dir.mkdir(parents=True, exist_ok=True)

        metric_params = {
            "n_bootstrap": args.n_bootstrap,
            "ci_level": 95,
            "random_state": 42,
            "neighborhood_mode": args.neighborhood_mode,
        }
        method_params = {
            "scVI": {"max_epochs": args.scvi_max_epochs},
        }
        bench = SingleCellBatchBenchmark(
            h5ad_path=str(f),
            method_params=method_params,
            metric_params=metric_params,
        )
        bench.run_full_benchmark(output_dir=str(this_run_dir), generate_plots=False)

        for metric_key, mr in bench.metrics.items():
            for method, row in mr.scores.iterrows():
                metric_rows.append(
                    {
                        "dataset": f.name,
                        "scenario": scenario,
                        "seed": seed,
                        "metric_key": metric_key,
                        "metric_name": mr.name,
                        "higher_is_better": mr.higher_is_better,
                        "method": method,
                        "score": float(row["score"]),
                        "ci_lower": float(row["ci_lower"]),
                        "ci_upper": float(row["ci_upper"]),
                    }
                )

        run_scores = compute_overall_ranking(bench.adata_dict, bench.metrics)
        for method, score in run_scores.items():
            rank_rows.append(
                {
                    "dataset": f.name,
                    "scenario": scenario,
                    "seed": seed,
                    "method": method,
                    "overall_score": score,
                }
            )

    metric_df = pd.DataFrame(metric_rows)
    rank_df = pd.DataFrame(rank_rows)

    metric_df.to_csv(out_root / "table_metrics_long.csv", index=False)
    rank_df.to_csv(out_root / "table_overall_ranking_long.csv", index=False)

    metric_summary = (
        metric_df.groupby(["scenario", "metric_key", "metric_name", "method"], as_index=False)
        .agg(score_mean=("score", "mean"), score_std=("score", "std"), n_runs=("score", "count"))
    )
    metric_summary["score_mean_std"] = metric_summary.apply(
        lambda x: format_mean_std(x["score_mean"], x["score_std"]), axis=1
    )
    metric_summary.to_csv(out_root / "table_metrics_summary.csv", index=False)

    rank_summary = (
        rank_df.groupby(["scenario", "method"], as_index=False)
        .agg(overall_mean=("overall_score", "mean"), overall_std=("overall_score", "std"), n_runs=("overall_score", "count"))
    )
    rank_summary["overall_mean_std"] = rank_summary.apply(
        lambda x: format_mean_std(x["overall_mean"], x["overall_std"]), axis=1
    )
    rank_summary = rank_summary.sort_values(["scenario", "overall_mean"], ascending=[True, False])
    rank_summary.to_csv(out_root / "table_overall_ranking_summary.csv", index=False)

    wide = metric_summary.pivot_table(
        index=["scenario", "method"],
        columns="metric_key",
        values="score_mean_std",
        aggfunc="first",
    ).reset_index()
    wide.to_csv(out_root / "table_rebuttal_main.csv", index=False)

    print("Done. Reused original metrics successfully.")
    print(f"Output directory: {out_root}")


if __name__ == "__main__":
    main()

