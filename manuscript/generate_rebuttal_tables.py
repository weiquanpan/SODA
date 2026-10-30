#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate rebuttal-ready tables and narrative for reviewer comments on batch correction.

Workflow
--------
1) Generate scDesign3 datasets under two scenarios:
   - no_batch:       mu_formula = "CellType" (calibration-like)
   - technical_batch: mu_formula = "CellType + Method" (true correction setting)
2) Run benchmark methods on each dataset (multiple seeds).
3) Aggregate metrics into CSV tables and create a markdown explanation draft.
"""

from __future__ import annotations

import argparse
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import pandas as pd

from soda_benchmark.core.benchmark import SingleCellBatchBenchmark


@dataclass(frozen=True)
class Scenario:
    name: str
    mu_formula: str
    other_covariates: str
    condition_label: str


SCENARIOS = [
    Scenario(
        name="no_batch",
        mu_formula="CellType",
        other_covariates="Method",
        condition_label="no_batch",
    ),
    Scenario(
        name="technical_batch",
        mu_formula="CellType + Method",
        other_covariates="Method",
        condition_label="technical_batch",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate rebuttal tables for scDesign3-based reviewer response."
    )
    parser.add_argument("--seeds", default="41,42,43", help="Comma-separated seeds.")
    parser.add_argument("--output-dir", default="results/rebuttal", help="Output root directory.")
    parser.add_argument("--nfeatures", type=int, default=800)
    parser.add_argument("--cells-per-group", type=int, default=200)
    parser.add_argument("--min-cells-per-group", type=int, default=120)
    parser.add_argument("--n-cores", type=int, default=4)
    parser.add_argument("--n-bootstrap", type=int, default=200)
    parser.add_argument("--scvi-max-epochs", type=int, default=50)
    parser.add_argument(
        "--copula",
        default="gaussian",
        help="scDesign3 copula setting, e.g. gaussian or vine (if supported).",
    )
    parser.add_argument(
        "--skip-generation",
        action="store_true",
        help="Skip data generation and only run benchmark aggregation.",
    )
    return parser.parse_args()


def parse_seeds(text: str) -> List[int]:
    seeds = []
    for x in text.split(","):
        x = x.strip()
        if not x:
            continue
        seeds.append(int(x))
    if not seeds:
        raise ValueError("No valid seeds provided.")
    return seeds


def run_cmd(cmd: List[str], cwd: Path) -> None:
    print("\n[RUN]", " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd), check=True)


def generate_dataset(
    project_root: Path,
    scenario: Scenario,
    seed: int,
    data_dir: Path,
    nfeatures: int,
    cells_per_group: int,
    min_cells_per_group: int,
    n_cores: int,
    copula: str,
) -> Path:
    out_path = data_dir / f"scdesign3_{scenario.name}_seed{seed}.h5ad"
    cmd = [
        "Rscript",
        "simulate_sc3desigin_data.r",
        f"--output={out_path.as_posix()}",
        f"--seed={seed}",
        f"--nfeatures={nfeatures}",
        f"--cells-per-group={cells_per_group}",
        f"--min-cells-per-group={min_cells_per_group}",
        f"--n-cores={n_cores}",
        f"--mu-formula={scenario.mu_formula}",
        f"--copula={copula}",
        f"--other-covariates={scenario.other_covariates}",
        f"--condition-label={scenario.condition_label}",
    ]
    run_cmd(cmd, cwd=project_root)
    return out_path


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


def run_single_benchmark(
    h5ad_path: Path,
    output_dir: Path,
    n_bootstrap: int,
    scvi_max_epochs: int,
):
    metric_params = {
        "n_bootstrap": n_bootstrap,
        "ci_level": 95,
        "random_state": 42,
        "neighborhood_mode": "representation",
    }
    method_params = {
        "scVI": {"max_epochs": scvi_max_epochs},
    }
    benchmark = SingleCellBatchBenchmark(
        h5ad_path=str(h5ad_path),
        method_params=method_params,
        metric_params=metric_params,
    )
    benchmark.run_full_benchmark(output_dir=str(output_dir), generate_plots=False)
    return benchmark


def format_mean_std(mean_val: float, std_val: float) -> str:
    if pd.isna(std_val):
        return f"{mean_val:.4f}"
    return f"{mean_val:.4f} ± {std_val:.4f}"


def main() -> None:
    args = parse_args()
    seeds = parse_seeds(args.seeds)
    project_root = Path(__file__).resolve().parent
    out_root = project_root / args.output_dir
    data_dir = out_root / "data"
    bench_dir = out_root / "benchmarks"
    out_root.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    bench_dir.mkdir(parents=True, exist_ok=True)

    all_metric_rows = []
    all_rank_rows = []

    for scenario in SCENARIOS:
        for seed in seeds:
            h5ad_path = data_dir / f"scdesign3_{scenario.name}_seed{seed}.h5ad"
            if not args.skip_generation or not h5ad_path.exists():
                h5ad_path = generate_dataset(
                    project_root=project_root,
                    scenario=scenario,
                    seed=seed,
                    data_dir=data_dir,
                    nfeatures=args.nfeatures,
                    cells_per_group=args.cells_per_group,
                    min_cells_per_group=args.min_cells_per_group,
                    n_cores=args.n_cores,
                    copula=args.copula,
                )

            run_dir = bench_dir / f"{scenario.name}_seed{seed}"
            run_dir.mkdir(parents=True, exist_ok=True)
            bench = run_single_benchmark(
                h5ad_path=h5ad_path,
                output_dir=run_dir,
                n_bootstrap=args.n_bootstrap,
                scvi_max_epochs=args.scvi_max_epochs,
            )

            # Long-form metric table.
            for metric_key, mr in bench.metrics.items():
                for method, row in mr.scores.iterrows():
                    all_metric_rows.append(
                        {
                            "scenario": scenario.name,
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

            # Per-run overall ranking.
            run_scores = compute_overall_ranking(bench.adata_dict, bench.metrics)
            for method, overall_score in run_scores.items():
                all_rank_rows.append(
                    {
                        "scenario": scenario.name,
                        "seed": seed,
                        "method": method,
                        "overall_score": overall_score,
                    }
                )

    metric_df = pd.DataFrame(all_metric_rows)
    rank_df = pd.DataFrame(all_rank_rows)

    metric_csv = out_root / "table_metrics_long.csv"
    rank_csv = out_root / "table_overall_ranking_long.csv"
    metric_df.to_csv(metric_csv, index=False)
    rank_df.to_csv(rank_csv, index=False)

    # Summary tables.
    metric_summary = (
        metric_df.groupby(["scenario", "metric_key", "metric_name", "method"], as_index=False)
        .agg(score_mean=("score", "mean"), score_std=("score", "std"), n_runs=("score", "count"))
    )
    metric_summary["score_mean_std"] = metric_summary.apply(
        lambda x: format_mean_std(x["score_mean"], x["score_std"]), axis=1
    )
    metric_summary_csv = out_root / "table_metrics_summary.csv"
    metric_summary.to_csv(metric_summary_csv, index=False)

    rank_summary = (
        rank_df.groupby(["scenario", "method"], as_index=False)
        .agg(overall_mean=("overall_score", "mean"), overall_std=("overall_score", "std"), n_runs=("overall_score", "count"))
    )
    rank_summary["overall_mean_std"] = rank_summary.apply(
        lambda x: format_mean_std(x["overall_mean"], x["overall_std"]), axis=1
    )
    rank_summary = rank_summary.sort_values(["scenario", "overall_mean"], ascending=[True, False])
    rank_summary_csv = out_root / "table_overall_ranking_summary.csv"
    rank_summary.to_csv(rank_summary_csv, index=False)

    # Wide table for manuscript/rebuttal copy-paste.
    wide = metric_summary.pivot_table(
        index=["scenario", "method"],
        columns="metric_key",
        values="score_mean_std",
        aggfunc="first",
    ).reset_index()
    wide_csv = out_root / "table_rebuttal_main.csv"
    wide.to_csv(wide_csv, index=False)

    # Draft explanation.
    top_by_scenario = (
        rank_summary.sort_values(["scenario", "overall_mean"], ascending=[True, False])
        .groupby("scenario", as_index=False)
        .first()
    )
    top_map = {row["scenario"]: (row["method"], row["overall_mean_std"]) for _, row in top_by_scenario.iterrows()}

    md_path = out_root / "rebuttal_explanation_draft.md"
    lines = [
        "# Reviewer Response Draft (scDesign3-based)",
        "",
        "## What we changed",
        "- We replaced pseudo-batch random splitting with scDesign3 simulation using real technical covariates (`Method`).",
        "- We evaluated both calibration-like and true correction settings:",
        "  - `no_batch`: `mu_formula = CellType`",
        "  - `technical_batch`: `mu_formula = CellType + Method`",
        f"- We repeated experiments across seeds: {', '.join(map(str, seeds))}.",
        "",
        "## Main evidence tables",
        f"- Long metrics table: `{metric_csv.as_posix()}`",
        f"- Metric summary (mean/std): `{metric_summary_csv.as_posix()}`",
        f"- Overall ranking (long): `{rank_csv.as_posix()}`",
        f"- Overall ranking summary: `{rank_summary_csv.as_posix()}`",
        f"- Rebuttal main table (wide): `{wide_csv.as_posix()}`",
        "",
        "## Suggested wording",
        "We agree that pseudo-batch splitting primarily evaluates calibration, not true correction.",
        "To address this concern, we generated datasets with scDesign3 under two controlled settings:",
        "(i) no technical batch effect (`mu_formula = CellType`) and",
        "(ii) technical batch effect (`mu_formula = CellType + Method`, with real platform covariate `Method`).",
        "Across multiple random seeds, we benchmarked all available integration methods and summarized both per-metric",
        "performance and overall normalized rankings (Tables generated above).",
        "",
        f"- Top method in `no_batch`: {top_map.get('no_batch', ('N/A', 'N/A'))[0]} ({top_map.get('no_batch', ('N/A', 'N/A'))[1]})",
        f"- Top method in `technical_batch`: {top_map.get('technical_batch', ('N/A', 'N/A'))[0]} ({top_map.get('technical_batch', ('N/A', 'N/A'))[1]})",
        "",
        "These results directly evaluate correction under real technical covariates rather than pseudo-batch calibration.",
    ]
    md_path.write_text("\n".join(lines), encoding="utf-8")

    print("\nDone.")
    print("Generated files:")
    print(f"  - {metric_csv}")
    print(f"  - {metric_summary_csv}")
    print(f"  - {rank_csv}")
    print(f"  - {rank_summary_csv}")
    print(f"  - {wide_csv}")
    print(f"  - {md_path}")


if __name__ == "__main__":
    main()

