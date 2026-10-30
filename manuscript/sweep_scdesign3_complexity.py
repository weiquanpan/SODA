#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fixed-seed scDesign3 complexity sweep for nonlinear/local-distortion analysis.

What this script does
---------------------
1) Use ONE fixed seed to generate multiple simulated datasets with different
   scDesign3 formula settings (from mild to stronger local distortion).
2) Run the existing benchmark framework (original methods + original metrics).
3) Output comparison tables so you can inspect how method rankings (especially
   SODA) change under different complexity settings.
"""

from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import pandas as pd

from soda_benchmark.core.benchmark import SingleCellBatchBenchmark


@dataclass(frozen=True)
class SweepCase:
    name: str
    mu_formula: str
    sigma_formula: str
    other_covariates: str = "Method"


# Keep this list ordered from simple -> complex.
DEFAULT_CASES = [
    SweepCase(
        name="linear_batch_mean",
        mu_formula="CellType + Method",
        sigma_formula="1",
    ),
    SweepCase(
        name="local_mean_interaction",
        mu_formula="CellType * Method",
        sigma_formula="1",
    ),
    SweepCase(
        name="local_var_interaction",
        mu_formula="CellType + Method",
        sigma_formula="CellType * Method",
    ),
    SweepCase(
        name="local_mean_var_interaction",
        mu_formula="CellType * Method",
        sigma_formula="CellType * Method",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fixed-seed parameter sweep for scDesign3 complexity."
    )
    parser.add_argument("--seed", type=int, default=42, help="Single fixed seed.")
    parser.add_argument("--output-dir", default="results/complexity_sweep", help="Output root directory.")
    parser.add_argument("--nfeatures", type=int, default=800)
    parser.add_argument("--cells-per-group", type=int, default=200)
    parser.add_argument("--min-cells-per-group", type=int, default=120)
    parser.add_argument("--n-cores", type=int, default=4)
    parser.add_argument("--n-bootstrap", type=int, default=200)
    parser.add_argument("--scvi-max-epochs", type=int, default=50)
    parser.add_argument(
        "--copula",
        default="gaussian",
        help="scDesign3 copula setting, e.g. gaussian or vine (if supported by your scDesign3 version).",
    )
    parser.add_argument(
        "--neighborhood-mode",
        choices=["expression", "representation"],
        default="representation",
    )
    parser.add_argument(
        "--skip-generation",
        action="store_true",
        help="Skip data generation and only rerun benchmark from existing h5ad files.",
    )
    return parser.parse_args()


def run_cmd(cmd: List[str], cwd: Path) -> None:
    print("\n[RUN]", " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd), check=True)


def generate_dataset(
    root: Path,
    case: SweepCase,
    seed: int,
    data_dir: Path,
    nfeatures: int,
    cells_per_group: int,
    min_cells_per_group: int,
    n_cores: int,
    copula: str,
) -> Path:
    out_path = data_dir / f"scdesign3_{case.name}_seed{seed}.h5ad"
    cmd = [
        "Rscript",
        "simulate_sc3desigin_data.r",
        f"--output={out_path.as_posix()}",
        f"--seed={seed}",
        f"--nfeatures={nfeatures}",
        f"--cells-per-group={cells_per_group}",
        f"--min-cells-per-group={min_cells_per_group}",
        f"--n-cores={n_cores}",
        f"--mu-formula={case.mu_formula}",
        f"--sigma-formula={case.sigma_formula}",
        f"--copula={copula}",
        f"--other-covariates={case.other_covariates}",
        f"--condition-label={case.name}",
    ]
    run_cmd(cmd, cwd=root)
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


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parent
    out_root = root / args.output_dir / f"seed_{args.seed}"
    data_dir = out_root / "data"
    bench_dir = out_root / "benchmarks"
    out_root.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    bench_dir.mkdir(parents=True, exist_ok=True)

    metric_rows = []
    rank_rows = []

    for case in DEFAULT_CASES:
        h5ad_path = data_dir / f"scdesign3_{case.name}_seed{args.seed}.h5ad"
        if not args.skip_generation or not h5ad_path.exists():
            h5ad_path = generate_dataset(
                root=root,
                case=case,
                seed=args.seed,
                data_dir=data_dir,
                nfeatures=args.nfeatures,
                cells_per_group=args.cells_per_group,
                min_cells_per_group=args.min_cells_per_group,
                n_cores=args.n_cores,
                copula=args.copula,
            )

        run_dir = bench_dir / case.name
        run_dir.mkdir(parents=True, exist_ok=True)
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
            h5ad_path=str(h5ad_path),
            method_params=method_params,
            metric_params=metric_params,
        )
        bench.run_full_benchmark(output_dir=str(run_dir), generate_plots=False)

        for metric_key, mr in bench.metrics.items():
            for method, row in mr.scores.iterrows():
                metric_rows.append(
                    {
                        "case": case.name,
                        "mu_formula": case.mu_formula,
                        "sigma_formula": case.sigma_formula,
                        "metric_key": metric_key,
                        "metric_name": mr.name,
                        "higher_is_better": mr.higher_is_better,
                        "copula": args.copula,
                        "method": method,
                        "score": float(row["score"]),
                        "ci_lower": float(row["ci_lower"]),
                        "ci_upper": float(row["ci_upper"]),
                    }
                )

        overall = compute_overall_ranking(bench.adata_dict, bench.metrics)
        for method, score in overall.items():
            rank_rows.append(
                {
                    "case": case.name,
                    "mu_formula": case.mu_formula,
                    "sigma_formula": case.sigma_formula,
                    "copula": args.copula,
                    "method": method,
                    "overall_score": score,
                }
            )

    metric_df = pd.DataFrame(metric_rows)
    rank_df = pd.DataFrame(rank_rows)

    metric_long_csv = out_root / "table_metrics_long.csv"
    rank_long_csv = out_root / "table_overall_ranking_long.csv"
    metric_df.to_csv(metric_long_csv, index=False)
    rank_df.to_csv(rank_long_csv, index=False)

    rank_summary = rank_df.sort_values(["case", "overall_score"], ascending=[True, False]).copy()
    rank_summary["rank"] = rank_summary.groupby("case")["overall_score"].rank(
        method="min", ascending=False
    ).astype(int)
    rank_summary_csv = out_root / "table_overall_ranking_by_case.csv"
    rank_summary.to_csv(rank_summary_csv, index=False)

    # Focus table: how SODA rank/score changes across complexity cases.
    soda_focus = rank_summary[rank_summary["method"] == "SODA"][
        ["case", "mu_formula", "sigma_formula", "copula", "rank", "overall_score"]
    ].sort_values("case")
    soda_focus_csv = out_root / "table_soda_focus.csv"
    soda_focus.to_csv(soda_focus_csv, index=False)

    # Wide ranking table for quick visual comparison.
    wide_rank = rank_df.pivot_table(
        index="method",
        columns="case",
        values="overall_score",
        aggfunc="first",
    ).reset_index()
    wide_rank_csv = out_root / "table_overall_ranking_wide.csv"
    wide_rank.to_csv(wide_rank_csv, index=False)

    print("\nDone.")
    print(f"Seed: {args.seed}")
    print(f"Copula: {args.copula}")
    print("Cases:")
    for c in DEFAULT_CASES:
        print(f"  - {c.name}: mu={c.mu_formula}; sigma={c.sigma_formula}")
    print("\nGenerated:")
    print(f"  - {metric_long_csv}")
    print(f"  - {rank_long_csv}")
    print(f"  - {rank_summary_csv}")
    print(f"  - {wide_rank_csv}")
    print(f"  - {soda_focus_csv}")


if __name__ == "__main__":
    main()

