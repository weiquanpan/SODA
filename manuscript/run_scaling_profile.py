# -*- coding: utf-8 -*-
"""
Runtime + memory scaling benchmark for reviewer response.

This script benchmarks all available integration methods across increasing
dataset sizes and generates:
1) Raw per-run table
2) Aggregated table (mean/std by method and cell size)
3) Two-panel figure (runtime and memory vs #cells)

Usage examples
--------------
python run_scaling_profile.py --data-files pbmc4k.h5ad
python run_scaling_profile.py --data-files pbmc4k.h5ad pancreas.h5ad --cell-sizes 2000,4000,8000
"""

import os
import gc
import csv
import sys
import time
import argparse
import threading
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

from soda_benchmark.methods import get_registered_methods

try:
    import psutil
    HAS_PSUTIL = True
except Exception:
    HAS_PSUTIL = False


if sys.stdout and hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
if sys.stderr and hasattr(sys.stderr, "reconfigure"):
    try:
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


def parse_int_list(raw):
    if raw is None:
        return None
    values = [int(x.strip()) for x in raw.split(",") if x.strip()]
    values = sorted(set(v for v in values if v > 0))
    return values if values else None


def _total_rss_mb(process):
    rss = 0
    try:
        rss += process.memory_info().rss
    except Exception:
        return 0.0
    for child in process.children(recursive=True):
        try:
            rss += child.memory_info().rss
        except Exception:
            continue
    return rss / (1024.0 ** 2)


def run_with_resource_monitor(fn, sample_interval=0.05):
    """
    Run function and measure elapsed time + peak RSS delta.

    Returns
    -------
    result, elapsed_sec, peak_delta_mb, error_message
    """
    err = None
    result = None
    t0 = time.perf_counter()

    if not HAS_PSUTIL:
        try:
            result = fn()
        except Exception as e:
            err = str(e)
        elapsed = time.perf_counter() - t0
        return result, elapsed, None, err

    process = psutil.Process(os.getpid())
    baseline_mb = _total_rss_mb(process)
    peak_mb = baseline_mb
    stop_event = threading.Event()

    def sampler():
        nonlocal peak_mb
        while not stop_event.is_set():
            current = _total_rss_mb(process)
            if current > peak_mb:
                peak_mb = current
            time.sleep(sample_interval)

    thread = threading.Thread(target=sampler, daemon=True)
    thread.start()

    try:
        result = fn()
    except Exception as e:
        err = str(e)
    finally:
        stop_event.set()
        thread.join(timeout=1.0)

    elapsed = time.perf_counter() - t0
    peak_delta_mb = max(0.0, peak_mb - baseline_mb)
    return result, elapsed, peak_delta_mb, err


def choose_available_methods(selected_methods=None):
    registered = get_registered_methods()
    available = {}
    for name, method_cls in registered.items():
        ok, _ = method_cls.check_dependencies()
        if ok:
            available[name] = method_cls
    if selected_methods:
        available = {
            k: v for k, v in available.items() if k in set(selected_methods)
        }
    return available


def stratified_subsample_by_batch(adata, target_cells, random_state=42):
    if target_cells >= adata.n_obs:
        return adata.copy()
    if "batch" not in adata.obs:
        raise ValueError("Input AnnData must contain 'batch' in adata.obs.")

    rng = np.random.default_rng(random_state)
    batch_series = adata.obs["batch"].astype(str)
    batch_counts = batch_series.value_counts()
    batch_names = list(batch_counts.index)

    proportions = batch_counts.values / batch_counts.values.sum()
    alloc = np.floor(proportions * target_cells).astype(int)
    alloc = np.maximum(alloc, 1)

    while alloc.sum() > target_cells:
        idx = np.argmax(alloc)
        alloc[idx] -= 1
    while alloc.sum() < target_cells:
        idx = np.argmin(alloc)
        alloc[idx] += 1

    selected_idx = []
    for i, batch_name in enumerate(batch_names):
        idx_all = np.where(batch_series.values == batch_name)[0]
        take = min(int(alloc[i]), len(idx_all))
        pick = rng.choice(idx_all, size=take, replace=False)
        selected_idx.extend(pick.tolist())

    if len(selected_idx) < target_cells:
        remain = target_cells - len(selected_idx)
        all_idx = np.arange(adata.n_obs)
        mask = np.ones(adata.n_obs, dtype=bool)
        mask[selected_idx] = False
        extra_pool = all_idx[mask]
        if remain > 0 and len(extra_pool) > 0:
            extra = rng.choice(extra_pool, size=min(remain, len(extra_pool)), replace=False)
            selected_idx.extend(extra.tolist())

    selected_idx = np.array(selected_idx[:target_cells], dtype=int)
    selected_idx.sort()
    return adata[selected_idx].copy()


def benchmark_one_dataset(
    adata,
    dataset_name,
    cell_sizes,
    methods,
    method_params,
    repeats,
    random_state,
):
    records = []
    for n_cells in cell_sizes:
        if n_cells > adata.n_obs:
            continue

        for rep in range(repeats):
            ad_sub = stratified_subsample_by_batch(
                adata=adata,
                target_cells=n_cells,
                random_state=random_state + rep,
            )
            print(f"\n[Dataset={dataset_name}] cells={n_cells:,}, repeat={rep + 1}/{repeats}")

            for method_name, method_cls in methods.items():
                print(f"  - Running {method_name} ...")
                params = method_params.get(method_name, {})
                gc.collect()

                def _run():
                    instance = method_cls(**params)
                    return instance.apply(ad_sub.copy())

                result, elapsed_sec, peak_delta_mb, err = run_with_resource_monitor(_run)
                ok = result is not None and err is None
                rec = {
                    "dataset": dataset_name,
                    "n_cells": int(n_cells),
                    "repeat": int(rep + 1),
                    "method": method_name,
                    "runtime_sec": float(elapsed_sec),
                    "peak_rss_delta_mb": (
                        float(peak_delta_mb) if peak_delta_mb is not None else np.nan
                    ),
                    "status": "ok" if ok else "failed",
                    "error": "" if err is None else err,
                }
                records.append(rec)
                status_msg = "ok" if ok else "failed"
                mem_msg = (
                    f"{rec['peak_rss_delta_mb']:.2f} MB"
                    if not np.isnan(rec["peak_rss_delta_mb"]) else "N/A"
                )
                print(
                    f"    status={status_msg}, time={elapsed_sec:.2f}s, peak_mem={mem_msg}"
                )
    return records


def aggregate_records(raw_df):
    ok_df = raw_df[raw_df["status"] == "ok"].copy()
    agg = (
        ok_df.groupby(["method", "n_cells"], as_index=False)
        .agg(
            runtime_sec_mean=("runtime_sec", "mean"),
            runtime_sec_std=("runtime_sec", "std"),
            peak_rss_delta_mb_mean=("peak_rss_delta_mb", "mean"),
            peak_rss_delta_mb_std=("peak_rss_delta_mb", "std"),
            n_runs=("runtime_sec", "count"),
        )
        .sort_values(["method", "n_cells"])
    )
    return agg


def plot_scaling_figure(agg_df, output_path):
    if agg_df.empty:
        print("No successful runs. Skip plotting.")
        return

    sns.set_style("whitegrid")
    methods = sorted(agg_df["method"].unique())
    palette = sns.color_palette("tab10", n_colors=len(methods))
    color_map = dict(zip(methods, palette))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharex=False)
    ax_time, ax_mem = axes

    for method in methods:
        sub = agg_df[agg_df["method"] == method].sort_values("n_cells")
        x = sub["n_cells"].values

        y_time = sub["runtime_sec_mean"].values / 60.0
        y_time_err = np.nan_to_num(sub["runtime_sec_std"].values / 60.0, nan=0.0)
        ax_time.plot(x, y_time, marker="o", lw=2, label=method, color=color_map[method])
        ax_time.fill_between(
            x,
            np.maximum(0, y_time - y_time_err),
            y_time + y_time_err,
            alpha=0.18,
            color=color_map[method],
        )

        y_mem = sub["peak_rss_delta_mb_mean"].values / 1024.0
        y_mem_err = np.nan_to_num(sub["peak_rss_delta_mb_std"].values / 1024.0, nan=0.0)
        ax_mem.plot(x, y_mem, marker="o", lw=2, label=method, color=color_map[method])
        ax_mem.fill_between(
            x,
            np.maximum(0, y_mem - y_mem_err),
            y_mem + y_mem_err,
            alpha=0.18,
            color=color_map[method],
        )

    ax_time.set_xlabel("# Cells")
    ax_time.set_ylabel("Runtime (minutes)")
    ax_time.set_title("a  Runtime Scaling")

    ax_mem.set_xlabel("# Cells")
    ax_mem.set_ylabel("Peak RSS Delta (GB)")
    ax_mem.set_title("b  Memory Scaling")

    handles, labels = ax_time.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=min(6, len(methods)), frameon=False)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {output_path}")


def resolve_data_paths(data_files):
    if not data_files:
        data_dir = "data"
        if not os.path.exists(data_dir):
            return []
        files = [
            os.path.join(data_dir, f)
            for f in os.listdir(data_dir)
            if f.endswith(".h5ad")
        ]
        return sorted(files)

    resolved = []
    for item in data_files:
        if os.path.exists(item):
            resolved.append(item)
            continue
        candidate = os.path.join("data", item)
        if os.path.exists(candidate):
            resolved.append(candidate)
    return sorted(resolved)


def infer_cell_sizes(adatas, explicit_sizes):
    if explicit_sizes:
        return explicit_sizes
    max_cells = min(ad.n_obs for ad in adatas)
    # Default increasing sizes for one figure line that can be shared by all datasets.
    fractions = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0]
    sizes = sorted(set(max(500, int(max_cells * f)) for f in fractions))
    sizes = [s for s in sizes if s <= max_cells]
    if not sizes:
        sizes = [max_cells]
    return sizes


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark runtime and memory scaling across increasing dataset sizes."
    )
    parser.add_argument(
        "--data-files",
        nargs="*",
        default=None,
        help=(
            "Input h5ad files. Can be absolute/relative paths, or file names under data/. "
            "If omitted, all *.h5ad under data/ are used."
        ),
    )
    parser.add_argument(
        "--cell-sizes",
        default=None,
        help=(
            "Comma-separated target cell sizes, e.g. 2000,5000,10000. "
            "If omitted, sizes are auto-inferred from smallest dataset."
        ),
    )
    parser.add_argument(
        "--methods",
        nargs="*",
        default=None,
        help="Optional method names to run. If omitted, run all available methods.",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=1,
        help="Number of repeats per (dataset, size, method). Default: 1",
    )
    parser.add_argument(
        "--output-root",
        default="results/runtime_memory_scaling",
        help="Output root directory.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed for subsampling.",
    )
    parser.add_argument(
        "--soda-tau",
        type=float,
        default=0.05,
        help="SODA tau/divergence threshold.",
    )
    parser.add_argument(
        "--soda-alpha",
        type=float,
        default=0.7,
        help="SODA alpha.",
    )
    parser.add_argument(
        "--soda-epsilon",
        type=float,
        default=0.05,
        help="SODA epsilon.",
    )
    args = parser.parse_args()

    data_paths = resolve_data_paths(args.data_files)
    if not data_paths:
        print("No input h5ad files found. Please provide --data-files or put files into data/.")
        return

    methods = choose_available_methods(args.methods)
    if not methods:
        print("No available methods to run (dependency check failed or invalid --methods).")
        return

    print("Methods to run:", ", ".join(methods.keys()))
    print("Memory monitor:", "psutil RSS sampler" if HAS_PSUTIL else "disabled (psutil not installed)")

    adatas = []
    meta = []
    for path in data_paths:
        adata = sc.read_h5ad(path)
        if "batch" not in adata.obs:
            print(f"Skip {path}: missing 'batch' in adata.obs.")
            continue
        adatas.append(adata)
        meta.append((os.path.splitext(os.path.basename(path))[0], path, adata.n_obs, adata.n_vars))

    if not adatas:
        print("No valid datasets (all missing required columns).")
        return

    explicit_sizes = parse_int_list(args.cell_sizes)
    cell_sizes = infer_cell_sizes(adatas, explicit_sizes)

    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = os.path.join(args.output_root, f"run_{ts}")
    os.makedirs(out_dir, exist_ok=True)

    method_params = {
        "SODA": {
            "tau": args.soda_tau,
            "alpha": args.soda_alpha,
            "epsilon": args.soda_epsilon,
        }
    }

    print("\nDatasets:")
    for name, path, n_obs, n_vars in meta:
        print(f"  - {name}: cells={n_obs:,}, genes={n_vars:,}, file={path}")
    print("Target cell sizes:", ", ".join(f"{x:,}" for x in cell_sizes))

    records = []
    for (name, _path, _n_obs, _n_vars), adata in zip(meta, adatas):
        dataset_sizes = [s for s in cell_sizes if s <= adata.n_obs]
        if not dataset_sizes:
            print(f"Skip dataset={name}: all target sizes are larger than dataset size.")
            continue
        records.extend(
            benchmark_one_dataset(
                adata=adata,
                dataset_name=name,
                cell_sizes=dataset_sizes,
                methods=methods,
                method_params=method_params,
                repeats=args.repeats,
                random_state=args.random_state,
            )
        )

    if not records:
        print("No benchmark records generated.")
        return

    raw_df = pd.DataFrame(records)
    raw_csv = os.path.join(out_dir, "runtime_memory_raw.csv")
    raw_df.to_csv(raw_csv, index=False, encoding="utf-8")

    agg_df = aggregate_records(raw_df)
    agg_csv = os.path.join(out_dir, "runtime_memory_aggregated.csv")
    agg_df.to_csv(agg_csv, index=False, encoding="utf-8")

    fig_path = os.path.join(out_dir, "figure_runtime_memory_scaling.png")
    plot_scaling_figure(agg_df, fig_path)

    summary_txt = os.path.join(out_dir, "summary.txt")
    with open(summary_txt, "w", encoding="utf-8") as f:
        f.write("Runtime + Memory Scaling Benchmark\n")
        f.write("=" * 60 + "\n")
        f.write(f"Generated at: {datetime.now().isoformat()}\n")
        f.write(f"Datasets: {len(meta)}\n")
        f.write(f"Methods: {', '.join(methods.keys())}\n")
        f.write(f"Cell sizes: {cell_sizes}\n")
        f.write(f"Repeats: {args.repeats}\n")
        f.write(f"Memory monitor: {'psutil RSS sampler' if HAS_PSUTIL else 'not available'}\n")
        f.write("\nOutput files:\n")
        f.write(f"- {raw_csv}\n")
        f.write(f"- {agg_csv}\n")
        f.write(f"- {fig_path}\n")

    print("\nDone.")
    print(f"Output directory: {os.path.abspath(out_dir)}")
    print(f"  - Raw table: {raw_csv}")
    print(f"  - Aggregated table: {agg_csv}")
    print(f"  - Figure: {fig_path}")
    print(f"  - Summary: {summary_txt}")


if __name__ == "__main__":
    main()
