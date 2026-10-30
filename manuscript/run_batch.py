# -*- coding: utf-8 -*-
"""
批量运行脚本：对 data/ 目录下多个 H5AD 文件依次执行基准流程，
并在 results/<dataset_name>/ 下保存对应输出。

用法
----
    python run_batch.py
"""

import os
import argparse
import csv
import gc
import importlib
from datetime import datetime
import matplotlib.pyplot as plt

from soda_benchmark import SingleCellBatchBenchmark
from soda_benchmark.metrics import get_registered_metrics


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


def release_memory():
    """在关键阶段主动回收 Python/R/深度学习缓存内存。"""
    gc.collect()
    try:
        import matplotlib.pyplot as _plt
        _plt.close("all")
    except Exception:
        pass
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    except Exception:
        pass
    try:
        tf = importlib.import_module("tensorflow")
        tf.keras.backend.clear_session()
    except Exception:
        pass
    try:
        import rpy2.robjects as ro
        ro.r("gc()")
    except Exception:
        pass
    gc.collect()


def _select_metric_keys(metrics, metric_keys=None):
    if metric_keys is None:
        return list(metrics.keys())
    return [k for k in metric_keys if k in metrics]


def compute_overall_ranking(adata_dict, metrics, metric_keys=None, include_original=False):
    """按指定指标计算综合排名。默认不将 Original 纳入方法排行榜。"""
    scores = {}
    selected_metric_keys = _select_metric_keys(metrics, metric_keys)
    if not selected_metric_keys:
        return []

    for method in adata_dict.keys():
        if not include_original and method == "Original":
            continue
        total_score = 0.0
        n_counted = 0
        for key in selected_metric_keys:
            mr = metrics[key]
            if method not in mr.scores.index:
                continue

            method_scores = mr.scores["score"]
            score_range = method_scores.max() - method_scores.min()
            if score_range <= 0:
                continue

            normalized = (method_scores[method] - method_scores.min()) / score_range
            if not mr.higher_is_better:
                normalized = 1 - normalized
            total_score += float(normalized)
            n_counted += 1

        if n_counted > 0:
            scores[method] = total_score / n_counted

    return sorted(scores.items(), key=lambda x: x[1], reverse=True)


def extract_method_overall_score(ranking, method_name="SODA"):
    """从综合排名中提取指定方法的分数和名次。"""
    for rank, (method, score) in enumerate(ranking, 1):
        if method == method_name:
            return float(score), rank
    return None, None


def extract_method_metric_scores(metrics, method_name="SODA"):
    """提取指定方法在各指标下的原始分数。"""
    result = {}
    for metric_name, metric_result in metrics.items():
        score = None
        if (
            hasattr(metric_result, "scores")
            and method_name in metric_result.scores.index
            and "score" in metric_result.scores.columns
        ):
            score = float(metric_result.scores.loc[method_name, "score"])
        result[metric_name] = score
    return result


def parse_float_list(raw):
    return [float(x.strip()) for x in raw.split(",") if x.strip()]


def calculate_metrics_for_adata_dict(adata_original, adata_dict, metric_params):
    """给定 adata_dict 重新计算全套评估指标。"""
    metrics = {}
    registered = get_registered_metrics()
    for key, metric_cls in registered.items():
        metric_instance = metric_cls(**metric_params)
        result = metric_instance.calculate(
            adata_original=adata_original,
            adata_dict=adata_dict,
        )
        metrics[key] = result
    return metrics


def prepare_non_soda_cache(h5ad_path, neighborhood_mode, base_soda_params,
                           batch_key="batch", batch_label_map=None,
                           recompute_original_pca=False, label_key=None,
                           lisi_n_neighbors=30):
    """每个数据集只运行一次非 SODA 方法，返回可复用结果。"""
    metric_params = {
        "n_bootstrap": 1000,
        "ci_level": 95,
        "random_state": 42,
        "neighborhood_mode": neighborhood_mode,
        "label_key": label_key,
        "lisi_n_neighbors": lisi_n_neighbors,
    }
    method_params = {"SODA": base_soda_params}
    probe = SingleCellBatchBenchmark(
        h5ad_path=h5ad_path,
        method_params=method_params,
        metric_params=metric_params,
        batch_key=batch_key,
        batch_label_map=batch_label_map,
        recompute_original_pca=recompute_original_pca,
    )
    non_soda_methods = [name for name in probe.available_methods.keys() if name != "SODA"]

    benchmark = SingleCellBatchBenchmark(
        h5ad_path=h5ad_path,
        method_params=method_params,
        metric_params=metric_params,
        selected_methods=non_soda_methods,
        batch_key=batch_key,
        batch_label_map=batch_label_map,
        recompute_original_pca=recompute_original_pca,
    )
    benchmark.load_preprocessed_data()
    benchmark._apply_all_methods()
    return benchmark.adata, benchmark.adata_dict, metric_params


def run_soda_only(h5ad_path, neighborhood_mode, soda_params,
                  batch_key="batch", batch_label_map=None,
                  recompute_original_pca=False, label_key=None,
                  lisi_n_neighbors=30):
    """仅运行一次 SODA，供敏感性扫描复用。"""
    metric_params = {
        "n_bootstrap": 1000,
        "ci_level": 95,
        "random_state": 42,
        "neighborhood_mode": neighborhood_mode,
        "label_key": label_key,
        "lisi_n_neighbors": lisi_n_neighbors,
    }
    method_params = {"SODA": soda_params}
    benchmark = SingleCellBatchBenchmark(
        h5ad_path=h5ad_path,
        method_params=method_params,
        metric_params=metric_params,
        selected_methods=["SODA"],
        batch_key=batch_key,
        batch_label_map=batch_label_map,
        recompute_original_pca=recompute_original_pca,
    )
    benchmark.load_preprocessed_data()
    benchmark._apply_all_methods()
    return benchmark.adata_dict.get("SODA")


def write_batch_ranking_summary(batch_results, output_root, neighborhood_mode):
    """写出批量任务总榜文件（每个数据集 + 跨数据集平均）。"""
    os.makedirs(output_root, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    save_path = os.path.join(output_root, f"batch_overall_ranking_{ts}.txt")

    with open(save_path, "w", encoding="utf-8") as f:
        f.write("=" * 100 + "\n")
        f.write("BATCH BENCHMARK OVERALL RANKING SUMMARY\n")
        f.write("=" * 100 + "\n")
        f.write(f"Neighborhood mode: {neighborhood_mode}\n")
        f.write(f"Datasets: {len(batch_results)}\n\n")

        method_to_scores = {}

        for item in batch_results:
            dataset_name = item["dataset_name"]
            ranking = item["ranking"]
            ranking_with_original = item.get("ranking_with_original", [])
            batch_ranking = item.get("batch_removal_ranking", [])
            bio_ranking = item.get("bio_conservation_ranking", [])
            report_path = os.path.join(item["output_dir"], "benchmark_report.txt")

            f.write(f"DATASET: {dataset_name}\n")
            f.write("-" * 100 + "\n")
            baseline_map = {m: s for m, s in ranking_with_original}
            if "Original" in baseline_map:
                f.write(
                    f"  Baseline (Original, integration score): "
                    f"{baseline_map['Original']:.4f}\n"
                )
            f.write("  [Integration ranking | methods only]\n")
            for rank, (method, score) in enumerate(ranking, 1):
                stars = "★" * min(5, int(score * 5) + 1)
                f.write(f"  {rank}. {method:20s}: {score:.4f}  {stars}\n")
                method_to_scores.setdefault(method, []).append(score)
            if ranking:
                f.write(f"\n  ★★★ TOP PERFORMER: {ranking[0][0]} ★★★\n")
            if batch_ranking:
                f.write(
                    f"  Top (batch removal): {batch_ranking[0][0]} "
                    f"({batch_ranking[0][1]:.4f})\n"
                )
            if bio_ranking:
                f.write(
                    f"  Top (bio conservation): {bio_ranking[0][0]} "
                    f"({bio_ranking[0][1]:.4f})\n"
                )
            f.write(f"  Report: {report_path}\n\n")

        f.write("=" * 100 + "\n")
        f.write("GLOBAL AVERAGE RANKING ACROSS DATASETS\n")
        f.write("-" * 100 + "\n")

        avg_scores = []
        for method, values in method_to_scores.items():
            avg_scores.append((method, sum(values) / len(values), len(values)))
        avg_scores.sort(key=lambda x: x[1], reverse=True)

        for rank, (method, avg_score, n) in enumerate(avg_scores, 1):
            stars = "★" * min(5, int(avg_score * 5) + 1)
            f.write(f"  {rank}. {method:20s}: {avg_score:.4f}  {stars}  (n={n})\n")

        if avg_scores:
            f.write(f"\n  ★★★ GLOBAL TOP PERFORMER: {avg_scores[0][0]} ★★★\n")

    print(f"\n✓ Batch ranking summary saved to: {save_path}")
    return save_path


def write_batch_gain_vs_original_summary(batch_results, output_root):
    """按数据集与跨数据集写出相对 Original baseline 的综合分提升。"""
    os.makedirs(output_root, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    save_path = os.path.join(output_root, f"batch_gain_vs_original_{ts}.txt")

    baseline_name = "Original"
    with open(save_path, "w", encoding="utf-8") as f:
        f.write("=" * 100 + "\n")
        f.write("BATCH RELATIVE GAIN VS ORIGINAL BASELINE (OVERALL SCORE)\n")
        f.write("=" * 100 + "\n\n")

        method_to_improvements = {}
        for item in batch_results:
            dataset_name = item["dataset_name"]
            ranking = item["ranking"]
            ranking_with_original = item.get("ranking_with_original", [])
            rank_map = {m: s for m, s in ranking_with_original}
            if baseline_name not in rank_map:
                continue
            baseline_score = rank_map[baseline_name]

            f.write(f"DATASET: {dataset_name}\n")
            f.write("-" * 100 + "\n")
            f.write(f"  Baseline ({baseline_name}) overall score: {baseline_score:.4f}\n")
            rows = []
            for method, score in ranking:
                if method == baseline_name:
                    continue
                improvement = score - baseline_score
                if abs(baseline_score) > 1e-12:
                    pct = (improvement / abs(baseline_score)) * 100.0
                    pct_text = f"{pct:+.2f}%"
                else:
                    pct_text = "n/a"
                rows.append((method, score, improvement, pct_text))
                method_to_improvements.setdefault(method, []).append(improvement)

            rows.sort(key=lambda x: x[2], reverse=True)
            for method, score, improvement, pct_text in rows:
                f.write(
                    f"  - {method:20s}: score={score:.4f}, "
                    f"improvement={improvement:+.4f} ({pct_text})\n"
                )
            f.write("\n")

        f.write("=" * 100 + "\n")
        f.write("GLOBAL MEAN IMPROVEMENT VS ORIGINAL\n")
        f.write("-" * 100 + "\n")
        avg_rows = []
        for method, improvements in method_to_improvements.items():
            if not improvements:
                continue
            avg_rows.append((method, sum(improvements) / len(improvements), len(improvements)))
        avg_rows.sort(key=lambda x: x[1], reverse=True)
        for rank, (method, mean_imp, n) in enumerate(avg_rows, 1):
            f.write(f"  {rank}. {method:20s}: {mean_imp:+.4f}  (n={n})\n")

    print(f"\n✓ Batch gain-vs-original summary saved to: {save_path}")
    return save_path


def run_for_dataset(h5ad_path, output_root="results",
                    generate_plots=False, neighborhood_mode="expression",
                    soda_params=None, selected_methods=None,
                    batch_key="batch", batch_label_map=None,
                    recompute_original_pca=False, label_key=None,
                    lisi_n_neighbors=30):
    dataset_name = os.path.splitext(os.path.basename(h5ad_path))[0]
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    run_folder = with_v1_prefix(dataset_name + "_" + timestamp)
    output_dir = os.path.join(output_root, run_folder)
    os.makedirs(output_dir, exist_ok=True)

    print("\n" + "=" * 100)
    print(f"DATASET: {dataset_name}")
    print("=" * 100)

    metric_params = {
        "n_bootstrap": 1000,
        "ci_level": 95,
        "random_state": 42,
        "neighborhood_mode": neighborhood_mode,
        "label_key": label_key,
        "lisi_n_neighbors": lisi_n_neighbors,
    }
    method_params = {"SODA": soda_params} if soda_params else None
    benchmark = SingleCellBatchBenchmark(
        h5ad_path=h5ad_path,
        method_params=method_params,
        metric_params=metric_params,
        selected_methods=selected_methods,
        batch_key=batch_key,
        batch_label_map=batch_label_map,
        recompute_original_pca=recompute_original_pca,
    )
    benchmark.run_full_benchmark(
        output_dir=output_dir,
        generate_plots=generate_plots,
    )
    integration_metric_keys = ["ilisi", "clisi", "asw_batch"]
    batch_removal_metric_keys = ["ilisi", "asw_batch"]
    bio_conservation_metric_keys = ["clisi"]

    ranking = compute_overall_ranking(
        benchmark.adata_dict,
        benchmark.metrics,
        metric_keys=integration_metric_keys,
        include_original=False,
    )
    ranking_with_original = compute_overall_ranking(
        benchmark.adata_dict,
        benchmark.metrics,
        metric_keys=integration_metric_keys,
        include_original=True,
    )
    batch_removal_ranking = compute_overall_ranking(
        benchmark.adata_dict,
        benchmark.metrics,
        metric_keys=batch_removal_metric_keys,
        include_original=False,
    )
    bio_conservation_ranking = compute_overall_ranking(
        benchmark.adata_dict,
        benchmark.metrics,
        metric_keys=bio_conservation_metric_keys,
        include_original=False,
    )
    soda_overall_score, soda_rank = extract_method_overall_score(
        ranking, method_name="SODA"
    )
    soda_metric_scores = extract_method_metric_scores(
        benchmark.metrics, method_name="SODA"
    )
    result = {
        "dataset_name": dataset_name,
        "output_dir": output_dir,
        "ranking": ranking,
        "ranking_with_original": ranking_with_original,
        "batch_removal_ranking": batch_removal_ranking,
        "bio_conservation_ranking": bio_conservation_ranking,
        "soda_overall_score": soda_overall_score,
        "soda_rank": soda_rank,
        "soda_metric_scores": soda_metric_scores,
    }
    del benchmark
    release_memory()
    return result


def run_sensitivity_sweep(h5_files, output_root, neighborhood_mode,
                          sweep_param, sweep_values, base_soda_params,
                          generate_sensitivity_plots=False,
                          batch_key="batch", batch_label_map=None,
                          recompute_original_pca=False, label_key=None,
                          lisi_n_neighbors=30):
    """对单个 SODA 超参数执行批量敏感性分析。"""
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    sensitivity_folder = with_v1_prefix(f"sensitivity_{sweep_param}_{ts}")
    sensitivity_root = os.path.join(output_root, sensitivity_folder)
    os.makedirs(sensitivity_root, exist_ok=True)

    records = []
    for h5ad_path in h5_files:
        dataset_name = os.path.splitext(os.path.basename(h5ad_path))[0]
        print("\n" + "#" * 100)
        print(f"CACHING NON-SODA METHODS | dataset={dataset_name}")
        print("#" * 100)
        adata_original, non_soda_adata_dict, metric_params = prepare_non_soda_cache(
            h5ad_path=h5ad_path,
            neighborhood_mode=neighborhood_mode,
            base_soda_params=base_soda_params,
            batch_key=batch_key,
            batch_label_map=batch_label_map,
            recompute_original_pca=recompute_original_pca,
            label_key=label_key,
            lisi_n_neighbors=lisi_n_neighbors,
        )

        for value in sweep_values:
            soda_params = dict(base_soda_params)
            soda_params[sweep_param] = value
            print("\n" + "#" * 100)
            print(
                f"SENSITIVITY RUN | dataset={dataset_name} | "
                f"{sweep_param}={value}"
            )
            print("#" * 100)
            soda_adata = run_soda_only(
                h5ad_path=h5ad_path,
                neighborhood_mode=neighborhood_mode,
                soda_params=soda_params,
                batch_key=batch_key,
                batch_label_map=batch_label_map,
                recompute_original_pca=recompute_original_pca,
                label_key=label_key,
                lisi_n_neighbors=lisi_n_neighbors,
            )
            if soda_adata is None:
                continue

            combined_adata_dict = dict(non_soda_adata_dict)
            combined_adata_dict["SODA"] = soda_adata
            metrics = calculate_metrics_for_adata_dict(
                adata_original=adata_original,
                adata_dict=combined_adata_dict,
                metric_params=metric_params,
            )
            ranking = compute_overall_ranking(combined_adata_dict, metrics)
            soda_overall_score, soda_rank = extract_method_overall_score(
                ranking, method_name="SODA"
            )
            soda_metric_scores = extract_method_metric_scores(
                metrics, method_name="SODA"
            )

            row = {
                "dataset": dataset_name,
                "param": sweep_param,
                "value": value,
                "soda_overall_score": soda_overall_score,
                "soda_rank": soda_rank,
                "top_method": ranking[0][0] if ranking else "",
            }
            for metric_name, score in soda_metric_scores.items():
                row[f"metric_{metric_name}"] = score
            records.append(row)
            del combined_adata_dict, metrics, ranking, soda_metric_scores, soda_adata
            release_memory()

        del adata_original, non_soda_adata_dict, metric_params
        release_memory()

    long_csv = os.path.join(sensitivity_root, f"table_{sweep_param}_sensitivity_long.csv")
    with open(long_csv, "w", newline="", encoding="utf-8") as f:
        if records:
            fieldnames = list(records[0].keys())
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)

    # 生成适合论文补充材料的宽表：行=数据集，列=参数取值，值=SODA综合分数
    value_cols = [str(v) for v in sweep_values]
    wide_rows = []
    datasets = sorted({r["dataset"] for r in records})
    for ds in datasets:
        row = {"dataset": ds}
        for v in sweep_values:
            row[str(v)] = ""
        for r in records:
            if r["dataset"] == ds:
                row[str(r["value"])] = r["soda_overall_score"]
        wide_rows.append(row)

    wide_csv = os.path.join(sensitivity_root, f"table_{sweep_param}_sensitivity_wide.csv")
    with open(wide_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["dataset"] + value_cols)
        writer.writeheader()
        writer.writerows(wide_rows)

    # 跨数据集平均分表：可直接回应“整体上是否敏感”
    avg_rows = []
    for v in sweep_values:
        values = [
            r["soda_overall_score"]
            for r in records
            if r["value"] == v and r["soda_overall_score"] is not None
        ]
        avg_rows.append(
            {
                "param": sweep_param,
                "value": v,
                "mean_soda_overall_score": (sum(values) / len(values)) if values else None,
                "n_datasets": len(values),
            }
        )
    avg_rows.sort(
        key=lambda x: x["mean_soda_overall_score"]
        if x["mean_soda_overall_score"] is not None else float("-inf"),
        reverse=True,
    )

    avg_csv = os.path.join(sensitivity_root, f"table_{sweep_param}_sensitivity_mean.csv")
    with open(avg_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["param", "value", "mean_soda_overall_score", "n_datasets"]
        )
        writer.writeheader()
        writer.writerows(avg_rows)

    if generate_sensitivity_plots:
        generate_sensitivity_plots_figures(
            records=records,
            avg_rows=avg_rows,
            sweep_param=sweep_param,
            output_dir=sensitivity_root,
        )

    print("\nSensitivity tables saved:")
    print(f"  - {long_csv}")
    print(f"  - {wide_csv}")
    print(f"  - {avg_csv}")
    return {
        "output_dir": sensitivity_root,
        "long_csv": long_csv,
        "wide_csv": wide_csv,
        "avg_csv": avg_csv,
    }


def generate_sensitivity_plots_figures(records, avg_rows, sweep_param, output_dir):
    """根据敏感性扫描结果输出论文可用图。"""
    if not records:
        return

    # 1) 跨数据集平均分曲线图
    mean_points = [
        (row["value"], row["mean_soda_overall_score"])
        for row in avg_rows
        if row["mean_soda_overall_score"] is not None
    ]
    if mean_points:
        mean_points = sorted(mean_points, key=lambda x: x[0])
        x_vals = [p[0] for p in mean_points]
        y_vals = [p[1] for p in mean_points]
        plt.figure(figsize=(7, 5))
        plt.plot(x_vals, y_vals, marker="o", linewidth=2)
        plt.xlabel(sweep_param)
        plt.ylabel("Mean SODA overall score")
        plt.title(f"SODA sensitivity ({sweep_param}) - cross-dataset mean")
        plt.grid(alpha=0.3)
        plt.tight_layout()
        mean_fig = os.path.join(
            output_dir, f"figure_{sweep_param}_sensitivity_mean.png"
        )
        plt.savefig(mean_fig, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  - {mean_fig}")

    # 2) 数据集分层曲线图
    dataset_to_points = {}
    for row in records:
        score = row.get("soda_overall_score")
        if score is None:
            continue
        dataset_to_points.setdefault(row["dataset"], []).append(
            (row["value"], score)
        )

    if dataset_to_points:
        plt.figure(figsize=(8, 5))
        for dataset, points in sorted(dataset_to_points.items()):
            points = sorted(points, key=lambda x: x[0])
            x_vals = [p[0] for p in points]
            y_vals = [p[1] for p in points]
            plt.plot(x_vals, y_vals, marker="o", linewidth=2, label=dataset)
        plt.xlabel(sweep_param)
        plt.ylabel("SODA overall score")
        plt.title(f"SODA sensitivity ({sweep_param}) - by dataset")
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()
        by_dataset_fig = os.path.join(
            output_dir, f"figure_{sweep_param}_sensitivity_by_dataset.png"
        )
        plt.savefig(by_dataset_fig, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  - {by_dataset_fig}")


def main():
    parser = argparse.ArgumentParser(
        description="Run benchmark on all h5ad files under data/."
    )
    parser.add_argument(
        "--neighborhood-mode",
        choices=["expression", "representation"],
        default="representation",
        help=(
            "Neighborhood metric mode: "
            "'expression' compares adata.X; "
            "'representation' compares method embeddings in obsm."
        ),
    )
    parser.add_argument(
        "--data-files",
        nargs="*",
        default=None,
        help=(
            "Optional dataset file names under data/. "
            "Example: --data-files pbmc4k.h5ad pancreas.h5ad"
        ),
    )
    parser.add_argument(
        "--soda-tau",
        type=float,
        default=0.01,
        help="Default SODA tau (divergence threshold). Lower = more correction.",
    )
    parser.add_argument(
        "--soda-alpha",
        type=float,
        default=0.7,
        help="Default SODA soft-update alpha.",
    )
    parser.add_argument(
        "--soda-epsilon",
        type=float,
        default=0.05,
        help="Default SODA entropy regularization epsilon.",
    )
    parser.add_argument(
        "--sensitivity",
        action="store_true",
        help="Enable SODA hyperparameter sensitivity sweep mode.",
    )
    parser.add_argument(
        "--sensitivity-param",
        choices=["tau", "alpha", "epsilon"],
        default="tau",
        help="Which SODA parameter to sweep when --sensitivity is enabled.",
    )
    parser.add_argument(
        "--tau-values",
        default="0.01,0.02,0.05,0.1,0.2",
        help="Comma-separated tau values for sensitivity sweep.",
    )
    parser.add_argument(
        "--alpha-values",
        default="0.3,0.5,0.7,0.9",
        help="Comma-separated alpha values for sensitivity sweep.",
    )
    parser.add_argument(
        "--epsilon-values",
        default="0.01,0.02,0.05,0.1,0.2",
        help="Comma-separated epsilon values for sensitivity sweep.",
    )
    parser.add_argument(
        "--output-root",
        default="solution",
        help="Root directory for outputs (default: solution).",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Disable benchmark plots for each dataset.",
    )
    parser.add_argument(
        "--sensitivity-plots",
        action="store_true",
        help="Generate sensitivity curve figures in sensitivity mode.",
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
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Whether to recompute uncorrected PCA/UMAP for original baseline "
            "(default: True; pass --no-recompute-original-pca to disable)."
        ),
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
    batch_label_map = parse_label_map(args.batch_label_map)

    data_dir = "data"
    if args.data_files:
        h5_files = [os.path.join(data_dir, f) for f in args.data_files]
    else:
        h5_files = [
            os.path.join(data_dir, f)
            for f in os.listdir(data_dir)
            if f.endswith(".h5ad")
        ]

    if not h5_files:
        print("No H5AD files found in 'data/' directory.")
        return

    base_soda_params = {
        "tau": args.soda_tau,
        "alpha": args.soda_alpha,
        "epsilon": args.soda_epsilon,
    }

    if args.sensitivity:
        value_map = {
            "tau": parse_float_list(args.tau_values),
            "alpha": parse_float_list(args.alpha_values),
            "epsilon": parse_float_list(args.epsilon_values),
        }
        run_sensitivity_sweep(
            h5_files=h5_files,
            output_root=args.output_root,
            neighborhood_mode=args.neighborhood_mode,
            sweep_param=args.sensitivity_param,
            sweep_values=value_map[args.sensitivity_param],
            base_soda_params=base_soda_params,
            generate_sensitivity_plots=args.sensitivity_plots,
            batch_key=args.batch_key,
            batch_label_map=batch_label_map,
            recompute_original_pca=args.recompute_original_pca,
            label_key=args.label_key,
            lisi_n_neighbors=args.lisi_n_neighbors,
        )
        return

    batch_results = []
    for h5ad_path in h5_files:
        result = run_for_dataset(
            h5ad_path,
            output_root=args.output_root,
            generate_plots=not args.no_plots,
            neighborhood_mode=args.neighborhood_mode,
            soda_params=base_soda_params,
            batch_key=args.batch_key,
            batch_label_map=batch_label_map,
            recompute_original_pca=args.recompute_original_pca,
            label_key=args.label_key,
            lisi_n_neighbors=args.lisi_n_neighbors,
        )
        batch_results.append(result)
        release_memory()

    write_batch_ranking_summary(
        batch_results=batch_results,
        output_root=args.output_root,
        neighborhood_mode=args.neighborhood_mode,
    )
    write_batch_gain_vs_original_summary(
        batch_results=batch_results,
        output_root=args.output_root,
    )


if __name__ == "__main__":
    main()
