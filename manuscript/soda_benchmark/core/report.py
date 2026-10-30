# -*- coding: utf-8 -*-
"""
报告生成模块。

通过 MetricResult 标准化结构，可自动处理任意数量的已注册指标，
生成文本报告和综合排名，无需为每个新指标手动修改代码。
"""

import os

from soda_benchmark.metrics import get_registered_metrics


class ReportGenerator:
    """负责汇总指标并输出文本报告。"""

    def __init__(self, adata_original, adata_dict, metrics, output_dir='results'):
        """
        Parameters
        ----------
        adata_original : AnnData
            原始数据
        adata_dict : dict
            方法名 → AnnData
        metrics : dict
            指标键名 → MetricResult
        output_dir : str
            输出目录
        """
        self.adata_original = adata_original
        self.adata_dict = adata_dict
        self.metrics = metrics
        self.output_dir = output_dir

    def generate_report(self, save_path=None):
        if save_path is None:
            save_path = os.path.join(self.output_dir, 'benchmark_report.txt')
        os.makedirs(os.path.dirname(save_path), exist_ok=True)

        # 获取已注册指标的类映射，用于调用 grade()
        registered = get_registered_metrics()

        with open(save_path, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write(" " * 15 + "SODA: SINGLE-CELL BATCH CORRECTION BENCHMARK\n")
            f.write("=" * 80 + "\n\n")

            # ---------------------------------------------------------- #
            # 1. 数据集信息
            # ---------------------------------------------------------- #
            f.write("1. DATASET INFORMATION\n" + "-" * 80 + "\n")
            f.write(f"Total cells: {self.adata_original.n_obs:,}\n")
            f.write(f"Total genes (HVG): {self.adata_original.n_vars:,}\n")

            batch_col = self.adata_original.obs['batch'].astype(str)
            batch_counts = batch_col.value_counts()
            f.write(f"Number of batches: {len(batch_counts)}\n")
            for batch_name, n_cells in batch_counts.items():
                f.write(f"Batch {batch_name} cells: {int(n_cells):,}\n")
            f.write(f"Original clusters: {len(self.adata_original.obs['leiden'].unique())}\n\n")

            # ---------------------------------------------------------- #
            # 2. 评估方法
            # ---------------------------------------------------------- #
            f.write("2. EVALUATED METHODS\n" + "-" * 80 + "\n")
            for method in self.adata_dict.keys():
                n_clusters = len(self.adata_dict[method].obs['leiden'].unique())
                tag = " (baseline)" if method == "Original" else ""
                f.write(f"  ✓ {method:20s}: {n_clusters} clusters{tag}\n")
            f.write(f"\nTotal methods (including baseline): {len(self.adata_dict)}\n\n")

            # ---------------------------------------------------------- #
            # 3. 各指标详细结果（自动遍历所有已注册指标）
            # ---------------------------------------------------------- #
            if self.metrics:
                f.write("3. PERFORMANCE METRICS\n" + "-" * 80 + "\n\n")

                for section_idx, (key, mr) in enumerate(self.metrics.items(), 1):
                    direction = "Higher = Better" if mr.higher_is_better else "Lower = Better"
                    f.write(f"3.{section_idx} {mr.name} ({direction})\n")

                    # 排序：越好在前面
                    sorted_scores = mr.scores.sort_values(
                        'score', ascending=not mr.higher_is_better)

                    # 确定数值格式：分数型 (0‑1) 用 .4f，其他用 .2f
                    is_fractional = (mr.scores['score'].abs().max() <= 2.0)
                    fmt = ".4f" if is_fractional else ".2f"

                    # 获取对应的 grade 函数
                    metric_cls = registered.get(key)
                    grader = metric_cls() if metric_cls else None

                    for rank, (method, row) in enumerate(sorted_scores.iterrows(), 1):
                        value = row['score']
                        grade = grader.grade(value) if grader else ""
                        f.write(
                            f"    {rank}. {method:20s}: {value:8{fmt}}  "
                            f"[95% CI {row['ci_lower']:{fmt}}, {row['ci_upper']:{fmt}}]"
                            f"  [{grade}]\n"
                        )
                    best = sorted_scores.index[0]
                    f.write(f"\n    ★ Best: {best}\n\n")

                # ------------------------------------------------------ #
                # 4. 方法排名（Original 仅作对照，不参与方法榜）
                # ------------------------------------------------------ #
                def _rank_methods(metric_keys):
                    scores = {}
                    for method in self.adata_dict.keys():
                        if method == "Original":
                            continue
                        total_score, n_counted = 0.0, 0
                        for key in metric_keys:
                            if key not in self.metrics:
                                continue
                            mr = self.metrics[key]
                            method_scores = mr.scores['score']
                            score_range = method_scores.max() - method_scores.min()
                            if score_range > 0 and method in method_scores.index:
                                normalized = ((method_scores[method] - method_scores.min())
                                              / score_range)
                                if not mr.higher_is_better:
                                    normalized = 1 - normalized
                                total_score += normalized
                                n_counted += 1
                        if n_counted > 0:
                            scores[method] = total_score / n_counted
                    return sorted(scores.items(), key=lambda x: x[1], reverse=True)

                integration_keys = ["ilisi", "clisi", "asw_batch"]
                batch_keys = ["ilisi", "asw_batch"]
                bio_keys = ["clisi"]

                sorted_methods = _rank_methods(integration_keys)
                sorted_batch = _rank_methods(batch_keys)
                sorted_bio = _rank_methods(bio_keys)

                f.write("4. METHOD RANKING (INTEGRATION SCORE, EXCLUDING ORIGINAL)\n")
                f.write("-" * 80 + "\n")
                for rank, (method, score) in enumerate(sorted_methods, 1):
                    stars = "★" * min(5, int(score * 5) + 1)
                    f.write(f"  {rank}. {method:20s}: {score:.4f}  {stars}\n")
                if sorted_methods:
                    f.write(f"\n  ★★★ TOP PERFORMER: {sorted_methods[0][0]} ★★★\n")

                f.write("\n4.1 Batch-removal ranking (iLISI + ASW_batch)\n")
                for rank, (method, score) in enumerate(sorted_batch, 1):
                    f.write(f"  {rank}. {method:20s}: {score:.4f}\n")

                f.write("\n4.2 Bio-conservation ranking (cLISI)\n")
                for rank, (method, score) in enumerate(sorted_bio, 1):
                    f.write(f"  {rank}. {method:20s}: {score:.4f}\n")
                f.write("\n")

                # ------------------------------------------------------ #
                # 5. 相对 Original baseline 的提升（方向敏感）
                # ------------------------------------------------------ #
                f.write("5. RELATIVE PERFORMANCE VS ORIGINAL BASELINE\n" + "-" * 80 + "\n")
                baseline_name = "Original"
                for key, mr in self.metrics.items():
                    if baseline_name not in mr.scores.index:
                        continue
                    f.write(f"\n5.{key} {mr.name}\n")
                    baseline_score = float(mr.scores.loc[baseline_name, "score"])
                    direction = "higher better" if mr.higher_is_better else "lower better"
                    f.write(
                        f"    Baseline ({baseline_name}) = {baseline_score:.4f} "
                        f"({direction})\n"
                    )

                    sorted_scores = mr.scores.sort_values(
                        'score', ascending=not mr.higher_is_better)
                    for method, row in sorted_scores.iterrows():
                        if method == baseline_name:
                            continue
                        value = float(row["score"])
                        raw_delta = value - baseline_score
                        improvement = raw_delta if mr.higher_is_better else -raw_delta
                        if abs(baseline_score) > 1e-12:
                            improvement_pct = (improvement / abs(baseline_score)) * 100.0
                            pct_text = f"{improvement_pct:+.2f}%"
                        else:
                            pct_text = "n/a"
                        f.write(
                            f"    - {method:20s}: score={value:.4f}, "
                            f"improvement={improvement:+.4f} ({pct_text})\n"
                        )

            f.write("=" * 80 + "\n")
            f.write("SODA: Sinkhorn Optimal Transport with Divergence-Based Adaptation\n")
            f.write("=" * 80 + "\n")

        with open(save_path, 'r', encoding='utf-8') as f:
            print("\n" + f.read())

        print(f"✓ Report saved to: {save_path}")

