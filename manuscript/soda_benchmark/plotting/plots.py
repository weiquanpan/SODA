# -*- coding: utf-8 -*-
"""
绘图模块：UMAP 对比、聚类分布与指标概览。

现已适配 MetricResult 标准结构，可自动处理任意数量的已注册指标。
"""
import os
import gc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


class PlotGenerator:
    """统一生成全部可视化图表。"""

    def __init__(self, adata_dict, metrics, output_dir='results',
                 batch_key='batch', batch_label_map=None):
        """
        Parameters
        ----------
        adata_dict : dict
            方法名 → AnnData
        metrics : dict
            指标键名 → MetricResult
        output_dir : str
            输出目录
        """
        self.adata_dict = adata_dict
        self.metrics = metrics
        self.output_dir = output_dir
        self.batch_key = batch_key
        self.batch_label_map = {
            str(k): str(v) for k, v in (batch_label_map or {}).items()
        }
        os.makedirs(self.output_dir, exist_ok=True)
        self.max_methods_per_figure = 6
        self.max_metrics_per_figure = 3

    def _safe_savefig(self, fig, save_path):
        """稳健保存图片：固定 300dpi，不降分辨率。"""
        attempts = [
            ("tight",),
            (None,),
        ]
        for (bbox,) in attempts:
            try:
                kwargs = {"dpi": 300}
                if bbox is not None:
                    kwargs["bbox_inches"] = bbox
                fig.savefig(save_path, **kwargs)
                return True
            except MemoryError:
                gc.collect()
            except Exception as e:
                msg = str(e).lower()
                if "memory" in msg or "bad allocation" in msg:
                    gc.collect()
                    continue
                raise
        print(f"  ✗ Skip saving due to insufficient memory: {save_path}")
        return False

    def _get_batch_palette(self, n_batches):
        """
        获取批次配色（优先使用高对比、色盲友好颜色）。
        """
        if n_batches <= 0:
            return sns.color_palette("tab10", 1)
        if n_batches == 1:
            return ["#0072B2"]
        if n_batches == 2:
            # Okabe-Ito 高对比配色：蓝 + 橙
            return ["#0072B2", "#E69F00"]
        return sns.color_palette("tab20", n_batches)

    # ================================================================== #
    # UMAP 可视化
    # ================================================================== #

    def plot_umap_comparison(self, save_path=None):
        if save_path is None:
            save_path = os.path.join(self.output_dir, 'umap_comparison.png')

        valid_methods = {
            name: adata for name, adata in self.adata_dict.items()
            if 'X_umap' in adata.obsm
        }

        if not valid_methods:
            print("  ✗ No methods with UMAP available for plotting")
            return

        method_items = list(valid_methods.items())
        base_path, ext = os.path.splitext(save_path)
        n_chunks = max(
            1, (len(method_items) + self.max_methods_per_figure - 1) // self.max_methods_per_figure
        )

        for chunk_idx in range(n_chunks):
            start = chunk_idx * self.max_methods_per_figure
            end = min(len(method_items), (chunk_idx + 1) * self.max_methods_per_figure)
            chunk_items = method_items[start:end]
            n_methods = len(chunk_items)
            ncols, nrows = 3, (n_methods + 2) // 3

            fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
            axes = axes.flatten() if n_methods > 1 else [axes]

            for idx, (method_name, adata) in enumerate(chunk_items):
                ax = axes[idx]
                if self.batch_key not in adata.obs.columns:
                    raise KeyError(
                        f"Batch key '{self.batch_key}' not found in adata.obs for method '{method_name}'."
                    )

                batch_col = adata.obs[self.batch_key].astype(str)
                batch_categories = sorted(batch_col.unique())
                color_palette = self._get_batch_palette(len(batch_categories))
                marker_cycle = ['o', '^', 's', 'D', 'P', 'X', 'v', '<', '>']

                for batch_idx, batch in enumerate(batch_categories):
                    mask = batch_col == batch
                    label = self.batch_label_map.get(str(batch), str(batch))
                    ax.scatter(
                        adata.obsm['X_umap'][mask, 0],
                        adata.obsm['X_umap'][mask, 1],
                        color=color_palette[batch_idx % len(color_palette)],
                        label=label,
                        marker=marker_cycle[batch_idx % len(marker_cycle)],
                        s=10,
                        alpha=0.75,
                        edgecolors='white',
                        linewidths=0.15
                    )

                ax.set_title(f'{method_name}', fontsize=14, fontweight='bold', pad=10)
                ax.set_xlabel('UMAP 1', fontsize=11)
                ax.set_ylabel('UMAP 2', fontsize=11)
                ax.legend(markerscale=2, loc='upper right', fontsize=9, title=self.batch_key)
                ax.grid(True, alpha=0.2)
                ax.set_aspect('equal', 'box')

            for idx in range(n_methods, len(axes)):
                axes[idx].axis('off')

            plt.suptitle('UMAP Visualization - Batch Integration', fontsize=16, fontweight='bold', y=1.00)
            plt.tight_layout()
            chunk_path = f"{base_path}_part{chunk_idx + 1}{ext}" if n_chunks > 1 else save_path
            saved = self._safe_savefig(fig, chunk_path)
            plt.close()
            if saved:
                print(f"  ✓ Saved: {chunk_path}")

    def plot_umap_clusters(self, save_path=None):
        if save_path is None:
            save_path = os.path.join(self.output_dir, 'umap_clusters.png')
        valid_methods = {
            name: adata for name, adata in self.adata_dict.items()
            if 'X_umap' in adata.obsm
        }

        if not valid_methods:
            print("  ✗ No methods with UMAP available for plotting")
            return

        method_items = list(valid_methods.items())
        base_path, ext = os.path.splitext(save_path)
        n_chunks = max(
            1, (len(method_items) + self.max_methods_per_figure - 1) // self.max_methods_per_figure
        )

        for chunk_idx in range(n_chunks):
            start = chunk_idx * self.max_methods_per_figure
            end = min(len(method_items), (chunk_idx + 1) * self.max_methods_per_figure)
            chunk_items = method_items[start:end]
            n_methods = len(chunk_items)
            ncols, nrows = 3, (n_methods + 2) // 3

            fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
            axes = axes.flatten() if n_methods > 1 else [axes]

            for idx, (method_name, adata) in enumerate(chunk_items):
                ax = axes[idx]
                ax.scatter(
                    adata.obsm['X_umap'][:, 0],
                    adata.obsm['X_umap'][:, 1],
                    c=adata.obs['leiden'].astype(int),
                    cmap='tab20',
                    s=8,
                    alpha=0.6,
                    edgecolors='none'
                )

                n_clusters = len(adata.obs['leiden'].unique())
                ax.set_title(f'{method_name}\n({n_clusters} clusters)', fontsize=12, fontweight='bold', pad=10)
                ax.set_xlabel('UMAP 1', fontsize=11)
                ax.set_ylabel('UMAP 2', fontsize=11)
                ax.grid(True, alpha=0.2)
                ax.set_aspect('equal', 'box')

            for idx in range(n_methods, len(axes)):
                axes[idx].axis('off')

            plt.suptitle('UMAP Visualization - Biological Clusters', fontsize=16, fontweight='bold', y=1.00)
            plt.tight_layout()
            chunk_path = f"{base_path}_part{chunk_idx + 1}{ext}" if n_chunks > 1 else save_path
            saved = self._safe_savefig(fig, chunk_path)
            plt.close()
            if saved:
                print(f"  ✓ Saved: {chunk_path}")

    # ================================================================== #
    # 指标概览（通用化：自动适配任意数量的 MetricResult）
    # ================================================================== #

    def plot_metrics_summary(self, save_path=None):
        if save_path is None:
            save_path = os.path.join(self.output_dir, 'metrics_summary.png')
        if not self.metrics:
            return

        metric_keys = list(self.metrics.keys())
        base_path, ext = os.path.splitext(save_path)
        n_chunks = max(
            1, (len(metric_keys) + self.max_metrics_per_figure - 1) // self.max_metrics_per_figure
        )
        for chunk_idx in range(n_chunks):
            start = chunk_idx * self.max_metrics_per_figure
            end = min(len(metric_keys), (chunk_idx + 1) * self.max_metrics_per_figure)
            chunk_keys = metric_keys[start:end]
            n_metrics = len(chunk_keys)
            fig, axes = plt.subplots(1, max(n_metrics, 1), figsize=(7 * n_metrics, 5))
            if n_metrics == 1:
                axes = [axes]

            for ax_idx, key in enumerate(chunk_keys):
                mr = self.metrics[key]
                ax = axes[ax_idx]
                # 如果有 per_cell 分布数据，使用 boxplot
                if 'per_cell' in mr.extra:
                    per_cell = mr.extra['per_cell']
                    bp = ax.boxplot(
                        [per_cell[col].values for col in per_cell.columns],
                        labels=per_cell.columns,
                        patch_artist=True,
                        showfliers=True,
                        widths=0.6,
                    )
                    colors = sns.color_palette("Set2", len(per_cell.columns))
                    for patch, color in zip(bp['boxes'], colors):
                        patch.set_facecolor(color)
                        patch.set_alpha(0.7)
                    ax.set_ylabel('Median Rank Change', fontsize=11)
                    ax.grid(True, alpha=0.3, axis='y')
                    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

                    for i, method in enumerate(mr.scores.index, start=1):
                        score_val = mr.scores.loc[method, 'score']
                        ci_low = mr.scores.loc[method, 'ci_lower']
                        ci_up = mr.scores.loc[method, 'ci_upper']
                        ax.text(
                            i, score_val,
                            f'{score_val:.1f}\n[{ci_low:.1f}, {ci_up:.1f}]',
                            ha='center', va='bottom',
                            fontsize=9, fontweight='bold',
                        )
                else:
                    # 通用 bar chart
                    data = mr.scores
                    colors = sns.color_palette("Set2", len(data))
                    bars = ax.bar(range(len(data)), data['score'].values,
                                  color=colors, alpha=0.8)
                    ax.set_xticks(range(len(data)))
                    ax.set_xticklabels(data.index, rotation=45, ha='right')
                    ax.set_ylabel(f'{mr.name} Score', fontsize=11)
                    ax.grid(True, alpha=0.3, axis='y')

                    score_vals = data['score'].values
                    ci_lower_vals = data['ci_lower'].values
                    ci_upper_vals = data['ci_upper'].values
                    lower_bound = np.minimum(np.minimum(ci_lower_vals, ci_upper_vals), score_vals)
                    upper_bound = np.maximum(np.maximum(ci_lower_vals, ci_upper_vals), score_vals)
                    yerr = [
                        score_vals - lower_bound,
                        upper_bound - score_vals,
                    ]
                    ax.errorbar(range(len(data)), data['score'].values,
                                yerr=yerr, fmt='none', ecolor='black', capsize=4)

                    for bar, (_, row) in zip(bars, data.iterrows()):
                        value = row['score']
                        ax.text(bar.get_x() + bar.get_width() / 2, value,
                                f'{value:.4f}', ha='center', va='bottom', fontsize=9)

                    if mr.higher_is_better:
                        ax.set_ylim([0, 1.1])
                        ax.axhline(y=0.8, color='red', linestyle='--', alpha=0.5,
                                   label='Excellent (>0.8)')
                        ax.legend(fontsize=9)

                direction = "Higher = Better" if mr.higher_is_better else "Lower = Better"
                ax.set_title(f'{mr.name}\n({direction})',
                             fontsize=13, fontweight='bold', pad=10)

            plt.suptitle('SODA Benchmark: Performance Metrics',
                         fontsize=16, fontweight='bold', y=1.00)
            plt.tight_layout()
            chunk_path = f"{base_path}_part{chunk_idx + 1}{ext}" if n_chunks > 1 else save_path
            saved = self._safe_savefig(fig, chunk_path)
            plt.close()
            if saved:
                print(f"  ✓ Saved: {chunk_path}")

    # ================================================================== #
    # 邻域保持分布图（仅当有 per_cell 数据时生成）
    # ================================================================== #

    def plot_neighborhood_distribution(self, save_path=None):
        if save_path is None:
            save_path = os.path.join(self.output_dir, 'neighborhood_distribution.png')

        # 收集所有包含 per_cell 的指标
        dist_metrics = {
            k: mr for k, mr in self.metrics.items()
            if 'per_cell' in mr.extra
        }
        if not dist_metrics:
            return

        for metric_key, mr in dist_metrics.items():
            per_cell = mr.extra['per_cell']
            summary = mr.scores
            method_names = list(per_cell.columns)
            if not method_names:
                print(f"  ✗ No methods available for {mr.name} distribution")
                continue

            chunk_size = 6
            base_path, ext = os.path.splitext(save_path)
            n_chunks = max(1, (len(method_names) + chunk_size - 1) // chunk_size)

            for chunk_idx in range(n_chunks):
                start = chunk_idx * chunk_size
                end = min(len(method_names), (chunk_idx + 1) * chunk_size)
                chunk_methods = method_names[start:end]
                width = max(8, len(chunk_methods) * 1.3)
                fig, ax = plt.subplots(figsize=(width, 6))

                parts = ax.violinplot(
                    [per_cell[col].values for col in chunk_methods],
                    positions=range(len(chunk_methods)),
                    showmeans=True,
                    showmedians=True,
                    widths=0.7,
                )

                colors = sns.color_palette("Set2", len(chunk_methods))
                for i, pc in enumerate(parts['bodies']):
                    pc.set_facecolor(colors[i])
                    pc.set_alpha(0.7)

                ax.set_xticks(range(len(chunk_methods)))
                ax.set_xticklabels(chunk_methods, rotation=45, ha='right')
                ax.set_title(f'{mr.name} Distribution',
                             fontsize=14, fontweight='bold', pad=15)
                ax.set_ylabel('Rank Change (log scale)', fontsize=12)
                ax.set_yscale('log')
                ax.grid(True, alpha=0.3, axis='y', which='both')

                summary_subset = summary.loc[chunk_methods]
                medians = summary_subset['score'].values
                ax.plot(range(len(medians)), medians, 'r-', linewidth=2,
                        label='Median', marker='o', markersize=8)
                for idx, method in enumerate(chunk_methods):
                    ci_text = (
                        f"[{summary_subset.loc[method, 'ci_lower']:.1f}, "
                        f"{summary_subset.loc[method, 'ci_upper']:.1f}]"
                    )
                    ax.text(idx, medians[idx], ci_text,
                            fontsize=9, ha='center', va='bottom')
                ax.legend(fontsize=11)

                plt.tight_layout()
                chunk_path = (
                    f"{base_path}_part{chunk_idx + 1}{ext}" if n_chunks > 1 else save_path
                )
                saved = self._safe_savefig(fig, chunk_path)
                plt.close()
                if saved:
                    print(f"  ✓ Saved: {chunk_path}")

    # ================================================================== #
    # 一键生成全部图表
    # ================================================================== #

    def generate_all_plots(self):
        print("\n" + "=" * 80)
        print("GENERATING VISUALIZATIONS")
        print("=" * 80)

        plot_tasks = [
            ("UMAP comparison", self.plot_umap_comparison),
            ("UMAP clusters", self.plot_umap_clusters),
            ("metrics summary", self.plot_metrics_summary),
            ("distribution", self.plot_neighborhood_distribution),
        ]
        for task_name, task_func in plot_tasks:
            try:
                task_func()
            except MemoryError as e:
                print(f"  ✗ Plot '{task_name}' skipped due to memory error: {e}")
                gc.collect()
            except Exception as e:
                print(f"  ✗ Plot '{task_name}' failed: {e}")
                gc.collect()

        print("\n✓ All visualizations generated")

