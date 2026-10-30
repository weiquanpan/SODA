# -*- coding: utf-8 -*-
"""
SingleCellBatchBenchmark：基准测试主调度器。

通过注册机制自动发现所有方法与指标，无需硬编码方法/指标列表。
"""

import os
import sys
import warnings
import gc

# Windows 控制台默认使用 GBK 编码，无法输出 ✓✗★ 等字符，强制切换为 UTF-8
if sys.stdout and hasattr(sys.stdout, 'reconfigure'):
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass
if sys.stderr and hasattr(sys.stderr, 'reconfigure'):
    try:
        sys.stderr.reconfigure(encoding='utf-8')
    except Exception:
        pass

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
from scipy import sparse

from soda_benchmark.methods import get_registered_methods
from soda_benchmark.metrics import get_registered_metrics
from soda_benchmark.plotting import PlotGenerator
from soda_benchmark.core.report import ReportGenerator

warnings.filterwarnings('ignore')

# 全局绘图设置
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
sc.settings.verbosity = 0


class SingleCellBatchBenchmark:
    """
    Comprehensive benchmarking framework for single-cell batch correction.

    调度批次校正方法、评估指标与可视化模块。
    通过注册表自动发现方法与指标，新增方法/指标无需修改本类。
    """

    def __init__(self, h5ad_path="data/pbmc4k-0.h5ad",
                 method_params=None, metric_params=None,
                 selected_methods=None, batch_key='batch',
                 batch_label_map=None, recompute_original_pca=False):
        """
        Parameters
        ----------
        h5ad_path : str
            预处理数据文件路径
        method_params : dict, optional
            方法级参数，键为方法名，值为该方法的 kwargs。
            例如 {'SODA': {'n_clusters': 50}, 'scVI': {'max_epochs': 100}}
        metric_params : dict, optional
            所有指标共享的参数（n_bootstrap, ci_level, random_state）
        selected_methods : list[str], optional
            若提供，则仅运行该列表中的方法（例如 ['SODA']）。
        """
        self.h5ad_path = h5ad_path
        self.adata = None
        self.adata_dict = {}
        self.metrics = {}
        self.selected_methods = selected_methods
        self.batch_key = batch_key
        self.batch_label_map = {
            str(k): str(v) for k, v in (batch_label_map or {}).items()
        }
        self.recompute_original_pca = bool(recompute_original_pca)

        self.method_params = method_params or {
            'SODA': {
                'n_clusters': 50,
                'max_iter': 10,
                'epsilon': 0.05,
                'divergence_threshold': 0.05,
            },
            'scVI': {
                'max_epochs': 100,
            },
            'CellANOVA': {
                'integrate_key': 'batch',
                'control_batches': None,
            },
            'SMAI': {
                'integrate_key': 'batch',
                'nfeatures': 2000,
                't_iter': 3,
                'knn': 30,
                'r_max': 200,
                'prop_align': 0.3,
                'cutoff': 1.001,
                'cutoff_t': 1.5,
            },
        }
        self.metric_params = metric_params or {
            'n_bootstrap': 1000,
            'ci_level': 95,
            'random_state': 42,
        }

        self.available_methods = {}
        self._check_dependencies()

    # ------------------------------------------------------------------ #
    # 依赖检查
    # ------------------------------------------------------------------ #

    def _check_dependencies(self):
        """检查所有已注册方法的依赖是否可用。"""
        registered = get_registered_methods()
        for name, method_cls in registered.items():
            available, msg = method_cls.check_dependencies()
            if available:
                self.available_methods[name] = method_cls
            else:
                print(f"Warning: {msg}")

        if self.selected_methods is not None:
            filtered = {}
            for name in self.selected_methods:
                if name in self.available_methods:
                    filtered[name] = self.available_methods[name]
                else:
                    print(f"Warning: selected method '{name}' is unavailable and will be skipped.")
            self.available_methods = filtered

    # ------------------------------------------------------------------ #
    # 数据加载
    # ------------------------------------------------------------------ #

    def load_preprocessed_data(self):
        """从 H5AD 文件加载预处理好的单细胞数据。"""
        print(f"Loading data from {self.h5ad_path}...")

        if not os.path.exists(self.h5ad_path):
            raise FileNotFoundError(
                f"Data file not found: {self.h5ad_path}\n"
                f"Please ensure the preprocessed data is available."
            )

        adata = sc.read_h5ad(self.h5ad_path)

        if self.batch_key not in adata.obs.columns:
            raise KeyError(
                f"Batch key '{self.batch_key}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )

        # 统一将选定批次列映射到标准 batch 字段，保证下游方法兼容。
        adata.obs['batch'] = adata.obs[self.batch_key].astype(str)

        print("Data loaded successfully!")
        print(f"  - Cells: {adata.n_obs:,}")
        print(f"  - Genes (HVG): {adata.n_vars:,}")
        print(f"  - Batch key: {self.batch_key}")
        batch_counts = adata.obs['batch'].value_counts()
        print(f"  - Number of batches: {len(batch_counts)}")
        for raw_batch, count in batch_counts.items():
            label = self.batch_label_map.get(str(raw_batch), str(raw_batch))
            print(f"    * {label} ({raw_batch}): {int(count):,} cells")
        print(f"  - Original clusters: {adata.obs['leiden'].nunique()}")

        need_recompute_umap = self.recompute_original_pca or ('X_umap' not in adata.obsm)
        if need_recompute_umap:
            print("  - Computing UMAP for original baseline...")
            if self.recompute_original_pca or ('X_pca' not in adata.obsm):
                print("  - Recomputing uncorrected PCA from current expression matrix...")
                sc.tl.pca(adata, n_comps=50, random_state=42)
            sc.pp.neighbors(adata, n_pcs=30, random_state=42)
            sc.tl.leiden(adata, resolution=0.5, random_state=42)
            sc.tl.umap(adata, random_state=42)

        self.adata = adata.copy()
        return self.adata

    # ------------------------------------------------------------------ #
    # 方法执行
    # ------------------------------------------------------------------ #

    def _slim_method_output(self, adata):
        """在 representation 模式下压缩方法输出，降低驻留内存。"""
        mode = str(self.metric_params.get("neighborhood_mode", "expression")).lower()
        if mode != "representation":
            return adata

        try:
            adata.layers.clear()
        except Exception:
            pass
        try:
            adata.obsp.clear()
        except Exception:
            pass
        try:
            adata.varm.clear()
        except Exception:
            pass
        try:
            adata.varp.clear()
        except Exception:
            pass

        # representation 评估不依赖表达矩阵，置为零稀疏矩阵可显著节省内存
        try:
            adata.X = sparse.csr_matrix((adata.n_obs, adata.n_vars), dtype=np.float32)
        except Exception:
            pass
        return adata

    def _apply_all_methods(self):
        """依次执行所有可用的批次校正方法。"""
        print("\n" + "=" * 80)
        print("APPLYING BATCH CORRECTION METHODS")
        print("=" * 80)

        self.adata_dict = {'Original': self.adata.copy()}

        for name, method_cls in self.available_methods.items():
            params = self.method_params.get(name, {})
            method_instance = method_cls(**params)
            try:
                adata_input = self.adata.copy()
            except MemoryError as e:
                print(f"  ✗ {name} skipped: insufficient memory while copying input ({e})")
                gc.collect()
                continue

            try:
                result = method_instance.apply(adata_input)
            except MemoryError as e:
                print(f"  ✗ {name} failed due to insufficient memory: {e}")
                result = None
            finally:
                del adata_input
                gc.collect()

            if result is not None:
                result = self._slim_method_output(result)
                self.adata_dict[name] = result

        print(f"\n✓ Applied {len(self.adata_dict) - 1} batch correction methods")
        return self.adata_dict

    # ------------------------------------------------------------------ #
    # 指标计算
    # ------------------------------------------------------------------ #

    def _calculate_all_metrics(self):
        """依次计算所有已注册的评估指标。"""
        print("\n" + "=" * 80)
        print("CALCULATING EVALUATION METRICS")
        print("=" * 80)

        registered = get_registered_metrics()
        for key, metric_cls in registered.items():
            metric_instance = metric_cls(**self.metric_params)
            result = metric_instance.calculate(
                adata_original=self.adata_dict['Original'],
                adata_dict=self.adata_dict,
            )
            self.metrics[key] = result

        print("\n✓ All metrics calculated")
        return self.metrics

    # ------------------------------------------------------------------ #
    # 完整流程
    # ------------------------------------------------------------------ #

    def run_full_benchmark(self, output_dir='results', generate_plots=False):
        """运行完整的基准测试流程。"""
        os.makedirs(output_dir, exist_ok=True)

        print("\n" + "=" * 80)
        print("SODA: SINGLE-CELL BATCH CORRECTION BENCHMARK")
        print("=" * 80)

        self.load_preprocessed_data()
        self._apply_all_methods()
        self._calculate_all_metrics()

        if generate_plots:
            plotter = PlotGenerator(
                self.adata_dict,
                self.metrics,
                output_dir=output_dir,
                batch_key='batch',
                batch_label_map=self.batch_label_map,
            )
            try:
                plotter.generate_all_plots()
            except MemoryError as e:
                print(f"Warning: plotting skipped due to memory error: {e}")
            except Exception as e:
                print(f"Warning: plotting failed and will be skipped: {e}")

        reporter = ReportGenerator(
            adata_original=self.adata,
            adata_dict=self.adata_dict,
            metrics=self.metrics,
            output_dir=output_dir,
        )
        reporter.generate_report()

        print("\n" + "=" * 80)
        print(f"BENCHMARK COMPLETE! Results saved to {output_dir}/")
        print("=" * 80)

