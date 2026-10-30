# -*- coding: utf-8 -*-
"""
Seurat CCA vs RPCA vs LIGER 批次整合对比分析。

通过 rpy2 在 Python 中调用 R，分别使用：
  - CCA  (Canonical Correlation Analysis, Seurat)
  - RPCA (Reciprocal PCA, Seurat)
  - LIGER (integrative NMF, rliger)
对单细胞 RNA-seq 数据进行批次整合，并全面对比三种方法的效果。

依赖
----
Python : rpy2 (>= 3.4), scanpy, numpy, scipy, matplotlib, scikit-learn
R      : Seurat (>= 5.0.0), rliger (>= 2.0.0)

用法
----
    python run_seurat_comparison.py
    python run_seurat_comparison.py --data data/pbmc4k.h5ad --output results/seurat
"""

import os
import sys
import time
import argparse

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse
from sklearn.metrics import silhouette_score, adjusted_rand_score

# Windows 控制台 UTF-8
if sys.stdout and hasattr(sys.stdout, 'reconfigure'):
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass


# ================================================================
#  1. 环境检查
# ================================================================

def check_r_environment():
    """验证 rpy2、R Seurat 和 rliger 是否可用，打印版本信息。"""
    try:
        import rpy2
        import rpy2.robjects as ro
    except ImportError:
        print("错误: rpy2 未安装。请执行: pip install rpy2")
        sys.exit(1)

    print(f"  rpy2 版本  : {rpy2.__version__}")

    try:
        ro.r('suppressPackageStartupMessages(library(Seurat))')
        version_str = str(ro.r('as.character(packageVersion("Seurat"))')[0])
        print(f"  Seurat 版本: {version_str}")
    except Exception as e:
        print(f"错误: 无法加载 R 包 Seurat — {e}")
        print("请在 R 中执行: install.packages('Seurat')")
        sys.exit(1)

    try:
        ro.r('suppressPackageStartupMessages(library(rliger))')
        liger_ver = str(ro.r('as.character(packageVersion("rliger"))')[0])
        print(f"  rliger 版本: {liger_ver}")
    except Exception as e:
        print(f"错误: 无法加载 R 包 rliger — {e}")
        print("请在 R 中执行: install.packages('rliger')")
        sys.exit(1)


# ================================================================
#  2. 数据加载与准备
# ================================================================

def load_data(h5ad_path):
    """加载 h5ad 数据并做基础检查。"""
    if not os.path.exists(h5ad_path):
        print(f"错误: 数据文件不存在 — {h5ad_path}")
        sys.exit(1)

    adata = sc.read_h5ad(h5ad_path)

    if 'batch' not in adata.obs.columns:
        print("错误: 数据中缺少 'batch' 列，无法进行批次整合")
        sys.exit(1)

    n_batches = adata.obs['batch'].nunique()
    if n_batches < 2:
        print(f"错误: 数据中只有 {n_batches} 个批次，至少需要 2 个")
        sys.exit(1)

    print(f"\n数据概况:")
    print(f"  细胞数       : {adata.n_obs:,}")
    print(f"  基因数 (HVG) : {adata.n_vars:,}")
    batch_counts = adata.obs['batch'].value_counts()
    for b, cnt in batch_counts.items():
        print(f"  Batch {b}      : {cnt:,} cells")

    return adata


def extract_expression_matrix(adata):
    """
    提取表达矩阵，优先使用原始计数 (raw counts)。

    Returns
    -------
    X : np.ndarray  (n_cells, n_genes)
    gene_names : list[str]
    source : str  描述数据来源
    """
    gene_names = list(adata.var_names)

    # 优先: adata.raw
    if adata.raw is not None:
        try:
            raw_ad = adata.raw.to_adata()
            hvg = [g for g in gene_names if g in raw_ad.var_names]
            if len(hvg) > 0:
                X = raw_ad[:, hvg].X
                if sparse.issparse(X):
                    X = np.asarray(X.todense(), dtype=np.float64)
                else:
                    X = np.asarray(X, dtype=np.float64)
                return X, hvg, "adata.raw (原始计数)"
        except Exception:
            pass

    # 其次: layers['counts']
    if 'counts' in adata.layers:
        X = adata.layers['counts']
        if sparse.issparse(X):
            X = np.asarray(X.todense(), dtype=np.float64)
        else:
            X = np.asarray(X, dtype=np.float64)
        return X, gene_names, "layers['counts']"

    # 兜底: adata.X
    X = adata.X
    if sparse.issparse(X):
        X = np.asarray(X.todense(), dtype=np.float64)
    else:
        X = np.asarray(X, dtype=np.float64)
    return X, gene_names, "adata.X (可能已归一化)"


# ================================================================
#  3. Seurat 整合 (在 R 中执行)
# ================================================================

def _r_matrix_to_numpy(ro, r_varname):
    """从 R 全局环境取出矩阵并转为 numpy ndarray (cells × dims)。"""
    r_obj = ro.globalenv[r_varname]
    arr = np.array(r_obj, dtype=np.float64)
    if arr.ndim == 2:
        return arr
    dims = ro.r('dim')(r_obj)
    nrow, ncol = int(dims[0]), int(dims[1])
    return arr.reshape((nrow, ncol), order='F')


def run_r_integration(adata):
    """
    将数据传入 R，分别运行 CCA、RPCA (Seurat) 和 LIGER (iNMF) 整合。

    CCA 和 RPCA 共享 Seurat 标准预处理；LIGER 使用自己的预处理流水线。

    Returns
    -------
    dict 包含:
        cca_embeddings   : np.ndarray (n_cells, n_dims)
        rpca_embeddings  : np.ndarray (n_cells, n_dims)
        liger_embeddings : np.ndarray (n_cells, n_dims)
        preprocess_time  : float
        cca_time         : float
        rpca_time        : float
        liger_time       : float
    """
    import rpy2.robjects as ro

    # ---- 提取表达矩阵 ----
    X, gene_names, source = extract_expression_matrix(adata)
    print(f"  表达矩阵来源: {source}")

    X_t = np.ascontiguousarray(X.T, dtype=np.float64)   # genes × cells
    n_genes, n_cells = X_t.shape

    cell_names = list(adata.obs_names)
    batch_labels = list(adata.obs['batch'].astype(str))

    print(f"  传输矩阵到 R: {n_genes} genes × {n_cells} cells ...")

    # ---- 传输到 R 全局环境（显式构造 R 对象，无需 converter） ----
    r_vec = ro.FloatVector(X_t.flatten(order='F'))
    ro.globalenv['expr_mat'] = ro.r.matrix(r_vec, nrow=n_genes, ncol=n_cells)
    ro.globalenv['gene_names_vec'] = ro.StrVector(gene_names)
    ro.globalenv['cell_names_vec'] = ro.StrVector(cell_names)
    ro.globalenv['batch_vec'] = ro.StrVector(batch_labels)

    # ---- Seurat 预处理 (CCA / RPCA 共享) ----
    print("\n  [Seurat] 预处理 (Normalize → HVG → Scale → PCA) ...")
    t0 = time.time()

    ro.r('''
        suppressPackageStartupMessages({
            library(Seurat)
            library(rliger)
        })
        options(future.globals.maxSize = 8 * 1024^3)

        rownames(expr_mat) <- gene_names_vec
        colnames(expr_mat) <- cell_names_vec

        obj <- CreateSeuratObject(counts = expr_mat)
        obj$batch <- batch_vec
        obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)

        obj <- NormalizeData(obj, verbose = FALSE)
        obj <- FindVariableFeatures(obj, verbose = FALSE)
        obj <- ScaleData(obj, verbose = FALSE)
        obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
    ''')

    t_preprocess = time.time() - t0
    print(f"  预处理完成 ({t_preprocess:.1f}s)")

    # ---- CCA 整合 ----
    print("\n  [Seurat] CCA 整合 ...")
    t1 = time.time()

    ro.r('''
        obj_cca <- IntegrateLayers(
            object = obj,
            method = CCAIntegration,
            orig.reduction = "pca",
            new.reduction = "integrated.cca",
            verbose = FALSE
        )
        cca_emb <- Embeddings(obj_cca, reduction = "integrated.cca")
    ''')

    t_cca = time.time() - t1
    cca_embeddings = _r_matrix_to_numpy(ro, 'cca_emb')
    print(f"  CCA 完成 ({t_cca:.1f}s)  嵌入维度: {cca_embeddings.shape}")

    # ---- RPCA 整合 ----
    print("\n  [Seurat] RPCA 整合 ...")
    t2 = time.time()

    ro.r('''
        obj_rpca <- IntegrateLayers(
            object = obj,
            method = RPCAIntegration,
            orig.reduction = "pca",
            new.reduction = "integrated.rpca",
            verbose = FALSE
        )
        rpca_emb <- Embeddings(obj_rpca, reduction = "integrated.rpca")
    ''')

    t_rpca = time.time() - t2
    rpca_embeddings = _r_matrix_to_numpy(ro, 'rpca_emb')
    print(f"  RPCA 完成 ({t_rpca:.1f}s)  嵌入维度: {rpca_embeddings.shape}")

    # ---- LIGER (iNMF) 整合 ----
    # LIGER 有自己的预处理管线，需要从原始 Seurat 对象重新构建
    print("\n  [LIGER] iNMF 整合 ...")
    t3 = time.time()

    ro.r('''
        obj_liger <- CreateSeuratObject(counts = expr_mat)
        obj_liger$batch <- batch_vec
        obj_liger[["RNA"]] <- split(obj_liger[["RNA"]], f = obj_liger$batch)

        # LIGER 专用预处理: normalize (无 log1p) → selectGenes → scaleNotCenter
        obj_liger <- normalize(obj_liger)
        obj_liger <- selectGenes(obj_liger)
        obj_liger <- scaleNotCenter(obj_liger)

        # iNMF 整合 + quantile normalization
        obj_liger <- runINMF(obj_liger, k = 20)
        obj_liger <- quantileNorm(obj_liger)

        liger_emb <- Embeddings(obj_liger, reduction = "inmfNorm")
    ''')

    t_liger = time.time() - t3
    liger_embeddings = _r_matrix_to_numpy(ro, 'liger_emb')
    print(f"  LIGER 完成 ({t_liger:.1f}s)  嵌入维度: {liger_embeddings.shape}")

    # 清理 R 内存
    ro.r('rm(list = ls()); invisible(gc())')

    return {
        'cca_embeddings': cca_embeddings,
        'rpca_embeddings': rpca_embeddings,
        'liger_embeddings': liger_embeddings,
        'preprocess_time': t_preprocess,
        'cca_time': t_cca,
        'rpca_time': t_rpca,
        'liger_time': t_liger,
    }


# ================================================================
#  4. Python 后处理
# ================================================================

def build_adata_variants(adata, integration_results):
    """
    基于 R 返回的整合嵌入，构建四份 AnnData:
    Original / Seurat_CCA / Seurat_RPCA / LIGER，各自包含 leiden 聚类和 UMAP。

    Returns
    -------
    dict[str, AnnData]
    """
    adata_dict = {}

    # ---- Original (无整合) ----
    ad_orig = adata.copy()
    if 'X_pca' not in ad_orig.obsm:
        sc.tl.pca(ad_orig, n_comps=30, random_state=42)
    if 'X_umap' not in ad_orig.obsm:
        sc.pp.neighbors(ad_orig, n_pcs=30, random_state=42)
        sc.tl.umap(ad_orig, random_state=42)
    if 'leiden' not in ad_orig.obs:
        sc.pp.neighbors(ad_orig, n_pcs=30, random_state=42)
        sc.tl.leiden(ad_orig, resolution=0.5, random_state=42)
    adata_dict['Original'] = ad_orig

    # ---- CCA ----
    ad_cca = adata.copy()
    ad_cca.obsm['X_seurat_cca'] = integration_results['cca_embeddings']
    sc.pp.neighbors(ad_cca, use_rep='X_seurat_cca', random_state=42)
    sc.tl.leiden(ad_cca, resolution=0.5, random_state=42)
    sc.tl.umap(ad_cca, random_state=42)
    adata_dict['Seurat_CCA'] = ad_cca

    # ---- RPCA ----
    ad_rpca = adata.copy()
    ad_rpca.obsm['X_seurat_rpca'] = integration_results['rpca_embeddings']
    sc.pp.neighbors(ad_rpca, use_rep='X_seurat_rpca', random_state=42)
    sc.tl.leiden(ad_rpca, resolution=0.5, random_state=42)
    sc.tl.umap(ad_rpca, random_state=42)
    adata_dict['Seurat_RPCA'] = ad_rpca

    # ---- LIGER ----
    ad_liger = adata.copy()
    ad_liger.obsm['X_liger'] = integration_results['liger_embeddings']
    sc.pp.neighbors(ad_liger, use_rep='X_liger', random_state=42)
    sc.tl.leiden(ad_liger, resolution=0.5, random_state=42)
    sc.tl.umap(ad_liger, random_state=42)
    adata_dict['LIGER'] = ad_liger

    return adata_dict


# ================================================================
#  5. 评估指标
# ================================================================

def compute_metrics(adata_dict):
    """
    计算三组关键指标：
    - Batch Silhouette (在整合嵌入空间，越接近 0 表示批次混合越好)
    - N Clusters (leiden 聚类数)
    - ARI vs Original (与原始聚类的 Adjusted Rand Index)

    Returns
    -------
    dict[str, dict]
    """
    rep_keys = {
        'Original': 'X_pca',
        'Seurat_CCA': 'X_seurat_cca',
        'Seurat_RPCA': 'X_seurat_rpca',
        'LIGER': 'X_liger',
    }

    metrics = {}
    orig_leiden = adata_dict['Original'].obs.get('leiden')

    for name, ad in adata_dict.items():
        m = {}

        # Batch Silhouette (在嵌入空间中计算)
        rep_key = rep_keys.get(name, 'X_pca')
        if rep_key in ad.obsm and ad.obs['batch'].nunique() > 1:
            try:
                m['batch_silhouette'] = silhouette_score(
                    ad.obsm[rep_key],
                    ad.obs['batch'].astype(str).values,
                )
            except Exception:
                m['batch_silhouette'] = float('nan')
        else:
            m['batch_silhouette'] = float('nan')

        # Leiden 聚类数
        m['n_clusters'] = int(ad.obs['leiden'].nunique()) if 'leiden' in ad.obs else 0

        # ARI vs Original
        if orig_leiden is not None and 'leiden' in ad.obs:
            m['ari_vs_original'] = adjusted_rand_score(
                orig_leiden.values, ad.obs['leiden'].values
            )
        else:
            m['ari_vs_original'] = float('nan')

        metrics[name] = m

    return metrics


# ================================================================
#  6. 可视化
# ================================================================

def _batch_color_map(adata):
    """为批次标签构建颜色映射。"""
    batches = sorted(adata.obs['batch'].astype(str).unique())
    palette = plt.cm.Set1.colors
    return {b: palette[i % len(palette)] for i, b in enumerate(batches)}


def generate_plots(adata_dict, metrics, timings, output_dir):
    """生成三张对比图: UMAP-batch / UMAP-cluster / 指标柱状图。"""
    os.makedirs(output_dir, exist_ok=True)
    methods = list(adata_dict.keys())
    n = len(methods)
    bar_colors = ['#7f8c8d', '#2980b9', '#27ae60', '#8e44ad'][:n]

    # 统一的批次颜色
    cmap_batch = _batch_color_map(adata_dict['Original'])

    # ---- 图 1: UMAP 按批次着色 ----
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))
    if n == 1:
        axes = [axes]
    for ax, name in zip(axes, methods):
        ad = adata_dict[name]
        batch_str = ad.obs['batch'].astype(str)
        for b, color in cmap_batch.items():
            mask = batch_str == b
            ax.scatter(
                ad.obsm['X_umap'][mask, 0],
                ad.obsm['X_umap'][mask, 1],
                c=[color], label=f'Batch {b}', s=2, alpha=0.4, rasterized=True,
            )
        ax.set_title(name, fontsize=14, fontweight='bold')
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.legend(markerscale=4, fontsize=8, loc='best')
    fig.suptitle('Batch Distribution (colored by batch)', fontsize=16, y=1.02)
    plt.tight_layout()
    path1 = os.path.join(output_dir, 'comparison_umap_batch.png')
    fig.savefig(path1, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  保存: {path1}")

    # ---- 图 2: UMAP 按聚类着色 ----
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))
    if n == 1:
        axes = [axes]
    for ax, name in zip(axes, methods):
        ad = adata_dict[name]
        clusters = ad.obs['leiden'].astype(int) if 'leiden' in ad.obs else np.zeros(ad.n_obs)
        sc_plot = ax.scatter(
            ad.obsm['X_umap'][:, 0],
            ad.obsm['X_umap'][:, 1],
            c=clusters, cmap='tab20', s=2, alpha=0.4, rasterized=True,
        )
        n_clust = ad.obs['leiden'].nunique() if 'leiden' in ad.obs else 0
        ax.set_title(f'{name}  ({n_clust} clusters)', fontsize=14, fontweight='bold')
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle('Cluster Distribution (colored by leiden)', fontsize=16, y=1.02)
    plt.tight_layout()
    path2 = os.path.join(output_dir, 'comparison_umap_clusters.png')
    fig.savefig(path2, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  保存: {path2}")

    # ---- 图 3: 指标对比柱状图 ----
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Batch Silhouette
    vals = [metrics[m].get('batch_silhouette', 0) for m in methods]
    axes[0].bar(methods, vals, color=bar_colors)
    axes[0].set_title('Batch Silhouette Score\n(closer to 0 → better mixing)')
    axes[0].set_ylabel('Silhouette')
    axes[0].axhline(y=0, color='gray', ls='--', alpha=0.5)
    for i, v in enumerate(vals):
        axes[0].text(i, v + 0.005, f'{v:.3f}', ha='center', fontsize=9)

    # N Clusters
    vals = [metrics[m].get('n_clusters', 0) for m in methods]
    axes[1].bar(methods, vals, color=bar_colors)
    axes[1].set_title('Number of Leiden Clusters')
    axes[1].set_ylabel('# Clusters')
    for i, v in enumerate(vals):
        axes[1].text(i, v + 0.3, str(v), ha='center', fontsize=9)

    # Runtime
    times = [timings.get(m, 0) for m in methods]
    axes[2].bar(methods, times, color=bar_colors)
    axes[2].set_title('Integration Runtime')
    axes[2].set_ylabel('Time (seconds)')
    for i, v in enumerate(times):
        axes[2].text(i, v + 0.2, f'{v:.1f}s', ha='center', fontsize=9)

    plt.tight_layout()
    path3 = os.path.join(output_dir, 'comparison_metrics.png')
    fig.savefig(path3, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  保存: {path3}")


# ================================================================
#  7. 对比报告
# ================================================================

def print_report(metrics, timings):
    """打印终端可读的对比报告。"""
    methods = list(metrics.keys())
    col_width = max(15, max(len(m) for m in methods) + 2)
    table_width = 28 + col_width * len(methods)

    print("\n" + "=" * table_width)
    print("  CCA vs RPCA vs LIGER 整合对比报告")
    print("=" * table_width)

    header = f"{'指标':<28}"
    for m in methods:
        header += f"{m:>{col_width}}"
    print(header)
    print("-" * table_width)

    row = f"{'Batch Silhouette':<28}"
    for m in methods:
        v = metrics[m].get('batch_silhouette', float('nan'))
        row += f"{v:>{col_width}.4f}"
    print(row)

    row = f"{'N Clusters (leiden)':<28}"
    for m in methods:
        v = metrics[m].get('n_clusters', 0)
        row += f"{v:>{col_width}d}"
    print(row)

    row = f"{'ARI vs Original':<28}"
    for m in methods:
        v = metrics[m].get('ari_vs_original', float('nan'))
        row += f"{v:>{col_width}.4f}"
    print(row)

    row = f"{'Runtime (s)':<28}"
    for m in methods:
        v = timings.get(m, 0)
        row += f"{v:>{col_width}.1f}"
    print(row)

    print("=" * table_width)

    # 对比总结
    integration_methods = [m for m in methods if m != 'Original']
    sils = {m: abs(metrics[m].get('batch_silhouette', 99)) for m in integration_methods}
    runtimes = {m: timings.get(m, 999) for m in integration_methods}

    best_mix = min(sils, key=sils.get)
    fastest = min(runtimes, key=runtimes.get)

    print("\n对比总结:")
    sil_detail = ", ".join(f"{m}={metrics[m].get('batch_silhouette', float('nan')):.4f}"
                           for m in integration_methods)
    print(f"  批次混合 : {best_mix} 最优 (Batch Silhouette: {sil_detail})")

    time_detail = ", ".join(f"{m}={runtimes[m]:.1f}s" for m in integration_methods)
    print(f"  运行速度 : {fastest} 最快 ({time_detail})")

    clust_detail = ", ".join(f"{m}={metrics[m].get('n_clusters', '?')}"
                             for m in integration_methods)
    print(f"  聚类数量 : {clust_detail}")
    print()


# ================================================================
#  主流程
# ================================================================

def main():
    parser = argparse.ArgumentParser(
        description='CCA vs RPCA vs LIGER 批次整合对比分析')
    parser.add_argument(
        '--data', default='data/pbmc4k.h5ad',
        help='输入 h5ad 数据路径 (默认: data/pbmc4k.h5ad)')
    parser.add_argument(
        '--output', default='results/seurat_comparison',
        help='输出目录 (默认: results/seurat_comparison)')
    args = parser.parse_args()

    print("=" * 72)
    print("  CCA vs RPCA vs LIGER 整合对比分析")
    print("=" * 72)

    # ---- Step 1: 环境检查 ----
    print("\n[1/6] 检查运行环境 ...")
    check_r_environment()

    # ---- Step 2: 加载数据 ----
    print("\n[2/6] 加载数据 ...")
    adata = load_data(args.data)

    # ---- Step 3: R 整合 (CCA + RPCA + LIGER) ----
    print("\n[3/6] 运行 R 整合 (CCA + RPCA + LIGER) ...")
    results = run_r_integration(adata)

    # ---- Step 4: Python 后处理 ----
    print("\n[4/6] 构建对比数据 (neighbors → leiden → UMAP) ...")
    adata_dict = build_adata_variants(adata, results)

    timings = {
        'Original': 0.0,
        'Seurat_CCA': results['cca_time'],
        'Seurat_RPCA': results['rpca_time'],
        'LIGER': results['liger_time'],
    }

    # ---- Step 5: 计算评估指标 ----
    print("\n[5/6] 计算评估指标 ...")
    metrics = compute_metrics(adata_dict)

    # ---- Step 6: 可视化与报告 ----
    print("\n[6/6] 生成可视化图表 ...")
    generate_plots(adata_dict, metrics, timings, args.output)
    print_report(metrics, timings)

    print(f"所有结果已保存至: {os.path.abspath(args.output)}/")
    print("分析完成!\n")


if __name__ == '__main__':
    main()
