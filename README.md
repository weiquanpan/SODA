# SODA
This repository contains the source code and reproduction scripts for the paper "SODA: Sinkhorn Optimal Transport with Divergence-Based Adaptation for Well-Calibrated Single-Cell Integration".

## 📖 论文信息

**SODA: Sinkhorn Optimal Transport with Divergence-Based Adaptation for Well-Calibrated Single-Cell Integration**

本仓库包含论文中 SODA 算法的完整实现及与其他主流方法的性能对比基准测试。

## 🌟 主要特性

- **SODA 批次校正算法**：基于 Sinkhorn 最优传输和散度自适应的批次效应校正
- **全面的基准测试**：与 Harmony、ComBat、BBKNN、scVI 等方法的系统性比较
- **多维度评估指标**：
  - 邻域保持性 (Neighborhood Preservation)
  - 聚类一致性 (Cluster Concordance)  
  - 调整兰德指数 (Adjusted Rand Index)
- **可视化结果**：自动生成 UMAP 可视化和性能指标图表

## 🚀 快速开始

### 环境要求

- Python 3.12
- 推荐使用 conda 或 virtualenv 创建独立环境



### 依赖包

主要依赖包包括：

- `scanpy` - 单细胞数据分析
- `numpy`, `pandas` - 数据处理
- `POT` - Python Optimal Transport 库
- `scikit-learn` - 机器学习工具
- `matplotlib`, `seaborn` - 可视化
- `bbknn` (可选) - BBKNN 批次校正
- `scvi-tools` (可选) - scVI 批次校正

## 📊 数据准备

### 数据格式

输入数据需要为 `.h5ad` 格式（AnnData 对象），并包含以下信息：

- **adata.X**: 基因表达矩阵
- **adata.obs['batch']**: 批次标签列
- **adata.obs['leiden']**: 原始聚类标签（可选，用于评估）

### 数据预处理

数据下载和预处理过程请参考：

- **数据来源**：见附件
- **预处理流程**：https://github.com/pmelsted/AM_2024

将原始数据转换为 h5ad 格式后，即可使用本仓库的算法进行批次校正。

### 示例数据结构

```python
import scanpy as sc

# 加载数据示例
adata = sc.read_h5ad('data/your_data.h5ad')
print(adata)
# AnnData object with n_obs × n_vars
#     obs: 'batch', 'leiden', ...
#     obsm: 'X_pca', 'X_umap'
```

## 💻 使用方法

### 基本用法

```bash
# 运行 SODA 批次校正及基准测试
python soda_benchmark.py --input data/your_data.h5ad
```

### 命令行参数

```bash
python soda_benchmark.py --help
```

**参数说明：**

| 参数           | 简写 | 类型 | 默认值    | 说明                       |
| -------------- | ---- | ---- | --------- | -------------------------- |
| `--input`      | `-i` | str  | 必需      | 输入 H5AD 文件路径         |
| `--output`     | `-o` | str  | `results` | 输出结果目录               |
| `--epochs`     | -    | int  | `100`     | scVI 训练的最大轮数        |
| `--batch-key`  | -    | str  | `batch`   | adata.obs 中批次信息的列名 |
| `--n-clusters` | -    | int  | `50`      | SODA 聚类数量              |
| `--max-iter`   | -    | int  | `10`      | SODA 最大迭代次数          |

### 使用示例

```bash
# 示例 1: 基本使用
python soda_benchmark.py --input data/simulated_data.h5ad

# 示例 2: 自定义输出目录和 scVI 训练轮数
python soda_benchmark.py --input data/pbmc4k-0.h5ad --epochs 200 --output results_pbmc

# 示例 3: 使用自定义批次列名
python soda_benchmark.py --input data/my_data.h5ad --batch-key batch_id

# 示例 4: 调整 SODA 参数
python soda_benchmark.py --input data/my_data.h5ad --n-clusters 100 --max-iter 15
```

## 📈 输出结果

运行完成后，在输出目录（默认 `results/`）中会生成：

### 可视化图表

- `umap_comparison.png` - 各方法的 UMAP 可视化（按批次着色）
- `umap_clusters.png` - 各方法的 UMAP 可视化（按聚类着色）
- `metrics_summary.png` - 性能指标综合对比
- `neighborhood_distribution.png` - 邻域保持性分布图

### 评估报告

- `benchmark_report.txt` - 详细的性能评估报告，包含：
  - 数据集信息
  - 评估方法列表
  - 各项指标得分
  - 综合排名

## 🔬 算法原理

SODA 采用 Harmony 的迭代聚类范式，但使用熵正则化最优传输替代线性质心校正：

1. **Sinkhorn 散度量化**：使用 Sinkhorn 散度 SD(P‖Q) = W²₂(P,Q) − ½W²₂(P,P) − ½W²₂(Q,Q) 量化批次效应
2. **自适应校正**：仅对散度超过阈值的聚类进行基于最优传输的对齐
3. **迭代优化**：通过迭代聚类和校正直至收敛


## 🙏 致谢

- 数据预处理参考：https://github.com/pmelsted/AM_2024
- 感谢 scanpy、POT 等开源项目的支持

---

**注意**：本仓库实现了论文中提出的 SODA 算法，数据预处理部分请参考上述链接。
