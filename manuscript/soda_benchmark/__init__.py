# -*- coding: utf-8 -*-
"""
SODA: Sinkhorn Optimal Transport with Divergence-Based Adaptation
for Well-Calibrated Single-Cell Integration

项目入口包，导出核心类供外部使用。
"""

from soda_benchmark.core.benchmark import SingleCellBatchBenchmark
from soda_benchmark.core.report import ReportGenerator

__all__ = [
    "SingleCellBatchBenchmark",
    "ReportGenerator",
]

