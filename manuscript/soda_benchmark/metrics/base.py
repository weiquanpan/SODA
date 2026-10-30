# -*- coding: utf-8 -*-
"""
评估指标的抽象基类、标准化结果结构与注册表。

扩展方式
--------
1. 在 soda_benchmark/metrics/ 目录下新建 Python 文件
2. 继承 BaseMetric 实现 calculate / grade
3. 使用 @register_metric("指标键名") 装饰器注册

框架会自动发现并调用注册的指标，无需修改其他文件。
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

# --------------- 全局注册表 ---------------
_METRIC_REGISTRY: Dict[str, type] = {}


def register_metric(key: str):
    """装饰器：将指标类注册到全局注册表。"""

    def decorator(cls):
        _METRIC_REGISTRY[key] = cls
        return cls

    return decorator


def get_registered_metrics() -> Dict[str, type]:
    """返回当前所有已注册的指标（键名 → 类）。"""
    return dict(_METRIC_REGISTRY)


# --------------- 标准化结果结构 ---------------

@dataclass
class MetricResult:
    """统一的指标计算结果。

    Attributes
    ----------
    name : str
        指标的显示名称，例如 "Neighborhood Preservation"
    key : str
        注册键名，例如 "neighborhood_preservation"
    higher_is_better : bool
        指标方向：True = 越大越好，False = 越小越好
    scores : pd.DataFrame
        索引为方法名，列固定为 ['score', 'ci_lower', 'ci_upper']
    extra : dict
        附加数据（如 neighborhood 指标的 per_cell 分布）
    """
    name: str
    key: str
    higher_is_better: bool
    scores: pd.DataFrame
    extra: dict = field(default_factory=dict)


# --------------- 抽象基类 ---------------

class BaseMetric(ABC):
    """所有评估指标的抽象基类。"""

    # 子类必须定义
    name: str = ""
    higher_is_better: bool = True

    def __init__(self, n_bootstrap: int = 200, ci_level: int = 95,
                 random_state: int = 42, **kwargs):
        self.n_bootstrap = n_bootstrap
        self.ci_level = ci_level
        self.rng = np.random.default_rng(random_state)

    @abstractmethod
    def calculate(self, adata_original, adata_dict) -> MetricResult:
        """
        计算指标。

        Parameters
        ----------
        adata_original : AnnData
            原始（未校正）数据
        adata_dict : dict
            方法名 → 校正后 AnnData 的字典

        Returns
        -------
        MetricResult
        """
        pass

    @abstractmethod
    def grade(self, value: float) -> str:
        """
        为单个数值打等级标签。

        Returns
        -------
        "Excellent" / "Good" / "Fair" / "Poor"
        """
        pass

    # ------------------------------------------------------------------ #
    # 通用 bootstrap 工具
    # ------------------------------------------------------------------ #

    def _bootstrap_ci(self, values: np.ndarray,
                      reducer=np.mean) -> Tuple[float, float]:
        """对一维数组做 bootstrap 置信区间。"""
        stats = []
        for _ in range(self.n_bootstrap):
            sample = self.rng.choice(values, size=len(values), replace=True)
            stats.append(reducer(sample))
        return self._ci_bounds(np.array(stats))

    def _bootstrap_paired_ci(self, a: np.ndarray, b: np.ndarray,
                             func) -> Tuple[float, float]:
        """对成对数据做 bootstrap 置信区间。"""
        stats = []
        n = len(a)
        for _ in range(self.n_bootstrap):
            idx = self.rng.integers(0, n, size=n)
            stats.append(func(a[idx], b[idx]))
        return self._ci_bounds(np.array(stats))

    def _ci_bounds(self, samples: np.ndarray) -> Tuple[float, float]:
        alpha = (100 - self.ci_level) / 2
        lower = np.percentile(samples, alpha)
        upper = np.percentile(samples, 100 - alpha)
        return float(lower), float(upper)

