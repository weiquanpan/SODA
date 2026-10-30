# -*- coding: utf-8 -*-
"""
评估指标子包。

采用自动发现机制：只需在本目录下添加新的指标模块并使用 @register_metric 装饰器，
即可被框架自动识别和调度，无需修改任何其他文件。
"""

import importlib
import pkgutil

from soda_benchmark.metrics.base import (
    BaseMetric,
    MetricResult,
    get_registered_metrics,
    register_metric,
)

# 自动导入本包内所有模块（除 base），触发 @register_metric 注册
for _finder, _module_name, _ispkg in pkgutil.iter_modules(__path__):
    if _module_name != "base":
        importlib.import_module(f".{_module_name}", __package__)

__all__ = [
    "BaseMetric",
    "MetricResult",
    "get_registered_metrics",
    "register_metric",
]

