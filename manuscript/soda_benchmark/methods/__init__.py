# -*- coding: utf-8 -*-
"""
批次校正方法子包。

采用自动发现机制：只需在本目录下添加新的方法模块并使用 @register_method 装饰器，
即可被框架自动识别和调度，无需修改任何其他文件。
"""

import importlib
import pkgutil

from soda_benchmark.methods.base import (
    BaseBatchMethod,
    get_registered_methods,
    register_method,
)

# 自动导入本包内所有模块（除 base），触发 @register_method 注册
for _finder, _module_name, _ispkg in pkgutil.iter_modules(__path__):
    if _module_name != "base":
        importlib.import_module(f".{_module_name}", __package__)

__all__ = [
    "BaseBatchMethod",
    "get_registered_methods",
    "register_method",
]

