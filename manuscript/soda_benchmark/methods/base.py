# -*- coding: utf-8 -*-
"""
批次校正方法的抽象基类与注册表。

扩展方式
--------
1. 在 soda_benchmark/methods/ 目录下新建 Python 文件
2. 继承 BaseBatchMethod 实现 check_dependencies / apply
3. 使用 @register_method("方法名") 装饰器注册

框架会自动发现并调用注册的方法，无需修改其他文件。
"""

from abc import ABC, abstractmethod
from typing import Dict, Tuple

# --------------- 全局注册表 ---------------
_METHOD_REGISTRY: Dict[str, type] = {}


def register_method(name: str):
    """装饰器：将方法类注册到全局注册表。"""

    def decorator(cls):
        _METHOD_REGISTRY[name] = cls
        return cls

    return decorator


def get_registered_methods() -> Dict[str, type]:
    """返回当前所有已注册的方法（名称 → 类）。"""
    return dict(_METHOD_REGISTRY)


# --------------- 抽象基类 ---------------
class BaseBatchMethod(ABC):
    """所有批次校正方法的抽象基类。"""

    @classmethod
    @abstractmethod
    def check_dependencies(cls) -> Tuple[bool, str]:
        """
        检查该方法所需的依赖是否可用。

        Returns
        -------
        (available, message)
            available: 是否可用
            message: 不可用时的提示信息
        """
        return True, ""

    @abstractmethod
    def apply(self, adata):
        """
        执行批次校正。

        Parameters
        ----------
        adata : AnnData
            输入的 AnnData 对象（已经是 copy）。

        Returns
        -------
        AnnData or None
            校正后的 AnnData（含 leiden + UMAP），失败时返回 None。
        """
        pass

