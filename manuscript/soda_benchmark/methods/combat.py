# -*- coding: utf-8 -*-
"""ComBat 批次校正方法。"""

import scanpy as sc

from soda_benchmark.methods.base import BaseBatchMethod, register_method


@register_method("ComBat")
class ComBatMethod(BaseBatchMethod):

    def __init__(self, **kwargs):
        pass

    @classmethod
    def check_dependencies(cls):
        # ComBat 内置于 scanpy，始终可用
        return True, ""

    def apply(self, adata):
        print("\nApplying ComBat...")
        try:
            sc.pp.combat(adata, key='batch')
            sc.pp.neighbors(adata, n_pcs=30, random_state=42)
            sc.tl.leiden(adata, resolution=0.5, random_state=42)
            sc.tl.umap(adata, random_state=42)
            print("  ✓ ComBat completed")
            return adata
        except Exception as e:
            print(f"  ✗ ComBat failed: {e}")
            return None

