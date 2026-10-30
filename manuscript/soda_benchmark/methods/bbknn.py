# -*- coding: utf-8 -*-
"""BBKNN 批次校正方法。"""

import scanpy as sc

from soda_benchmark.methods.base import BaseBatchMethod, register_method


@register_method("BBKNN")
class BBKNNMethod(BaseBatchMethod):

    def __init__(self, **kwargs):
        pass

    @classmethod
    def check_dependencies(cls):
        try:
            import bbknn  # noqa: F401
            return True, ""
        except ImportError:
            return False, "bbknn not installed. Install with: pip install bbknn"

    def apply(self, adata):
        print("\nApplying BBKNN...")
        try:
            import bbknn
            sc.tl.pca(adata, n_comps=50, random_state=42)
            bbknn.bbknn(adata, batch_key='batch', neighbors_within_batch=3)
            sc.tl.leiden(adata, resolution=0.5, random_state=42)
            sc.tl.umap(adata, random_state=42)
            print("  ✓ BBKNN completed")
            return adata
        except Exception as e:
            print(f"  ✗ BBKNN failed: {e}")
            return None

