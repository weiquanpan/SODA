# -*- coding: utf-8 -*-
"""scVI 批次校正方法。"""

import scanpy as sc
from scipy import sparse

from soda_benchmark.methods.base import BaseBatchMethod, register_method


@register_method("scVI")
class ScVIMethod(BaseBatchMethod):

    def __init__(self, max_epochs=100, **kwargs):
        self.max_epochs = max_epochs

    @classmethod
    def check_dependencies(cls):
        try:
            import scvi  # noqa: F401
            return True, ""
        except ImportError:
            return False, "scvi-tools not installed. Install with: pip install scvi-tools"

    def apply(self, adata):
        print("\nApplying scVI...")
        try:
            import scvi

            adata.obs['batch'] = adata.obs['batch'].astype(str)

            if adata.raw is not None:
                adata_full = adata.raw.to_adata()
                hvg_genes = adata.var_names
                if 'counts' in adata_full.layers:
                    adata.layers['counts'] = adata_full[:, hvg_genes].layers['counts']
                else:
                    adata.layers['counts'] = adata_full[:, hvg_genes].X.copy()
            elif 'counts' not in adata.layers:
                print("  Warning: No raw counts available, using approximation")
                adata.layers['counts'] = adata.X.copy()

            if not sparse.issparse(adata.layers['counts']):
                adata.layers['counts'] = sparse.csr_matrix(adata.layers['counts'])

            scvi.settings.seed = 42
            scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='batch')
            model = scvi.model.SCVI(adata, n_latent=30, n_layers=2)

            try:
                model.train(max_epochs=self.max_epochs, train_size=0.9,
                            early_stopping=True, check_val_every_n_epoch=10,
                            accelerator='auto')
            except Exception:
                model.train(max_epochs=self.max_epochs, train_size=0.9,
                            early_stopping=True, check_val_every_n_epoch=10,
                            use_gpu=False)

            adata.obsm['X_scVI'] = model.get_latent_representation()
            sc.pp.neighbors(adata, use_rep='X_scVI', random_state=42)
            sc.tl.leiden(adata, resolution=0.5, random_state=42)
            sc.tl.umap(adata, random_state=42)
            print("  ✓ scVI completed")
            return adata
        except Exception as e:
            print(f"  ✗ scVI failed: {e}")
            return None

