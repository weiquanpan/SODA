# -*- coding: utf-8 -*-
"""CellANOVA 批次校正方法。"""

import functools as ft

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse
from scipy.sparse.linalg import svds
from sklearn import linear_model

from soda_benchmark.methods.base import BaseBatchMethod, register_method


@register_method("CellANOVA")
class CellANOVAMethod(BaseBatchMethod):
    """
    CellANOVA integration with automatic control pool construction.

    Notes
    -----
    - integrate_key 默认使用 'batch'
    - 若未提供 control_batches，则默认将所有批次作为 control pool
    """

    def __init__(
        self,
        integrate_key="batch",
        control_batches=None,
        control_pool_name="pool1",
        **kwargs,
    ):
        self.integrate_key = integrate_key
        self.control_batches = control_batches
        self.control_pool_name = control_pool_name

    @classmethod
    def check_dependencies(cls):
        try:
            import cellanova  # noqa: F401
            return True, ""
        except ImportError:
            return (
                False,
                "cellanova not installed. Install with: pip install cellanova",
            )

    def _build_control_dict(self, adata):
        all_batches = list(adata.obs[self.integrate_key].astype(str).unique())

        if self.control_batches is None:
            selected = all_batches
        else:
            selected = [str(x) for x in self.control_batches]
            missing = [b for b in selected if b not in all_batches]
            if missing:
                raise ValueError(
                    f"control_batches contains unknown batches: {missing}. "
                    f"Available batches: {all_batches}"
                )

        if len(selected) == 0:
            raise ValueError("control pool is empty. Please provide valid control batches.")

        return {self.control_pool_name: selected}

    def _to_contiguous_dense(self, x, dtype=np.float64):
        arr = np.asarray(x, dtype=dtype)
        # Torch/harmonypy cannot handle negative strides from some numpy views.
        if np.any(np.array(arr.strides) < 0) or not arr.flags["C_CONTIGUOUS"]:
            arr = np.ascontiguousarray(arr.copy(), dtype=dtype)
        return arr

    def _ensure_no_negative_strides(self, adata):
        if sparse.issparse(adata.X):
            adata.X = adata.X.tocsr(copy=True)
        else:
            adata.X = self._to_contiguous_dense(adata.X, dtype=np.float64)

        for layer_key in list(adata.layers.keys()):
            layer_val = adata.layers[layer_key]
            if sparse.issparse(layer_val):
                adata.layers[layer_key] = layer_val.tocsr(copy=True)
            else:
                adata.layers[layer_key] = self._to_contiguous_dense(layer_val, dtype=np.float64)

        for key in list(adata.obsm.keys()):
            val = adata.obsm[key]
            if sparse.issparse(val):
                adata.obsm[key] = val.tocsr(copy=True)
            else:
                adata.obsm[key] = self._to_contiguous_dense(val, dtype=np.float64)
        return adata

    def _calc_me_compatible(self, adata, k_harmony=70, dim_me=70):
        """
        Compute C/M/main_effect with harmonypy API compatibility.

        The upstream cellanova implementation uses scanpy's harmony wrapper.
        New harmonypy versions can return Z_corr with shape (n_cells, n_pcs),
        while old versions return (n_pcs, n_cells). We normalize both cases.
        """
        rna_temp = ad.AnnData(
            adata[:, adata.var.highly_variable].X,
            dtype=np.float32,
        )
        rna_temp.var_names = adata[:, adata.var.highly_variable].var_names
        rna_temp.obs = adata.obs.copy()

        n_comps = max(2, min(k_harmony, rna_temp.n_vars, rna_temp.n_obs - 1))
        sc.tl.pca(rna_temp, n_comps=n_comps, random_state=42)
        rna_temp.obsm["X_pca"] = self._to_contiguous_dense(
            rna_temp.obsm["X_pca"], dtype=np.float64
        )

        import harmonypy as hm

        ho = hm.run_harmony(
            rna_temp.obsm["X_pca"],
            rna_temp.obs,
            self.integrate_key,
            max_iter_harmony=30,
        )
        z_corr = ho.Z_corr
        if hasattr(z_corr, "cpu"):
            z_corr = z_corr.cpu().detach().numpy()
        z_corr = self._to_contiguous_dense(z_corr, dtype=np.float64)

        if z_corr.shape[0] == rna_temp.n_obs:
            x_pca_harmony = z_corr
        else:
            x_pca_harmony = z_corr.T

        adata_harmony = ad.AnnData(
            x_pca_harmony @ rna_temp.varm["PCs"].T,
            dtype=np.float32,
        )
        adata_harmony.var_names = rna_temp.var_names
        adata_harmony.obs = rna_temp.obs.copy()

        dim_use = max(1, min(dim_me, adata_harmony.n_obs - 1, adata_harmony.n_vars - 1))
        cmat, _, _ = svds(adata_harmony.X, k=dim_use)

        dataidx = np.unique(adata.obs[self.integrate_key])
        bloc_coef = ()
        for x in dataidx:
            mask = (adata.obs[self.integrate_key] == x).to_numpy()
            bloc_temp = linear_model.LinearRegression(fit_intercept=False)
            y_scale = self._to_contiguous_dense(adata.layers["scale"][mask, :], dtype=np.float64)
            bloc_temp.fit(cmat[mask, :], y_scale)
            bloc_coef = bloc_coef + (bloc_temp.coef_.T,)

        mmat = ft.reduce(np.add, bloc_coef) / len(dataidx)

        adata.obsm["Cmat"] = cmat
        adata.varm["Mmat"] = mmat.T
        adata.layers["main_effect"] = cmat @ mmat
        return adata

    def _prepare_cellanova_input(self, adata):
        """
        Prepare a non-negative, finite matrix for CellANOVA preprocessing.

        CellANOVA expects raw count-like input. For robustness:
        1) prefer `adata.raw`
        2) then `adata.layers['counts']`
        3) fallback to `adata.X`
        and sanitize NaN/Inf/negative values.
        """
        adata_in = adata.copy()

        # Prefer raw counts aligned to current features.
        if adata.raw is not None:
            try:
                raw_adata = adata.raw.to_adata()
                common = [g for g in adata.var_names if g in raw_adata.var_names]
                if len(common) > 0:
                    adata_in = adata[:, common].copy()
                    adata_in.X = raw_adata[:, common].X.copy()
                    adata_in.raw = None
                    return self._sanitize_matrix(adata_in)
            except Exception:
                pass

        # Next prefer explicit counts layer.
        if "counts" in adata.layers:
            adata_in.X = adata.layers["counts"].copy()
            adata_in.raw = None
            return self._sanitize_matrix(adata_in)

        # Fallback: use adata.X directly.
        adata_in.raw = None
        return self._sanitize_matrix(adata_in)

    def _sanitize_matrix(self, adata):
        X = adata.X
        if sparse.issparse(X):
            X = X.tocsr(copy=True)
            X.data = np.nan_to_num(X.data, nan=0.0, posinf=0.0, neginf=0.0)
            X.data[X.data < 0] = 0.0

            row_sums = np.asarray(X.sum(axis=1)).ravel()
            zero_rows = np.where(row_sums <= 0)[0]
            if len(zero_rows) > 0 and X.shape[1] > 0:
                X = X.tolil()
                X[zero_rows, 0] = 1.0
                X = X.tocsr()
            adata.X = X
        else:
            X = np.asarray(X, dtype=np.float64)
            X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
            X[X < 0] = 0.0
            row_sums = X.sum(axis=1)
            zero_rows = row_sums <= 0
            if np.any(zero_rows) and X.shape[1] > 0:
                X[zero_rows, 0] = 1.0
            adata.X = X

        return adata

    def apply(self, adata):
        print("\nApplying CellANOVA...")
        try:
            import cellanova as cnova

            if self.integrate_key not in adata.obs:
                if self.integrate_key == "batch" and "batch" in adata.obs:
                    pass
                else:
                    raise KeyError(
                        f"'{self.integrate_key}' not found in adata.obs. "
                        "Please set integrate_key to an existing column."
                    )

            adata.obs[self.integrate_key] = adata.obs[self.integrate_key].astype(str)
            control_dict = self._build_control_dict(adata)

            adata_input = self._prepare_cellanova_input(adata)
            adata_input.obs[self.integrate_key] = adata_input.obs[self.integrate_key].astype(str)
            adata_input = self._ensure_no_negative_strides(adata_input)

            adata_prep = cnova.model.preprocess_data(
                adata_input,
                integrate_key=self.integrate_key,
            )
            adata_prep = self._ensure_no_negative_strides(adata_prep)
            adata_prep = self._calc_me_compatible(adata_prep)
            adata_prep = cnova.model.calc_BE(
                adata_prep, self.integrate_key, control_dict
            )
            adata_prep = cnova.model.calc_TE(adata_prep, self.integrate_key)

            if "denoised" not in adata_prep.layers:
                raise RuntimeError("CellANOVA output missing layer 'denoised'.")

            denoised = adata_prep.layers["denoised"]
            if sparse.issparse(denoised):
                denoised = denoised.tocsr()
            else:
                denoised = np.asarray(denoised, dtype=np.float64)

            # 使用 denoised 表达作为后续度量输入，同时构建低维嵌入用于邻居图。
            adata.layers["cellanova_denoised"] = denoised
            adata.X = denoised.copy()

            n_comps = max(2, min(50, adata.n_vars))
            sc.tl.pca(adata, n_comps=n_comps, random_state=42)
            adata.obsm["X_cellanova"] = adata.obsm["X_pca"].copy()

            n_pcs = min(30, adata.obsm["X_cellanova"].shape[1])
            sc.pp.neighbors(adata, use_rep="X_cellanova", n_pcs=n_pcs, random_state=42)
            sc.tl.leiden(adata, resolution=0.5, random_state=42)
            sc.tl.umap(adata, random_state=42)

            print("  ✓ CellANOVA completed")
            return adata
        except Exception as e:
            print(f"  ✗ CellANOVA failed: {e}")
            return None
