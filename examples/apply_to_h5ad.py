"""Minimal example: apply SODA to a real h5ad dataset."""

from pathlib import Path

import scanpy as sc

from soda import SODA, integrate_adata


INPUT_PATH = Path("path/to/real_dataset.h5ad")
OUTPUT_PATH = Path("path/to/real_dataset_soda.h5ad")


def main() -> None:
    adata = sc.read_h5ad(INPUT_PATH)

    integrator = SODA(
        n_clusters=50,
        max_iter=10,
        epsilon=0.05,
        tau=0.05,
        alpha=0.7,
    )

    corrected = integrate_adata(
        adata=adata,
        integrator=integrator,
        batch_key="batch",
        use_rep="X_pca",
        output_key="X_soda",
        copy=True,
    )
    corrected.write_h5ad(OUTPUT_PATH)


if __name__ == "__main__":
    main()
