# SODA

SODA is a Python package for single-cell batch integration based on Sinkhorn Optimal Transport with Divergence-Based Adaptation.

This repository provides a clean, reusable implementation of SODA intended for readers of the accompanying manuscript who want to install the method, run it on AnnData objects, and reproduce a lightweight demo without downloading large benchmark files.

## Repository Contents

- `src/soda/`: package source code
- `examples/`: minimal usage examples and a lightweight demo workflow
- `docs/`: short maintenance notes for the repository

The repository is intentionally focused on the reusable method implementation and example workflows. Large manuscript-specific datasets and benchmark outputs are not bundled here.

## Installation

Clone the repository and install it in editable mode:

```bash
git clone https://github.com/weiquanpan/SODA.git
cd SODA
pip install -e .
```

Optional development tools:

```bash
pip install -e .[dev]
```

## Expected Input

SODA works with AnnData objects stored in `.h5ad` format.

Typical input requirements are:

- a batch column in `adata.obs`, usually `batch`
- a low-dimensional representation in `adata.obsm`, usually `X_pca`

If `X_pca` is not present, the package can compute PCA automatically.

## Command-Line Usage

The package installs a command-line entry point called `soda-run`.

Minimal example:

```bash
soda-run --input path/to/data.h5ad --batch-key batch
```

This writes a corrected file next to the input by default:

```text
path/to/data_soda.h5ad
```

Common options:

- `--input`: input `.h5ad` file
- `--output`: output `.h5ad` file
- `--batch-key`: batch label column in `adata.obs`
- `--use-rep`: embedding key in `adata.obsm`
- `--tau`, `--alpha`, `--epsilon`: SODA hyperparameters
- `--no-neighbors`, `--no-leiden`, `--no-umap`: disable downstream Scanpy steps

## Python Usage

```python
import scanpy as sc

from soda import SODA, integrate_adata

adata = sc.read_h5ad("path/to/data.h5ad")

integrator = SODA(
    n_clusters=50,
    max_iter=10,
    epsilon=0.05,
    tau=0.05,
    alpha=0.7,
)

adata_soda = integrate_adata(
    adata=adata,
    integrator=integrator,
    batch_key="batch",
    use_rep="X_pca",
    output_key="X_soda",
    copy=True,
)

adata_soda.write_h5ad("path/to/data_soda.h5ad")
```

## Output

By default, the package stores the corrected representation in:

- `adata.obsm["X_soda"]`

It also records run metadata in:

- `adata.uns["soda"]["X_soda"]`

If neighbor graph construction, Leiden clustering, and UMAP are enabled, the output AnnData object also contains the corresponding derived results.

## Lightweight Demo

To keep the repository small, the demo does not ship large binary datasets.

Instead, a lightweight reproducible example is provided in [`examples/`](examples/):

1. Generate a compact demo `.h5ad` file from package-provided single-cell data:

```bash
Rscript examples/prepare_demo_data.R
```

2. Run SODA on the generated file:

```bash
python examples/run_real_data_demo.py --input outputs/demo_data/pbmc_small_demo.h5ad --batch-key batch --celltype-key cell_type
```

This workflow produces:

- a demo input dataset
- a corrected `.h5ad` file
- a short summary text file
- UMAP figures colored by batch
- optional UMAP figures colored by cell type

The demo preparation step uses package-provided data from `SeuratObject::pbmc_small`, so it is suitable for a GitHub repository where large files should be avoided.

## Using Your Own Data

If you already have a real `.h5ad` dataset, you can skip the `R` demo-preparation step and run either:

```bash
soda-run --input your_data.h5ad --batch-key batch
```

or:

```bash
python examples/run_real_data_demo.py --input your_data.h5ad --batch-key batch
```

Additional example notes are available in [`examples/README.md`](examples/README.md).

## Package API

The main public objects are:

- `soda.SODA`
- `soda.SODAResult`
- `soda.integrate_adata`

## Repository Structure

```text
SODA/
|-- pyproject.toml
|-- README.md
|-- CONTRIBUTING.md
|-- docs/
|   `-- GITHUB_RELEASE_CHECKLIST.md
|-- examples/
|   |-- README.md
|   |-- apply_to_h5ad.py
|   |-- prepare_demo_data.R
|   `-- run_real_data_demo.py
`-- src/
    `-- soda/
        |-- __init__.py
        |-- __main__.py
        |-- _version.py
        |-- cli.py
        |-- model.py
        `-- scanpy.py
```

## Relation to the Manuscript

This repository contains the reusable SODA implementation and example workflows intended for software use and method illustration. It does not attempt to mirror the full manuscript revision workspace, benchmarking outputs, or all paper-specific analysis files.

## Citation

If you use SODA in your work, please cite the accompanying manuscript once the final bibliographic information is available.
