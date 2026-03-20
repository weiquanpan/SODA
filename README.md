# SODA

SODA is a lightweight Python package for single-cell batch integration based on Sinkhorn Optimal Transport with Divergence-Based Adaptation.

This repository contains a clean, package-ready implementation extracted from the paper revision codebase, with a small public API, a command-line interface, and runnable examples for `.h5ad` datasets.

## Why This Repository

- Clean separation between package code and paper-specific benchmarking scripts
- Standard `src/` layout for packaging and GitHub release
- Minimal public API for Python and command-line use
- Ready-to-extend examples for lightweight demos and user-provided single-cell datasets

## Installation

Clone the repository and install it in editable mode:

```bash
git clone <your-repo-url>
cd SODA
pip install -e .
```

For development tools:

```bash
pip install -e .[dev]
```

## Quick Start

### Command line

Run SODA on any `.h5ad` file:

```bash
soda-run --input path/to/data.h5ad --batch-key batch
```

This writes a corrected file next to the input by default:

```text
path/to/data_soda.h5ad
```

### Python API

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

## Input Expectations

SODA expects:

- an AnnData object stored in `.h5ad` format
- a batch label column in `adata.obs`, typically `batch`
- a low-dimensional representation in `adata.obsm`, typically `X_pca`

If `X_pca` is missing, the package will compute PCA automatically.

## GitHub-Friendly Demo

A lightweight reproducible demo is included in [examples/README.md](/e:/E/W-test/compare/SODA/examples/README.md).

The recommended demo workflow is:

- [prepare_demo_data.R](/e:/E/W-test/compare/SODA/examples/prepare_demo_data.R)
- [run_real_data_demo.py](/e:/E/W-test/compare/SODA/examples/run_real_data_demo.py)

It is designed for GitHub use, so the repository does not need to ship large binary datasets. The `R` helper script creates a small `.h5ad` demo dataset from package-provided single-cell data, and the Python demo then runs SODA on that generated file.

For example:

```bash
Rscript examples/prepare_demo_data.R
python examples/run_real_data_demo.py --input outputs/demo_data/pbmc_small_demo.h5ad --batch-key batch --celltype-key cell_type
```

The repository also includes generic workflows for user-provided real `.h5ad` data. Together, these examples show how to:

- prepare or read an `.h5ad` dataset
- run SODA integration
- save the corrected AnnData object
- export UMAP figures before and after correction
- optionally compare batch and cell-type views

## Repository Layout

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

## Public API

The main public objects are:

- `soda.SODA`
- `soda.SODAResult`
- `soda.integrate_adata`

## Development

Basic local checks:

```bash
python -m compileall src examples
```

If you add tests later:

```bash
pytest
```

## Notes Before Public GitHub Release

This repository is now structurally ready for GitHub, but before public release you should still decide:

- the final repository name
- license choice
- author metadata
- project URL fields in `pyproject.toml`
- whether to add a `CITATION.cff` file

A short checklist is included in [docs/GITHUB_RELEASE_CHECKLIST.md](/e:/E/W-test/compare/SODA/docs/GITHUB_RELEASE_CHECKLIST.md).

The repository intentionally avoids committing large example datasets; instead, the demo data can be generated locally from the included example workflow.

## Citation

If you use SODA in this repository, please cite the accompanying manuscript once the final bibliographic information is fixed.
