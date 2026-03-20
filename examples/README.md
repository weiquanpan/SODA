# Examples

This folder contains lightweight examples that show how to apply SODA without committing large datasets to the repository.

## Files

- `prepare_demo_data.R`: creates a tiny GitHub-friendly `.h5ad` demo dataset from `SeuratObject::pbmc_small`
- `apply_to_h5ad.py`: smallest possible Python usage example
- `run_real_data_demo.py`: fuller demo runner with outputs and UMAP figures

## Before Running

Install the package from the repository root:

```bash
pip install -e .
```

## Lightweight Demo From Package Data

To avoid storing large data files in the repository, we provide a small demo dataset generated from `SeuratObject::pbmc_small`.

Step 1: generate the demo `.h5ad`

```bash
Rscript examples/prepare_demo_data.R
```

Step 2: run SODA on the generated file

```bash
python examples/run_real_data_demo.py --input outputs/demo_data/pbmc_small_demo.h5ad --batch-key batch --celltype-key cell_type
```

This produces a complete, reproducible demo without adding large binary data to GitHub.

## Generic `.h5ad` Workflow

Run:

```bash
python examples/run_real_data_demo.py --input path/to/real_dataset.h5ad
```

This script will:

- read the input AnnData object
- run SODA using the chosen batch column
- save a corrected `.h5ad` file
- save UMAP figures colored by batch
- optionally save UMAP figures colored by cell type

### Useful arguments

```bash
python examples/run_real_data_demo.py \
  --input path/to/real_dataset.h5ad \
  --output-dir outputs/pbmcsca_demo \
  --batch-key batch \
  --celltype-key CellType \
  --use-rep X_pca
```

## Expected Output

The demo writes files such as:

- `integrated.h5ad`
- `summary.txt`
- `umap_original_batch.png`
- `umap_soda_batch.png`
- `umap_original_celltype.png`
- `umap_soda_celltype.png`

## R Dependencies For Demo Preparation

The demo-data preparation script expects these R packages:

- `SeuratObject`
- `SingleCellExperiment`
- `zellkonverter`

These are only needed for generating the lightweight demo dataset. They are not required for the Python package itself.
