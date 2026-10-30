# Manuscript benchmark and analysis scripts

These scripts were supplied by the authors for the SODA manuscript. They are kept separate from the reusable package in `src/soda/` and retain their original directory structure. `source-manifest.json` records SHA-256 hashes of the supplied sources.

## Entry points

Run commands from this `manuscript/` directory, because several workflows use relative `data/` and `results/` paths.

- `run_benchmark.py`: run a single-dataset benchmark, for example `python run_benchmark.py --data data/pbmc4k.h5ad --output-dir results`.
- `run_batch.py` and `run_three_datasets.py`: batch benchmark workflows.
- `run_seurat_comparison.py`: Seurat comparison workflow.
- `run_scaling_profile.py`: runtime and memory profiling.
- `run_reviewer_all_sensitivity.py`: sensitivity analysis.
- `simulate_cellanova_data.py`, `simulate_sc3desigin_data.r`, and `sweep_scdesign3_complexity.py`: simulation workflows.
- `export_real_pbmcsca_subset_h5ad.R` and `run_reviewer_figure1_pbmcsca.py`: cross-technology PBMC preparation and analysis.
- `generate_rebuttal_tables.py`: assemble analysis summaries.
- `soda_benchmark/`: benchmark methods, metrics, plots, and report code used by these entry points.

## Environment

These are research scripts, not a separately packaged application. Python imports include NumPy, pandas, SciPy, scikit-learn, Scanpy, AnnData, matplotlib, seaborn, POT (`ot`), psutil, torch, and method-specific packages such as harmonypy, bbknn, cellanova, scvi-tools, and rpy2. R-based methods additionally require their R dependencies. Inspect the selected workflow for its exact inputs and optional method requirements. No verified environment lockfile was supplied with this source snapshot.

## Data and results

The original public dataset sources are listed in Supplementary Table S3 of the manuscript submission. This source upload does not include local `.h5ad` datasets, generated benchmark result directories, figures, training logs, or caches. A public location for manuscript-specific processed datasets and outputs remains to be specified. Do not interpret the presence of these scripts as verification that every published table and figure can already be reproduced from this repository alone.

## Validation

All included Python files passed a syntax parse. Full benchmarks were not rerun as part of this source upload, and output agreement with the manuscript has not been independently verified.
