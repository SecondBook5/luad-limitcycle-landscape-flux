# limitcycle-landscape

Public-facing tool for evaluating the cell cycle as a limit cycle and aligning scRNA-seq with ATAC-seq
signals along the inferred phase. The implementation migrates the core notebook logic into a
reusable, documented Python package with a CLI.

## What this repository provides

* **Ingestion** of the GSE154989 scRNA-seq timecourse into an AnnData object, including sample metadata
  parsing and author-provided embeddings.
* **Limit-cycle phase inference** from S and G2/M marker modules, with multi-basis circle fitting and
  persisted `theta`/`radius` values for downstream analyses.
* **ATAC alignment hooks** that aggregate accessibility (scATAC layer) by inferred phase bins.

## Quickstart

1. Create an environment and install dependencies.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

2. Ensure the raw data are placed under the paths in `configs/params.yaml`.

3. Run the pipeline:

```bash
limitcycle-landscape ingest-scrna
limitcycle-landscape fit-cycle
limitcycle-landscape bin-atac --layer X_atac --bins 60
```

Artifacts are saved under `outputs/` and `data/interim/`. Configuration lives in `configs/params.yaml`.

## CLI overview

| Command | Purpose |
| --- | --- |
| `ingest-scrna` | Build `data/interim/mm_timecourse.h5ad` with metadata + counts. |
| `fit-cycle` | Compute z-scored S/G2M modules and fit the principal circle to obtain θ. |
| `bin-atac` | Aggregate ATAC signal by θ bins and store in `adata.uns["atac_pseudobulk"]`. |

## Package layout

```
src/limitcycle_landscape/
  cli.py              CLI entry point
  config.py           configuration helpers
  ingest.py           scRNA-seq ingestion
  cycle.py            limit-cycle phase inference + ATAC binning
  cc_limitcycle.py    math utilities (ported from notebooks/_lib)
```

## Notes on ATAC alignment

The `bin-atac` command expects an AnnData layer (default `X_atac`) containing cell-by-peak
accessibility aligned to the same cells used for phase inference. The output payload contains
bin edges, counts per bin, and the binned matrix, enabling downstream correlation with phase-resolved
scRNA expression.

## HPCS dynamics pipeline (Marjanovic 2020 LUAD time course)

This repository now includes an end-to-end HPCS dynamics model that reproduces the transition graph
and adds landscape/flux plus ML-driven dynamics modeling.

### One-command run

```bash
python -m hpcs_model.cli run --config configs/hpcs_default.yaml
```

### Outputs

* Preprocessed AnnData: `outputs/adata_preprocessed.h5ad`
* OT couplings: `outputs/tables/couplings/`
* Cluster transition graph: `outputs/tables/cluster_transition_edges.csv`
* HPCS selection + bootstrap: `outputs/tables/hpcs_selection.json`, `outputs/tables/hpcs_bootstrap.csv`
* Landscape + flux: `outputs/figures/potential_contours.png`, `outputs/figures/drift_quiver.png`,
  `outputs/figures/flux_residual_magnitude.png`
* ML model + evaluation: `outputs/models/drift_model.pt`, `outputs/tables/ml_eval.json`

### CLI stages

| Command | Purpose |
| --- | --- |
| `validate` | Load and validate required metadata columns. |
| `preprocess` | QC, normalize, HVG, PCA, neighbors, embedding. |
| `couple` | Compute OT couplings between consecutive timepoints. |
| `graph` | Aggregate couplings into a transition graph and identify HPCS. |
| `landscape_flux` | Estimate potential, drift, and flux residuals. |
| `train_ml` | Train drift regression model from OT-derived displacements. |
| `run` | Run the full pipeline end-to-end. |
