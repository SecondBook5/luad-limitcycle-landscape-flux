# limitcycle-landscape

A small CLI for treating the cell cycle as a limit cycle in single-cell data, then aligning other modalities (for example scATAC-seq) to the inferred phase.

This repo is meant to be practical: it ingests data into AnnData, infers a circular phase coordinate θ using S and G2/M programs, and bins per-cell signals by phase.

## What it does

1. **Ingest scRNA-seq into AnnData**
   - Builds a single `adata` with counts plus metadata (`timepoint`, `sample_id`, optional `genotype`)
   - Optionally imports author-provided embeddings if present

2. **Infer a limit-cycle phase**
   - Scores S and G2/M programs per cell
   - Fits a circle in one or more bases (configurable) and computes:
     - `theta` in `[0, 2π)`
     - `radius` and fit diagnostics
   - Persists results into `adata.obs`

3. **Align and bin other signals by phase**
   - Given a per-cell matrix (for example ATAC accessibility), aggregates by θ bins
   - Saves a phase-resolved pseudobulk object into `adata.uns` and exports tables
---

## Requirements

- Python 3.10+ recommended
- Works on Linux, macOS, and WSL2
- If using large sparse matrices, you will want enough RAM for PCA and OT-style operations

---

## Installation

```bash
python -m venv .venv
source .venv/bin/activate

pip install -U pip wheel setuptools
pip install -r requirements.txt
pip install -e .
````

Sanity check:

```bash
limitcycle-landscape --help
```

---

## Data layout

By default the repo expects:

```
data/
  raw/
    <dataset_name>/
      ... your raw files ...
  interim/
  processed/
outputs/
  figures/
  tables/
  logs/
configs/
  params.yaml
```

You control paths via `configs/params.yaml`.

---

## Configuration

Edit `configs/params.yaml` to point to your raw data and choose options.

Minimal example:

```yaml
paths:
  raw_root: "data/raw"
  interim_root: "data/interim"
  processed_root: "data/processed"
  outputs_root: "outputs"

dataset:
  name: "gse154989"          # or "generic_h5ad"
  input_h5ad: null           # only used when name == "generic_h5ad"

preprocess:
  min_genes_per_cell: 200
  min_cells_per_gene: 3
  n_hvgs: 2000
  pca_n_components: 50

cycle:
  s_genes_file: "configs/genes_s.txt"
  g2m_genes_file: "configs/genes_g2m.txt"
  score_method: "scanpy"     # implementation-dependent
  circle_fit_basis: "S_G2M"  # options depend on code, keep default unless you know why
  theta_key: "theta"
  radius_key: "radius"

atac:
  default_layer: "X_atac"
  n_bins: 60
  output_uns_key: "atac_pseudobulk"
```

Notes:

* If you are not using the built-in ingestion adapter, set `dataset.name: generic_h5ad` and provide `dataset.input_h5ad` pointing to an existing `.h5ad` that already contains `obs` metadata and `X` counts.
* The ATAC binning step assumes your ATAC signal is already aligned to the same cells as the scRNA object. If it is in a separate AnnData, you need a merge step before `bin-atac`.

---

## Quickstart (end to end)

1. Put raw data where `paths.raw_root` points.
2. Run ingestion.
3. Fit the cycle.
4. Bin ATAC (or any other layer) by phase.

```bash
limitcycle-landscape ingest-scrna
limitcycle-landscape fit-cycle
limitcycle-landscape bin-atac --layer X_atac --bins 60
```

Expected artifacts:

* `data/interim/mm_timecourse.h5ad` (or dataset-specific name)
* `data/processed/mm_timecourse_with_theta.h5ad`
* `outputs/tables/atac_phase_bins.csv`
* `outputs/figures/cycle_fit.png`
* `outputs/figures/theta_umap.png`

If your CLI supports a config override, use:

```bash
limitcycle-landscape ingest-scrna --config configs/params.yaml
```

---

## CLI overview

| Command        | Output                             | Purpose                                      |
| -------------- | ---------------------------------- | -------------------------------------------- |
| `ingest-scrna` | `data/interim/*.h5ad`              | Build AnnData with metadata + counts.        |
| `fit-cycle`    | `data/processed/*_with_theta.h5ad` | Score S/G2M, fit circle, write `obs[theta]`. |
| `bin-atac`     | tables + `adata.uns[...]`          | Aggregate a layer by θ bins.                 |

Run `limitcycle-landscape <command> --help` to see the arguments supported by your local version.

---

## How phase inference works (high level)

1. Compute module scores for S genes and G2/M genes.
2. Optionally z-score the scores and/or project into a basis used for circle fitting.
3. Fit a circle to obtain a phase angle θ per cell.
4. Store:

   * `adata.obs["theta"]` in radians `[0, 2π)`
   * `adata.obs["radius"]`
   * Fit diagnostics in `adata.uns["cycle_fit"]` (if implemented)

This is meant to be compatible with a landscape and flux view of nonequilibrium dynamics:

* The inferred θ is a coordinate along a cyclic trajectory (limit cycle).
* Phase-binned signals can be interpreted as observables along that cyclic flow.

---

## Notes on ATAC alignment

`bin-atac` expects an AnnData **layer** containing a matrix shaped `(n_cells, n_features)` aligned to the same cells used for θ.

Common patterns:

* `adata.layers["X_atac"]` is sparse and has the same row order as `adata.X`
* Peaks or topics live in `adata.var_names` or in a separate feature index

The command will:

1. Bin cells by θ into `n_bins` equal-width bins.
2. Aggregate ATAC signal within each bin (mean or sum, depending on implementation).
3. Save:

   * bin edges
   * counts per bin
   * aggregated matrix

If your ATAC data are in a separate object, merge them first so that the row order and cell IDs match.

---

## Using this on the HPCS LUAD dataset (Marjanovic 2020)

This repo can be applied to any time course dataset that contains cycling cells, including the Marjanovic 2020 LUAD evolution dataset, but only if:

* you subset to epithelial tumor cells that are cycling, or at least not dominated by non-cycling programs
* you validate that S and G2/M programs show the expected periodic structure

Practical advice:

* Run `fit-cycle` on a cycling-enriched subset.
* Use the phase coordinate θ as an alignment axis for comparing accessibility, TF activity, or other modalities within the cycling compartment.

Do not claim “limit cycle” unless the fitted circle is supported by marker periodicity and stable across samples.

---

## Troubleshooting

### `ModuleNotFoundError: limitcycle_landscape`

* You likely forgot `pip install -e .` or your venv is not active.

### `KeyError: X_atac` when running `bin-atac`

* Confirm the layer exists:

  * In Python: `print(adata.layers.keys())`
* Or pass the correct layer name:

  * `limitcycle-landscape bin-atac --layer <your_layer>`

### θ looks random

* Check that your S/G2M gene lists match the organism and gene naming convention.
* Check you have enough cycling cells. If most cells are quiescent, circle fitting will fail silently or produce nonsense.
* Plot S score vs G2/M score before fitting. You should see a roughly circular or annular structure.

---

## References (conceptual)

Zhu and Wang (2025), landscape and flux perspective:

```bibtex
@article{zhu2025quantifying,
  title={Quantifying landscape and flux from single-cell omics: unraveling the physical mechanisms of cell function},
  author={Zhu, Ligang and Wang, Jin},
  journal={JACS Au},
  volume={5},
  number={8},
  pages={3738--3757},
  year={2025},
  publisher={ACS Publications}
}
```

Marjanovic et al. (2020), LUAD evolution and HPCS:

```bibtex
@article{marjanovic2020emergence,
  title={Emergence of a high-plasticity cell state during lung cancer evolution},
  author={Marjanovic, Nemanja Despot and Hofree, Matan and Chan, Jason E and others},
  journal={Cancer Cell},
  volume={38},
  number={2},
  pages={229--246},
  year={2020},
  publisher={Elsevier}
}

