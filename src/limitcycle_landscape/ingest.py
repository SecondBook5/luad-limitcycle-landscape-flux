from __future__ import annotations

from pathlib import Path
import re
from typing import Dict, Tuple

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from limitcycle_landscape.config import load_config, repo_root


SAMPLE_RE = re.compile(
    r"^(?P<genotype>KP|K)_(?P<week>\d+)w_(?P<treatment>[^_]+)_m(?P<mouse>\d+)"
    r"_T(?P<timepoint>\d+)_P(?P<plate>\d+)_S(?P<well>\d+)$"
)


def load_coo_triplets(h5_path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(h5_path, "r") as h5:
        i = np.ravel(h5["i"][()])
        j = np.ravel(h5["j"][()])
        v = np.ravel(h5["v"][()])
    return i.astype(np.int64), j.astype(np.int64), v.astype(np.float32)


def parse_week(text: str) -> int | None:
    match = re.search(r"(\d+)\s*w", str(text))
    return int(match.group(1)) if match else None


def first_existing(cols: pd.Index | list[str], keys: list[str], default: str | None = None):
    col_set = set(cols)
    for key in keys:
        if key in col_set:
            return key
    return default


def parse_sample_fields(sample_id: str) -> Dict[str, object]:
    match = SAMPLE_RE.match(str(sample_id))
    if not match:
        return {}
    data = match.groupdict()
    data["week"] = int(data["week"])
    data["timepoint"] = int(data["timepoint"])
    data["plate"] = int(data["plate"])
    data["well"] = int(data["well"])
    data["mouse"] = f"m{data['mouse']}"
    return data


def ingest_timecourse_scrna(
    root: Path | None = None,
    output_path: Path | None = None,
    *,
    config: Dict[str, object] | None = None,
) -> ad.AnnData:
    base = repo_root(root)
    params = config or load_config(base)
    raw_dir = base / params["paths"]["raw_timecourse"]

    gene_path = raw_dir / "GSE154989_mmLungPlate_fQC_geneTable.csv.gz"
    qcstat_path = raw_dir / "GSE154989_mmLungPlate_fQC_dZ_QCstat_smpTable.csv.gz"
    smp_path = raw_dir / "GSE154989_mmLungPlate_fQC_smpTable.csv.gz"
    annot_path = raw_dir / "GSE154989_mmLungPlate_fQC_dZ_annot_smpTable.csv.gz"

    for path in (gene_path, qcstat_path, smp_path):
        if not path.exists():
            raise FileNotFoundError(f"Missing required file: {path}")

    gene = pd.read_csv(gene_path)
    smp = pd.read_csv(qcstat_path)
    smp_sample = pd.read_csv(smp_path)
    annot = pd.read_csv(annot_path) if annot_path.exists() else pd.DataFrame()

    gene_id_col = first_existing(
        gene.columns,
        ["geneSymbol", "gene_symbol", "symbol", "gene_name", "GeneSymbol"],
    )
    sample_id_col = first_existing(smp.columns, ["sampleID", "SampleID", "sample_id"])
    if gene_id_col is None or sample_id_col is None:
        raise ValueError("Unable to infer gene symbol or sample ID columns.")

    cell_id_col = sample_id_col
    if smp[sample_id_col].duplicated().any():
        smp = smp.reset_index(drop=False).rename(columns={"index": "_row"})
        smp["cell_id"] = (
            smp[sample_id_col].astype(str) + "__" + smp["_row"].astype(str).str.zfill(6)
        )
        cell_id_col = "cell_id"
    else:
        smp[sample_id_col] = smp[sample_id_col].astype(str)

    parsed = smp[sample_id_col].map(parse_sample_fields)
    parsed_df = pd.DataFrame(list(parsed))
    for col in ("genotype", "treatment"):
        if col in parsed_df:
            parsed_df[col] = pd.Categorical(parsed_df[col])
    if "week" in parsed_df:
        parsed_df["week"] = parsed_df["week"].astype("Int64")
    smp = pd.concat([smp, parsed_df], axis=1)

    n_cells_expected = len(smp)
    n_genes_expected = len(gene)

    coo_h5 = None
    for cand in [
        raw_dir / "GSE154989_mmLungPlate_fQC_dSp_rawCount.h5",
        raw_dir / "GSE154989_mmLungPlate_fQC_dSp_rawCountOrig.h5",
    ]:
        if cand.exists():
            coo_h5 = cand
            break
    if coo_h5 is None:
        raise FileNotFoundError("Could not find a COO count file (rawCount/rawCountOrig).")

    i, j, v = load_coo_triplets(coo_h5)
    one_based_i = (i.min() == 1) or (i.max() == n_genes_expected)
    one_based_j = (j.min() == 1) or (j.max() == n_cells_expected)
    if one_based_i:
        i -= 1
    if one_based_j:
        j -= 1

    ok_cg = (i.max() < n_cells_expected) and (j.max() < n_genes_expected)
    ok_gc = (i.max() < n_genes_expected) and (j.max() < n_cells_expected)
    if ok_cg and not ok_gc:
        X = sp.coo_matrix((v, (i, j)), shape=(n_cells_expected, n_genes_expected)).tocsr()
    elif ok_gc and not ok_cg:
        X_gc = sp.coo_matrix((v, (i, j)), shape=(n_genes_expected, n_cells_expected))
        X = X_gc.T.tocsr()
    elif ok_gc and ok_cg:
        X_gc = sp.coo_matrix((v, (i, j)), shape=(n_genes_expected, n_cells_expected))
        X = X_gc.T.tocsr()
    else:
        raise ValueError("Indices exceed expected dimensions for count matrix.")

    X = X.astype(np.float32, copy=False)

    merge_key = "sampleID" if "sampleID" in smp.columns else sample_id_col
    keep_cols = [
        col
        for col in ["sampleID", "timesimple", "mouseID", "plateID", "typeID"]
        if col in smp_sample.columns
    ]
    smp_aug = smp.merge(smp_sample[keep_cols], on="sampleID", how="left") if keep_cols else smp.copy()

    if ("week" not in smp_aug.columns) or smp_aug["week"].isna().all():
        if "timesimple" in smp_aug.columns:
            smp_aug["timesimple"] = smp_aug["timesimple"].astype(str)
            smp_aug["week"] = smp_aug["timesimple"].map(parse_week).astype("Int64")
        else:
            smp_aug["week"] = pd.NA

    if "week" in smp_aug.columns:
        late_mask = smp_aug["week"].ge(12).fillna(False)
        stage_vals = np.where(late_mask, "late", "early")
        smp_aug["stage_label"] = pd.Categorical(
            stage_vals, categories=["early", "late"], ordered=True
        )
    else:
        smp_aug["stage_label"] = pd.Categorical(["unknown"] * len(smp_aug))

    order_timesimple = [
        "01_T_early_ND",
        "02_KorKP_early_ND",
        "04_K_12w_ND",
        "06_KP_12w_ND",
        "07_KP_20w_ND",
        "08_KP_30w_ND",
    ]
    if "timesimple" in smp_aug.columns:
        smp_aug["timesimple"] = pd.Categorical(
            smp_aug["timesimple"], categories=order_timesimple, ordered=True
        )

    smp_aug["sampleID"] = smp_aug["sampleID"].astype(str)
    if not annot.empty and "sampleID" in annot.columns:
        annot = annot.copy()
        annot["sampleID"] = annot["sampleID"].astype(str)
        embed_cols = [c for c in ["tSNE_1", "tSNE_2", "phate_1", "phate_2"] if c in annot.columns]
        if embed_cols:
            annot_small = annot[["sampleID"] + embed_cols].drop_duplicates("sampleID")
            smp_aug = smp_aug.merge(annot_small, on="sampleID", how="left")

    adata = ad.AnnData(X)
    adata.obs = smp_aug.copy()
    adata.var = gene.copy()

    adata.obs_names = adata.obs[cell_id_col].astype(str).values
    if adata.obs_names.has_duplicates:
        adata.obs_names_make_unique()
    adata.obs.index = adata.obs.index.astype(str)

    adata.var_names = adata.var[gene_id_col].astype(str).values
    if adata.var_names.has_duplicates:
        adata.var_names_make_unique()
    adata.var.index = adata.var.index.astype(str)

    adata.layers["counts"] = adata.X.copy()

    adata.uns.setdefault("sources", {})
    adata.uns["sources"]["scrna_timecourse_gse154989"] = str(raw_dir)
    adata.uns["sources"]["scrna_tigit_gse154978"] = "data/raw/scrna/tigit_gse154978"
    adata.uns["sources"]["scrna_unsorted_gse154977"] = "data/raw/scrna/tenx_unsorted_gse154977"
    adata.uns["sources"]["scatac_gse154965"] = "data/raw/scatac/gse154965"
    adata.uns["sources"]["bulk_atac_gse154966"] = "data/raw/bulk_atac/gse154966"

    if "week" in adata.obs:
        w = adata.obs["week"].dropna()
        if len(w) and ((w < 0).any() or (w > 60).any()):
            raise ValueError("Week outside expected range.")

    if output_path is None:
        output_path = base / params["paths"]["anndata"]
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if "tSNE_1" in adata.obs and "tSNE_2" in adata.obs:
        tsne = np.vstack(
            [
                pd.to_numeric(adata.obs["tSNE_1"], errors="coerce").to_numpy(),
                pd.to_numeric(adata.obs["tSNE_2"], errors="coerce").to_numpy(),
            ]
        ).T.astype(np.float32)
        adata.obsm["X_tsne"] = tsne
    if "phate_1" in adata.obs and "phate_2" in adata.obs:
        phate = np.vstack(
            [
                pd.to_numeric(adata.obs["phate_1"], errors="coerce").to_numpy(),
                pd.to_numeric(adata.obs["phate_2"], errors="coerce").to_numpy(),
            ]
        ).T.astype(np.float32)
        adata.obsm["X_phate"] = phate

    adata.write(output_path)
    return adata
