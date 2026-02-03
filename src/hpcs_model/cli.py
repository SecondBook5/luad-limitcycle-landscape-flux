from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch

from hpcs_model.config import load_config, output_paths
from hpcs_model.io.loaders import load_anndata
from hpcs_model.io.validate_inputs import validate_adata
from hpcs_model.preprocess import batch_correct, normalize_log1p, qc_filter, run_pca, select_hvg
from hpcs_model.clustering.cluster import run_leiden
from hpcs_model.clustering.markers import rank_markers
from hpcs_model.dynamics.ot_coupling import compute_couplings
from hpcs_model.dynamics.aggregate_graph import aggregate_to_cluster_graph
from hpcs_model.dynamics.plasticity import bootstrap_hpcs, compute_plasticity, select_hpcs
from hpcs_model.landscape_flux import (
    compute_potential,
    dominant_paths,
    drift_from_couplings,
    estimate_density,
    flux_residual,
)
from hpcs_model.ml.train import train_drift_model
from hpcs_model.ml.eval import evaluate_drift_model
from hpcs_model.viz.export import ensure_dirs
from hpcs_model.viz.plots import plot_drift, plot_embedding, plot_flux_magnitude, plot_potential


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="hpcs-model")
    parser.add_argument("--config", type=Path, required=True)
    sub = parser.add_subparsers(dest="command", required=True)
    for cmd in ["validate", "preprocess", "couple", "graph", "landscape_flux", "train_ml", "run"]:
        sub.add_parser(cmd)
    return parser.parse_args()


def cmd_validate(cfg):
    adata = load_anndata(cfg)
    validate_adata(adata)


def cmd_preprocess(cfg):
    adata = load_anndata(cfg)
    qc_filter(adata, cfg)
    normalize_log1p(adata)
    select_hvg(adata, cfg)
    run_pca(adata, cfg)
    batch_correct(adata, cfg)
    outputs = output_paths(cfg)
    ensure_dirs(outputs["root"], outputs["figures"], outputs["tables"])
    adata.write(outputs["root"] / "adata_preprocessed.h5ad")
    Z2 = adata.obsm["X_umap"] if "X_umap" in adata.obsm else adata.obsm["X_pca"][:, :2]
    plot_embedding(Z2, adata.obs["timepoint"].cat.codes, outputs["figures"] / "embedding_timepoint.png", "timepoint")


def cmd_couple(cfg):
    outputs = output_paths(cfg)
    adata = ad.read_h5ad(outputs["root"] / "adata_preprocessed.h5ad")
    coupling_dir = outputs["tables"] / "couplings"
    summaries = compute_couplings(adata, cfg, coupling_dir)
    with (outputs["tables"] / "coupling_summary.json").open("w") as fh:
        json.dump(summaries, fh, indent=2)


def cmd_graph(cfg):
    outputs = output_paths(cfg)
    adata = ad.read_h5ad(outputs["root"] / "adata_preprocessed.h5ad")
    coupling_dir = outputs["tables"] / "couplings"
    summary_path = outputs["tables"] / "coupling_summary.json"
    summaries = json.loads(summary_path.read_text()) if summary_path.exists() else []
    edge_df = aggregate_to_cluster_graph(adata, cfg, summaries, coupling_dir)
    plasticity_df = compute_plasticity(edge_df)
    if not plasticity_df.empty:
        plasticity_df.to_csv(outputs["tables"] / "cluster_plasticity.csv", index=False)
        hpcs_cluster = select_hpcs(plasticity_df)
        bootstrap_df = bootstrap_hpcs(edge_df, cfg)
        bootstrap_df.to_csv(outputs["tables"] / "hpcs_bootstrap.csv", index=False)
        with (outputs["tables"] / "hpcs_selection.json").open("w") as fh:
            json.dump({"hpcs_cluster": hpcs_cluster}, fh, indent=2)

    if cfg["clustering"]["cluster_key"] in adata.obs:
        Z2 = adata.obsm["X_umap"] if "X_umap" in adata.obsm else adata.obsm["X_pca"][:, :2]
        plot_embedding(Z2, adata.obs[cfg["clustering"]["cluster_key"]], outputs["figures"] / "embedding_clusters.png", "clusters")

    markers = rank_markers(adata, cfg)
    markers.to_csv(outputs["tables"] / "cluster_markers.csv", index=False)


def cmd_landscape_flux(cfg):
    outputs = output_paths(cfg)
    adata = ad.read_h5ad(outputs["root"] / "adata_preprocessed.h5ad")
    Z2 = adata.obsm["X_umap"] if "X_umap" in adata.obsm else adata.obsm["X_pca"][:, :2]
    density, radius = estimate_density(Z2, cfg)
    U = compute_potential(density, cfg)
    grid_size = int(cfg["landscape_flux"]["grid_size"])
    grid_x = np.linspace(Z2[:, 0].min(), Z2[:, 0].max(), grid_size)
    grid_y = np.linspace(Z2[:, 1].min(), Z2[:, 1].max(), grid_size)
    grid = np.stack(np.meshgrid(grid_x, grid_y, indexing="ij"), axis=-1).reshape(-1, 2)
    grid_density, _ = estimate_density(grid, cfg)
    U_grid = compute_potential(grid_density, cfg).reshape(grid_size, grid_size)

    coupling_dir = outputs["tables"] / "couplings"
    summaries = json.loads((outputs["tables"] / "coupling_summary.json").read_text())
    time_pairs = [(s["time_t"], s["time_t1"]) for s in summaries]
    files = [coupling_dir / s["file"] for s in summaries]
    drift = drift_from_couplings(Z2, adata.obs["timepoint"].to_numpy(), time_pairs, files)
    residual, mag, curl = flux_residual(Z2, U_grid, grid_x, grid_y, drift)

    flux_df = pd.DataFrame(
        {
            "cell_id": adata.obs_names,
            "drift_x": drift[:, 0],
            "drift_y": drift[:, 1],
            "flux_mag": mag,
        }
    )
    flux_df.to_csv(outputs["tables"] / "flux_residuals.csv", index=False)
    plot_potential(grid_x, grid_y, U_grid, outputs["figures"] / "potential_contours.png")
    plot_drift(Z2, drift, outputs["figures"] / "drift_quiver.png")
    plot_flux_magnitude(Z2, mag, outputs["figures"] / "flux_residual_magnitude.png")

    paths = dominant_paths(Z2, U, k=3)
    with (outputs["tables"] / "dominant_paths.json").open("w") as fh:
        json.dump(paths, fh, indent=2)


def cmd_train_ml(cfg):
    outputs = output_paths(cfg)
    adata = ad.read_h5ad(outputs["root"] / "adata_preprocessed.h5ad")
    Z2 = adata.obsm["X_pca"]
    summaries = json.loads((outputs["tables"] / "coupling_summary.json").read_text())
    time_pairs = [(s["time_t"], s["time_t1"]) for s in summaries]
    coupling_dir = outputs["tables"] / "couplings"
    files = [coupling_dir / s["file"] for s in summaries]
    drift = drift_from_couplings(Z2, adata.obs["timepoint"].to_numpy(), time_pairs, files)
    condition = adata.obs["timepoint"].cat.codes.to_numpy()
    model_path = train_drift_model(Z2, drift, condition, cfg, outputs["models"])
    model_state = torch.load(model_path, map_location="cpu")
    mse = evaluate_drift_model(Z2, drift, condition, cfg, model_state)
    with (outputs["tables"] / "ml_eval.json").open("w") as fh:
        json.dump({"drift_mse": mse}, fh, indent=2)


def cmd_run(cfg):
    cmd_validate(cfg)
    cmd_preprocess(cfg)
    cmd_couple(cfg)
    cmd_graph(cfg)
    cmd_landscape_flux(cfg)
    cmd_train_ml(cfg)


def main() -> None:
    args = _parse_args()
    cfg = load_config(args.config)
    command = args.command
    if command == "validate":
        cmd_validate(cfg)
    elif command == "preprocess":
        cmd_preprocess(cfg)
    elif command == "couple":
        cmd_couple(cfg)
    elif command == "graph":
        cmd_graph(cfg)
    elif command == "landscape_flux":
        cmd_landscape_flux(cfg)
    elif command == "train_ml":
        cmd_train_ml(cfg)
    elif command == "run":
        cmd_run(cfg)


if __name__ == "__main__":
    main()
