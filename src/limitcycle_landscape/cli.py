from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad

from limitcycle_landscape.config import load_config, repo_root
from limitcycle_landscape.cycle import bin_atac_by_theta, fit_cycle_from_markers
from limitcycle_landscape.ingest import ingest_timecourse_scrna


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="limitcycle-landscape",
        description="Limit-cycle landscape CLI for cell-cycle modeling.",
    )
    parser.add_argument("--root", type=Path, default=None, help="Repo root (defaults to auto-detect).")

    subparsers = parser.add_subparsers(dest="command", required=True)

    ingest_parser = subparsers.add_parser("ingest-scrna", help="Ingest scRNA timecourse into AnnData.")
    ingest_parser.add_argument(
        "--output", type=Path, default=None, help="Override output .h5ad path."
    )

    cycle_parser = subparsers.add_parser("fit-cycle", help="Fit cycle phase using marker genes.")
    cycle_parser.add_argument(
        "--output", type=Path, default=None, help="Override output .h5ad path."
    )
    cycle_parser.add_argument(
        "--top-frac", type=float, default=0.05, help="Top fraction of S/G2M anchors."
    )

    atac_parser = subparsers.add_parser(
        "bin-atac", help="Aggregate ATAC accessibility by theta bins."
    )
    atac_parser.add_argument(
        "--layer", type=str, default="X_atac", help="ATAC layer name in AnnData."
    )
    atac_parser.add_argument(
        "--bins", type=int, default=60, help="Number of theta bins."
    )
    atac_parser.add_argument(
        "--input", type=Path, default=None, help="AnnData input (defaults to config path)."
    )

    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    base = repo_root(args.root)
    config = load_config(base)

    if args.command == "ingest-scrna":
        ingest_timecourse_scrna(base, output_path=args.output, config=config)
        return

    if args.command == "fit-cycle":
        fit_cycle_from_markers(
            base, config=config, output_path=args.output, top_frac=args.top_frac
        )
        return

    if args.command == "bin-atac":
        input_path = args.input or (base / config["paths"]["anndata"])
        adata = ad.read_h5ad(input_path)
        payload = bin_atac_by_theta(adata, atac_layer=args.layer, bins=args.bins)
        adata.write(input_path)
        print(
            f"[saved] atac_pseudobulk bins={len(payload['counts'])} "
            f"layer={args.layer} → {input_path}"
        )


if __name__ == "__main__":
    main()
