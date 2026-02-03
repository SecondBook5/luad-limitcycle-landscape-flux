from __future__ import annotations

from pathlib import Path
from typing import Iterable, List


def load_markers(path: Path | str) -> List[str]:
    marker_path = Path(path)
    if not marker_path.exists():
        raise FileNotFoundError(f"Marker file not found: {marker_path}")
    return [
        line.strip()
        for line in marker_path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def filter_present(markers: Iterable[str], genes: Iterable[str]) -> list[str]:
    gene_set = set(map(str, genes))
    return [gene for gene in markers if str(gene) in gene_set]
