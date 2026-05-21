from __future__ import annotations

import hashlib
from pathlib import Path
import re

from bijux_phylogenetics.phylo.alignment.partitions import LocusPartition


def _sidecar(path: Path, label: str) -> Path:
    return path.parent / f"{path.name}.{label}"


def _prefix_path(out_dir: Path, prefix: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / prefix


def _manifest_path_from_output(path: Path) -> Path:
    return _sidecar(path, "manifest.json")


def _partition_support_path(prefix_path: Path, suffix: str) -> Path:
    return prefix_path.parent / f"{prefix_path.name}.{suffix}"


def _partition_alignment_file_name(partition: LocusPartition) -> str:
    normalized = re.sub(r"[^A-Za-z0-9._-]+", "-", partition.name.strip().lower())
    normalized = normalized.strip("-._") or "partition"
    digest = hashlib.sha256(partition.name.encode("utf-8")).hexdigest()[:8]
    return f"{normalized}-{digest}.fasta"
