from __future__ import annotations

import csv
import json
from pathlib import Path

_CONFIG_FILENAME = "workflow-config.resolved.json"
_SLURM_STORAGE_CATEGORIES_FILENAME = "slurm-storage-categories.tsv"
_SLURM_OUTPUT_EXPLOSION_SUMMARY_FILENAME = "slurm-output-explosion-report.json"
_MEBIBYTE = 1024 * 1024
_TREE_FILE_SUFFIXES = {
    ".boottrees",
    ".contree",
    ".newick",
    ".nwk",
    ".phy",
    ".tree",
    ".treefile",
    ".trees",
    ".tre",
}
_COMPRESSED_SUFFIXES = {".bz2", ".gz", ".xz"}
_THINNING_RECOMMENDED_TREE_COUNT = 2_000
_THINNING_REQUIRED_TREE_COUNT = 10_000
_THINNING_TARGET_TREE_COUNT = 1_000
_COMPRESSION_RECOMMENDED_TREE_COUNT = 2_000
_COMPRESSION_REQUIRED_TREE_COUNT = 10_000
_COMPRESSION_RECOMMENDED_BYTES = 4 * _MEBIBYTE
_COMPRESSION_REQUIRED_BYTES = 64 * _MEBIBYTE


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def _write_tsv(
    path: Path,
    *,
    fieldnames: tuple[str, ...],
    rows: list[dict[str, object]],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path
