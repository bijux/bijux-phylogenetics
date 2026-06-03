from __future__ import annotations

import csv
import json
from pathlib import Path

_EXECUTION_RECORD_FILENAME = "rabies-method-sensitivity-panel.run.json"
_CONFIG_FILENAME = "workflow-config.resolved.json"
_SLURM_ARRAY_PARTITIONS_FILENAME = "slurm-array-partitions.tsv"
_SLURM_ARRAY_MEMBERS_FILENAME = "slurm-array-members.tsv"


def collect_output_observations(output_root: Path) -> tuple[int, int]:
    if not output_root.exists():
        return 0, 0
    file_paths = sorted(path for path in output_root.rglob("*") if path.is_file())
    return len(file_paths), sum(path.stat().st_size for path in file_paths)


def missing_required_variant_outputs(
    output_root: Path,
    variant_id: str,
) -> tuple[str, ...]:
    return tuple(
        filename
        for filename in (
            f"{variant_id}.aln",
            f"{variant_id}.trimmed.aln",
            "fasttree.nwk",
            "iqtree-support.nwk",
            "rooted-engine-comparison.tsv",
            "rooted-fasttree.nwk",
            "rooted-iqtree-support.nwk",
            "rooting-summary.tsv",
            "unrooted-comparison.tsv",
            "unrooted-conclusions.tsv",
            "unrooted-conflicting-clades.tsv",
            "unrooted-shared-clades.tsv",
            "unrooted-stability-summary.tsv",
            "unrooted-support-weighted-conflicts.tsv",
        )
        if not (output_root / filename).is_file()
    )


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def relative_bundle_path(bundle_root: Path, value: Path) -> str:
    try:
        return value.relative_to(bundle_root).as_posix()
    except ValueError:
        return value.as_posix()


def parse_task_log(path: Path) -> dict[str, str]:
    payload: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        payload[key.strip()] = value.strip()
    return payload


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
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
            writer.writerow(
                {key: "" if value is None else value for key, value in row.items()}
            )
    return path
