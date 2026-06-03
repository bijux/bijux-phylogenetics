from __future__ import annotations

import csv
import hashlib
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.diagnostics.validation.structure import _load_tree
from bijux_phylogenetics.io.fasta.records import summarise_fasta

from .models import ReportInputLedgerEntry


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dataset_surface_taxa_count(path: Path, role: str) -> int:
    if role in {"tree", "inference_tree"}:
        return _load_tree(path).tip_count
    if role in {"alignment", "filtered_alignment"}:
        return summarise_fasta(path).sequence_count
    if role in {"metadata", "traits", "tip_dates", "reported_taxa"}:
        return load_taxon_table(path).row_count
    if role in {"synonym_table", "calibrations"}:
        delimiter = "," if path.suffix.lower() == ".csv" else "\t"
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            return sum(1 for _ in reader)
    raise ValueError(f"unsupported dataset ledger role: {role}")


def build_input_ledger(
    entries: list[tuple[Path, str, list[str]]],
) -> list[ReportInputLedgerEntry]:
    return [
        ReportInputLedgerEntry(
            path=path,
            role=role,
            checksum=sha256(path),
            taxa_count=dataset_surface_taxa_count(path, role),
            usage=usage,
        )
        for path, role, usage in entries
    ]


def serialize_input_ledger(
    entries: list[ReportInputLedgerEntry],
) -> list[dict[str, object]]:
    return [
        {
            "path": str(entry.path),
            "role": entry.role,
            "checksum": entry.checksum,
            "taxa_count": entry.taxa_count,
            "usage": entry.usage,
        }
        for entry in entries
    ]
