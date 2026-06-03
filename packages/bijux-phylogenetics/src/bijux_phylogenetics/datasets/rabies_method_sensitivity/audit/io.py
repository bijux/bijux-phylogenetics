from __future__ import annotations

import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.engines.common import file_sha256
from bijux_phylogenetics.io.fasta import load_fasta_alignment

from .contracts import RabiesMethodSensitivityReproducibilityAuditReport


def write_rabies_method_sensitivity_reproducibility_checks_table(
    path: Path, report: RabiesMethodSensitivityReproducibilityAuditReport
) -> Path:
    """Write one tabular ledger of top-level audit checks."""
    return _write_tsv(
        path,
        fieldnames=("check_id", "surface", "status", "expected", "observed", "detail"),
        rows=[asdict(row) for row in report.checks],
    )


def write_rabies_method_sensitivity_variant_audit_table(
    path: Path, report: RabiesMethodSensitivityReproducibilityAuditReport
) -> Path:
    """Write one per-variant audit ledger."""
    return _write_tsv(
        path,
        fieldnames=(
            "variant_id",
            "status",
            "output_file_count",
            "output_byte_count",
            "output_digest",
            "missing_required_files",
            "unexpected_files",
            "issues",
        ),
        rows=[
            {
                **asdict(row),
                "missing_required_files": "; ".join(row.missing_required_files),
                "unexpected_files": "; ".join(row.unexpected_files),
                "issues": "; ".join(row.issues),
            }
            for row in report.variants
        ],
    )


def write_rabies_method_sensitivity_reproducibility_audit_json(
    path: Path, report: RabiesMethodSensitivityReproducibilityAuditReport
) -> Path:
    """Write one machine-readable JSON summary for the bundle audit."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    return path


def _alignment_length(path: Path) -> int:
    records = load_fasta_alignment(path)
    return len(records[0].sequence)


def _directory_digest(path: Path) -> str:
    if not path.is_dir():
        return ""
    lines = [
        f"{entry.relative_to(path).as_posix()}\t{file_sha256(entry)}"
        for entry in sorted(path.rglob("*"))
        if entry.is_file()
    ]
    payload = "\n".join(lines).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _format_float(value: float) -> str:
    return format(value, ".12g")


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _parse_task_log(path: Path) -> dict[str, str]:
    payload: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        payload[key.strip()] = value.strip()
    return payload


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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
            writer.writerow(
                {key: "" if value is None else value for key, value in row.items()}
            )
    return path
