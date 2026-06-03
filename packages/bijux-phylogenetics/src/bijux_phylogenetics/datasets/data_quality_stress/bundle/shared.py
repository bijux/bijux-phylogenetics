from __future__ import annotations

from pathlib import Path
import shutil

from bijux_phylogenetics.diagnostics.validation import TreeValidationReport


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _substantive_alignment_warnings(warnings: list[str]) -> list[str]:
    ignored = {
        "automatic sequence type defaults to dna from nucleotide-like characters that remain protein-compatible by alphabet alone"
    }
    return [warning for warning in warnings if warning not in ignored]


def _tree_warning_nodes(
    report: TreeValidationReport,
    *,
    warning_code: str,
) -> list[str]:
    affected_nodes = [
        node
        for warning in report.warning_details
        if warning.code == warning_code
        for node in warning.affected_nodes
    ]
    return sorted(dict.fromkeys(affected_nodes))


def _format_number(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".12g")
