from __future__ import annotations

from pathlib import Path

from ..registry import ApeParityCase
from .normalization import (
    _canonical_newick,
    _load_json,
    _load_rows_table,
    _normalize_reference_summary,
)


def _load_reference_case_payload(
    case: ApeParityCase,
    execution_root: Path,
) -> tuple[dict[str, object], list[dict[str, object]] | None, str | None]:
    if case.operation in {
        "read-tree-structure",
        "write-tree-structure",
        "root-tree-outgroup",
        "unroot-tree",
        "drop-tree-taxa",
        "keep-tree-taxa",
        "extract-tree-clade",
    }:
        summary = _normalize_reference_summary(
            _load_json(execution_root / "summary.json")
        )
        expected_tip_labels = {str(label) for label in summary.get("tip_labels", [])}
        rows = _load_rows_table(
            execution_root / "clades.tsv",
            expected_tip_labels=expected_tip_labels,
            sort_rows=True,
        )
        normalized_text = _canonical_newick(
            execution_root / "normalized-tree.nwk",
            expected_tip_labels=expected_tip_labels,
        )
        return summary, rows, normalized_text
    if case.operation == "get-tree-mrca":
        summary = _normalize_reference_summary(
            _load_json(execution_root / "summary.json")
        )
        return summary, None, None
    if case.operation == "assess-tree-monophyly":
        summary = _normalize_reference_summary(
            _load_json(execution_root / "summary.json")
        )
        return summary, None, None
    if case.operation in {"read-tree-set-structure", "write-tree-set-structure"}:
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "clades.tsv", sort_rows=True)
        return summary, rows, None
    if case.operation == "tree-consensus":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "clade-frequencies.tsv")
        normalized_text = _canonical_newick(execution_root / "normalized-tree.nwk")
        return summary, rows, normalized_text
    if case.operation == "tree-clade-support":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "support-table.tsv")
        return summary, rows, None
    if case.operation == "tree-tip-distance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "tip-distance-long.tsv")
        return summary, rows, None
    if case.operation == "distance-matrix-neighbor-joining":
        summary = _normalize_reference_summary(
            _load_json(execution_root / "summary.json")
        )
        normalized_text = _canonical_newick(execution_root / "normalized-tree.nwk")
        return summary, None, normalized_text
    if case.operation == "tree-topology-distance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "split-table.tsv")
        return summary, rows, None
    if case.operation == "tree-brownian-covariance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "covariance-long.tsv")
        return summary, rows, None
    if case.operation == "tree-continuous-ancestral-states":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "continuous-ancestral.tsv")
        return summary, rows, None
    if case.operation == "tree-discrete-ancestral-states":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "discrete-ancestral.tsv")
        return summary, rows, None
    if case.operation == "tree-independent-contrasts":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "independent-contrasts.tsv")
        return summary, rows, None
    if case.operation == "tree-node-depth":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "node-depths.tsv")
        return summary, rows, None
    if case.operation == "tree-branching-times":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "branching-times.tsv")
        return summary, rows, None
    if case.operation == "tree-diversification-gamma-statistic":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "gamma-statistic.tsv")
        return summary, rows, None
    if case.operation == "tree-simulation-envelope":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "simulation-envelope.tsv")
        return summary, rows, None
    if case.operation == "tree-ultrametricity":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "ultrametric-diagnostics.tsv")
        return summary, rows, None
    if case.operation == "dna-dnabin-structure":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "dnabin.tsv")
        return summary, rows, None
    if case.operation == "dna-base-frequency":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "base-frequency.tsv")
        return summary, rows, None
    if case.operation == "dna-segregating-sites":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "segregating-sites.tsv")
        return summary, rows, None
    if case.operation == "dna-distance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "distance-matrix.tsv")
        return summary, rows, None
    if case.operation == "dna-translation":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "translation.tsv")
        return summary, rows, None
    raise ValueError(f"unsupported ape parity operation '{case.operation}'")
