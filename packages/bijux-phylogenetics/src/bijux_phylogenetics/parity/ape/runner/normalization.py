from __future__ import annotations

import csv
import json
from pathlib import Path
import re

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import TreeNode

from .tree_payloads import _sort_parity_rows


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _normalize_newick_label(label: str) -> str:
    if len(label) >= 2 and label.startswith("'") and label.endswith("'"):
        return label[1:-1].replace("''", "'")
    return label


def _normalize_expected_label(
    label: str,
    *,
    expected_tip_labels: set[str] | None,
) -> str:
    normalized = _normalize_newick_label(label)
    if (
        expected_tip_labels
        and normalized not in expected_tip_labels
        and normalized.replace("_", " ") in expected_tip_labels
    ):
        return normalized.replace("_", " ")
    return normalized


def _normalize_joined_labels(
    value: str,
    *,
    expected_tip_labels: set[str] | None,
) -> str:
    if value == "":
        return value
    labels = [
        _normalize_expected_label(label, expected_tip_labels=expected_tip_labels)
        for label in value.split("|")
    ]
    return "|".join(sorted(labels))


def _normalize_reference_summary(summary: dict[str, object]) -> dict[str, object]:
    normalized = dict(summary)
    tip_labels = normalized.get("tip_labels")
    if isinstance(tip_labels, list):
        expected_tip_labels = {
            _normalize_newick_label(str(label)) for label in tip_labels
        }
        normalized["tip_labels"] = [
            _normalize_expected_label(
                str(label),
                expected_tip_labels=expected_tip_labels,
            )
            for label in tip_labels
        ]
    return normalized


def _summary_rooted_flag(summary: dict[str, object]) -> bool:
    rooted = summary.get("rooted")
    if isinstance(rooted, bool):
        return rooted
    raise ValueError("reference summary must include a boolean rooted flag")


def _optional_payload_string(payload: dict[str, object], key: str) -> str | None:
    value = payload.get(key)
    return value if isinstance(value, str) else None


def _coerce_table_cell(value: str) -> object:
    if value == "":
        return ""
    if value.lower() in {"true", "false"}:
        return value.lower() == "true"
    if re.fullmatch(r"-?\d+", value):
        return int(value)
    if re.fullmatch(r"-?(?:\d+\.\d*|\d*\.\d+)(?:[eE][+-]?\d+)?", value):
        return float(value)
    return value


def _load_rows_table(
    path: Path,
    *,
    expected_tip_labels: set[str] | None = None,
    sort_rows: bool = False,
) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    normalized_rows = [
        {
            key: _normalize_expected_label(
                value, expected_tip_labels=expected_tip_labels
            )
            if key.endswith("label")
            else _normalize_joined_labels(
                value,
                expected_tip_labels=expected_tip_labels,
            )
            if key in {"clade_id", "taxa", "descendant_taxa", "shared_taxa"}
            else _coerce_table_cell(value)
            for key, value in row.items()
        }
        for row in rows
    ]
    if sort_rows:
        return _sort_parity_rows(normalized_rows)
    return normalized_rows


def _compare_scalar(expected: object, observed: object, *, tolerance: float) -> bool:
    if isinstance(expected, (int, float)) and isinstance(observed, (int, float)):
        return abs(float(expected) - float(observed)) <= tolerance
    return expected == observed


def _compare_json(expected: object, observed: object, *, tolerance: float) -> bool:
    if isinstance(expected, dict) and isinstance(observed, dict):
        if set(expected) != set(observed):
            return False
        return all(
            _compare_json(expected[key], observed[key], tolerance=tolerance)
            for key in expected
        )
    if isinstance(expected, list) and isinstance(observed, list):
        if len(expected) != len(observed):
            return False
        return all(
            _compare_json(left, right, tolerance=tolerance)
            for left, right in zip(expected, observed, strict=True)
        )
    return _compare_scalar(expected, observed, tolerance=tolerance)


def _normalize_tree_labels(
    node: TreeNode,
    *,
    expected_tip_labels: set[str] | None,
) -> None:
    if node.name is not None:
        normalized = _normalize_newick_label(node.name)
        if (
            expected_tip_labels
            and normalized not in expected_tip_labels
            and normalized.replace("_", " ") in expected_tip_labels
        ):
            normalized = normalized.replace("_", " ")
        node.name = normalized
    for child in node.children:
        _normalize_tree_labels(child, expected_tip_labels=expected_tip_labels)


def _clear_branch_lengths(node: TreeNode) -> None:
    node.branch_length = None
    for child in node.children:
        _clear_branch_lengths(child)


def _canonical_newick(
    path: Path,
    *,
    expected_tip_labels: set[str] | None = None,
) -> str:
    tree = load_tree(path)
    _normalize_tree_labels(tree.root, expected_tip_labels=expected_tip_labels)
    return dumps_newick(tree)
