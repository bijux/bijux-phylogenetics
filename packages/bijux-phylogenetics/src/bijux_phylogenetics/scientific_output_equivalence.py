from __future__ import annotations

import csv
from dataclasses import dataclass
from html.parser import HTMLParser
import json
import math
from pathlib import Path
from typing import Any

from Bio import SeqIO

from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_support_values,
    compare_tree_paths,
)

_TREE_SUFFIXES = {".newick", ".nwk", ".tree"}
_FASTA_SUFFIXES = {".aln", ".fa", ".faa", ".fasta", ".fna"}
_TABULAR_SUFFIXES = {".csv", ".tsv"}


@dataclass(slots=True)
class ScientificOutputEquivalenceIssue:
    """One scientific-output mismatch between expected and observed artifacts."""

    relative_path: Path
    kind: str
    detail: str


@dataclass(slots=True)
class ScientificOutputEquivalenceReport:
    """Directory- or file-level scientific-equivalence comparison result."""

    expected_root: Path
    observed_root: Path
    equivalent: bool
    compared_file_count: int
    issues: list[ScientificOutputEquivalenceIssue]


def compare_scientific_output(
    expected_root: Path,
    observed_root: Path,
    *,
    numeric_tolerance: float = 1e-9,
    branch_length_tolerance: float = 1e-9,
) -> ScientificOutputEquivalenceReport:
    """Compare scientific outputs semantically instead of byte-for-byte."""
    issues: list[ScientificOutputEquivalenceIssue] = []
    compared_file_count = 0

    if expected_root.is_file() and observed_root.is_file():
        compared_file_count = 1
        _compare_file(
            expected_path=expected_root,
            observed_path=observed_root,
            relative_path=Path(expected_root.name),
            issues=issues,
            numeric_tolerance=numeric_tolerance,
            branch_length_tolerance=branch_length_tolerance,
        )
        return ScientificOutputEquivalenceReport(
            expected_root=expected_root,
            observed_root=observed_root,
            equivalent=not issues,
            compared_file_count=compared_file_count,
            issues=issues,
        )

    expected_files = _relative_file_map(expected_root)
    observed_files = _relative_file_map(observed_root)
    expected_relatives = set(expected_files)
    observed_relatives = set(observed_files)

    for relative_path in sorted(expected_relatives - observed_relatives):
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="missing_output",
                detail="expected scientific output is missing from the observed bundle",
            )
        )
    for relative_path in sorted(observed_relatives - expected_relatives):
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="unexpected_output",
                detail="observed bundle contains an extra output not present in the governed reference bundle",
            )
        )

    shared_paths = sorted(expected_relatives & observed_relatives)
    compared_file_count = len(shared_paths)
    for relative_path in shared_paths:
        _compare_file(
            expected_path=expected_files[relative_path],
            observed_path=observed_files[relative_path],
            relative_path=relative_path,
            issues=issues,
            numeric_tolerance=numeric_tolerance,
            branch_length_tolerance=branch_length_tolerance,
        )

    return ScientificOutputEquivalenceReport(
        expected_root=expected_root,
        observed_root=observed_root,
        equivalent=not issues,
        compared_file_count=compared_file_count,
        issues=issues,
    )


def _compare_file(
    *,
    expected_path: Path,
    observed_path: Path,
    relative_path: Path,
    issues: list[ScientificOutputEquivalenceIssue],
    numeric_tolerance: float,
    branch_length_tolerance: float,
) -> None:
    suffix = expected_path.suffix.lower()
    if suffix in _TREE_SUFFIXES:
        _compare_tree_file(
            expected_path=expected_path,
            observed_path=observed_path,
            relative_path=relative_path,
            issues=issues,
            branch_length_tolerance=branch_length_tolerance,
        )
        return
    if suffix in _TABULAR_SUFFIXES:
        _compare_tabular_file(
            expected_path=expected_path,
            observed_path=observed_path,
            relative_path=relative_path,
            issues=issues,
            numeric_tolerance=numeric_tolerance,
        )
        return
    if suffix in _FASTA_SUFFIXES:
        _compare_fasta_file(
            expected_path=expected_path,
            observed_path=observed_path,
            relative_path=relative_path,
            issues=issues,
        )
        return
    if suffix == ".json":
        _compare_json_file(
            expected_path=expected_path,
            observed_path=observed_path,
            relative_path=relative_path,
            issues=issues,
        )
        return
    if suffix == ".html":
        _compare_html_file(
            expected_path=expected_path,
            observed_path=observed_path,
            relative_path=relative_path,
            issues=issues,
        )
        return
    _compare_text_or_binary_file(
        expected_path=expected_path,
        observed_path=observed_path,
        relative_path=relative_path,
        issues=issues,
    )


def _relative_file_map(root: Path) -> dict[Path, Path]:
    return {path.relative_to(root): path for path in root.rglob("*") if path.is_file()}


def _compare_tree_file(
    *,
    expected_path: Path,
    observed_path: Path,
    relative_path: Path,
    issues: list[ScientificOutputEquivalenceIssue],
    branch_length_tolerance: float,
) -> None:
    comparison = compare_tree_paths(
        expected_path,
        observed_path,
        rf_mode="rooted",
        taxon_overlap_policy="require-identical",
    )
    if not comparison.topology_equal:
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="tree_topology_mismatch",
                detail=(
                    "rooted topology changed relative to the governed output "
                    f"(rooted_rf={comparison.rooted_robinson_foulds_distance})"
                ),
            )
        )

    branch_lengths = compare_branch_lengths(
        expected_path,
        observed_path,
        taxon_overlap_policy="require-identical",
    )
    for split in branch_lengths.shared_splits:
        if split.left_length is None and split.right_length is None:
            continue
        if split.left_length is None or split.right_length is None:
            issues.append(
                ScientificOutputEquivalenceIssue(
                    relative_path=relative_path,
                    kind="branch_length_presence_mismatch",
                    detail=f"branch length presence changed for clade {split.split_id}",
                )
            )
            continue
        if not math.isclose(
            split.left_length,
            split.right_length,
            rel_tol=branch_length_tolerance,
            abs_tol=branch_length_tolerance,
        ):
            issues.append(
                ScientificOutputEquivalenceIssue(
                    relative_path=relative_path,
                    kind="branch_length_delta_exceeded",
                    detail=(
                        f"branch length changed beyond tolerance for clade {split.split_id} "
                        f"({split.left_length} -> {split.right_length})"
                    ),
                )
            )

    support = compare_support_values(expected_path, observed_path)
    if support.conflicting_clades:
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="support_value_mismatch",
                detail=(
                    "support-by-clade changed relative to the governed output "
                    f"({len(support.conflicting_clades)} conflicting clades)"
                ),
            )
        )


def _compare_tabular_file(
    *,
    expected_path: Path,
    observed_path: Path,
    relative_path: Path,
    issues: list[ScientificOutputEquivalenceIssue],
    numeric_tolerance: float,
) -> None:
    delimiter = "\t" if expected_path.suffix.lower() == ".tsv" else ","
    expected_rows, expected_fields = _read_tabular_rows(
        expected_path, delimiter=delimiter
    )
    observed_rows, observed_fields = _read_tabular_rows(
        observed_path, delimiter=delimiter
    )
    if expected_fields != observed_fields:
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="table_schema_mismatch",
                detail=(
                    f"table columns changed from {expected_fields!r} to {observed_fields!r}"
                ),
            )
        )
        return

    if len(expected_rows) != len(observed_rows):
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="table_row_count_mismatch",
                detail=(
                    f"table row count changed from {len(expected_rows)} to {len(observed_rows)}"
                ),
            )
        )
        return

    identity_fields = _resolve_tabular_identity_fields(
        expected_fields,
        expected_rows,
        observed_rows,
    )
    sort_fields = identity_fields or expected_fields
    expected_sorted = sorted(
        expected_rows, key=lambda row: _row_sort_key(row, sort_fields)
    )
    observed_sorted = sorted(
        observed_rows, key=lambda row: _row_sort_key(row, sort_fields)
    )
    for index, (expected_row, observed_row) in enumerate(
        zip(expected_sorted, observed_sorted, strict=True),
        start=1,
    ):
        for field in expected_fields:
            expected_value = expected_row[field]
            observed_value = observed_row[field]
            field_tolerance = _table_field_tolerance(
                field,
                relative_path=relative_path,
                base_tolerance=numeric_tolerance,
            )
            if _cells_match(
                expected_value,
                observed_value,
                numeric_tolerance=field_tolerance,
            ):
                continue
            if _field_is_tie_equivalent(
                field,
                expected_row,
                observed_row,
                numeric_tolerance=field_tolerance,
            ):
                continue
            issues.append(
                ScientificOutputEquivalenceIssue(
                    relative_path=relative_path,
                    kind="table_value_mismatch",
                    detail=(
                        f"row {index} field {field!r} changed from "
                        f"{expected_value!r} to {observed_value!r}"
                    ),
                )
            )


def _compare_fasta_file(
    *,
    expected_path: Path,
    observed_path: Path,
    relative_path: Path,
    issues: list[ScientificOutputEquivalenceIssue],
) -> None:
    expected_records = _load_fasta_records(expected_path)
    observed_records = _load_fasta_records(observed_path)
    if expected_records != observed_records:
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="alignment_sequence_mismatch",
                detail="alignment or FASTA sequence content changed relative to the governed output",
            )
        )


def _compare_json_file(
    *,
    expected_path: Path,
    observed_path: Path,
    relative_path: Path,
    issues: list[ScientificOutputEquivalenceIssue],
) -> None:
    expected_payload = json.loads(expected_path.read_text(encoding="utf-8"))
    observed_payload = json.loads(observed_path.read_text(encoding="utf-8"))
    if expected_payload != observed_payload:
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="json_payload_mismatch",
                detail="structured JSON content changed relative to the governed output",
            )
        )


def _compare_html_file(
    *,
    expected_path: Path,
    observed_path: Path,
    relative_path: Path,
    issues: list[ScientificOutputEquivalenceIssue],
) -> None:
    expected_report = _parse_html_contract(expected_path)
    observed_report = _parse_html_contract(observed_path)
    if expected_report.title != observed_report.title:
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="report_title_mismatch",
                detail=(
                    f"report title changed from {expected_report.title!r} "
                    f"to {observed_report.title!r}"
                ),
            )
        )
    if expected_report.headings != observed_report.headings:
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="report_heading_mismatch",
                detail="report heading contract changed relative to the governed output",
            )
        )
    if expected_report.local_links != observed_report.local_links:
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="report_artifact_link_mismatch",
                detail="report linked-artifact contract changed relative to the governed output",
            )
        )
    if _normalize_inline_report_manifest(
        expected_report.inline_manifest
    ) != _normalize_inline_report_manifest(observed_report.inline_manifest):
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="report_manifest_mismatch",
                detail="embedded report manifest content changed relative to the governed output",
            )
        )


def _compare_text_or_binary_file(
    *,
    expected_path: Path,
    observed_path: Path,
    relative_path: Path,
    issues: list[ScientificOutputEquivalenceIssue],
) -> None:
    if expected_path.read_bytes() != observed_path.read_bytes():
        issues.append(
            ScientificOutputEquivalenceIssue(
                relative_path=relative_path,
                kind="artifact_bytes_mismatch",
                detail="artifact bytes changed and no scientific-equivalence comparator is registered for this file type",
            )
        )


def _read_tabular_rows(
    path: Path, *, delimiter: str
) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fieldnames = list(reader.fieldnames or [])
        rows = [{field: row.get(field, "") for field in fieldnames} for row in reader]
    return rows, fieldnames


def _row_sort_key(row: dict[str, str], fields: list[str]) -> tuple[str, ...]:
    return tuple(row[field].strip() for field in fields)


def _table_field_tolerance(
    field: str,
    *,
    relative_path: Path,
    base_tolerance: float,
) -> float:
    if field == "confidence" or field.endswith(("_confidence", "_probabilities")):
        return max(base_tolerance, 0.05)
    if relative_path.name.endswith("-probabilities.tsv"):
        return max(base_tolerance, 0.05)
    return base_tolerance


def _resolve_tabular_identity_fields(
    fields: list[str],
    expected_rows: list[dict[str, str]],
    observed_rows: list[dict[str, str]],
) -> list[str] | None:
    preferred_fields = [
        "node",
        "node_id",
        "clade",
        "split_id",
        "taxon",
        "species",
        "variant_id",
        "label",
        "name",
        "tree_side",
        "source_label",
        "target_label",
    ]
    for field in preferred_fields:
        if field not in fields:
            continue
        expected_values = [row[field].strip() for row in expected_rows]
        observed_values = [row[field].strip() for row in observed_rows]
        if len(set(expected_values)) != len(expected_values):
            continue
        if len(set(observed_values)) != len(observed_values):
            continue
        if set(expected_values) != set(observed_values):
            continue
        return [field]
    return None


def _cells_match(
    expected_value: str, observed_value: str, *, numeric_tolerance: float
) -> bool:
    expected_text = expected_value.strip()
    observed_text = observed_value.strip()
    if expected_text == observed_text:
        return True
    if _structured_cell_match(
        expected_text,
        observed_text,
        numeric_tolerance=numeric_tolerance,
    ):
        return True
    expected_numeric = _maybe_float(expected_text)
    observed_numeric = _maybe_float(observed_text)
    if expected_numeric is None or observed_numeric is None:
        return False
    return math.isclose(
        expected_numeric,
        observed_numeric,
        rel_tol=numeric_tolerance,
        abs_tol=numeric_tolerance,
    )


def _maybe_float(value: str) -> float | None:
    if not value:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    if math.isnan(parsed) or math.isinf(parsed):
        return None
    return parsed


def _structured_cell_match(
    expected_text: str,
    observed_text: str,
    *,
    numeric_tolerance: float,
) -> bool:
    expected_payload = _maybe_json_cell(expected_text)
    observed_payload = _maybe_json_cell(observed_text)
    if expected_payload is None or observed_payload is None:
        return False
    return _json_like_equal(
        expected_payload,
        observed_payload,
        numeric_tolerance=numeric_tolerance,
    )


def _field_is_tie_equivalent(
    field: str,
    expected_row: dict[str, str],
    observed_row: dict[str, str],
    *,
    numeric_tolerance: float,
) -> bool:
    if field != "most_likely_state":
        return False
    expected_probabilities = _maybe_json_cell(
        expected_row.get("state_probabilities", "")
    )
    observed_probabilities = _maybe_json_cell(
        observed_row.get("state_probabilities", "")
    )
    if not isinstance(expected_probabilities, dict) or not isinstance(
        observed_probabilities, dict
    ):
        return False
    if not _json_like_equal(
        expected_probabilities,
        observed_probabilities,
        numeric_tolerance=numeric_tolerance,
    ):
        return False
    expected_state = expected_row.get(field, "").strip()
    observed_state = observed_row.get(field, "").strip()
    if (
        expected_state not in expected_probabilities
        or observed_state not in expected_probabilities
    ):
        return False
    maximum_probability = max(float(value) for value in expected_probabilities.values())
    return math.isclose(
        float(expected_probabilities[expected_state]),
        maximum_probability,
        rel_tol=numeric_tolerance,
        abs_tol=numeric_tolerance,
    ) and math.isclose(
        float(expected_probabilities[observed_state]),
        maximum_probability,
        rel_tol=numeric_tolerance,
        abs_tol=numeric_tolerance,
    )


def _maybe_json_cell(value: str) -> Any | None:
    if not value or value[0] not in {"{", "["}:
        return None
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return None


def _json_like_equal(
    expected_payload: Any,
    observed_payload: Any,
    *,
    numeric_tolerance: float,
) -> bool:
    if isinstance(expected_payload, dict) and isinstance(observed_payload, dict):
        if set(expected_payload) != set(observed_payload):
            return False
        return all(
            _json_like_equal(
                expected_payload[key],
                observed_payload[key],
                numeric_tolerance=numeric_tolerance,
            )
            for key in expected_payload
        )
    if isinstance(expected_payload, list) and isinstance(observed_payload, list):
        if len(expected_payload) != len(observed_payload):
            return False
        return all(
            _json_like_equal(
                left,
                right,
                numeric_tolerance=numeric_tolerance,
            )
            for left, right in zip(expected_payload, observed_payload, strict=True)
        )
    if isinstance(expected_payload, (int, float)) and isinstance(
        observed_payload, (int, float)
    ):
        return math.isclose(
            float(expected_payload),
            float(observed_payload),
            rel_tol=numeric_tolerance,
            abs_tol=numeric_tolerance,
        )
    return expected_payload == observed_payload


def _load_fasta_records(path: Path) -> dict[str, str]:
    records = {}
    for record in SeqIO.parse(path, "fasta"):
        records[record.id] = str(record.seq)
    return records


@dataclass(slots=True)
class _HtmlContract:
    title: str | None
    headings: list[str]
    local_links: list[str]
    inline_manifest: Any | None


class _ContractHtmlParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self._current_tag: str | None = None
        self._current_attrs: dict[str, str] = {}
        self._text_chunks: list[str] = []
        self.title: str | None = None
        self.headings: list[str] = []
        self.local_links: list[str] = []
        self.inline_manifest: Any | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self._flush_text()
        self._current_tag = tag
        self._current_attrs = {key: value or "" for key, value in attrs}
        for attribute in ("href", "src"):
            value = self._current_attrs.get(attribute)
            if value and _is_local_link(value):
                self.local_links.append(value)

    def handle_endtag(self, tag: str) -> None:
        self._flush_text(force_tag=tag)
        self._current_tag = None
        self._current_attrs = {}

    def handle_data(self, data: str) -> None:
        if self._current_tag is not None:
            self._text_chunks.append(data)

    def _flush_text(self, force_tag: str | None = None) -> None:
        if not self._text_chunks:
            return
        active_tag = force_tag or self._current_tag
        text = "".join(self._text_chunks).strip()
        self._text_chunks.clear()
        if not text or active_tag is None:
            return
        if active_tag == "title":
            self.title = text
        elif active_tag in {"h1", "h2", "h3"}:
            self.headings.append(text)
        elif (
            active_tag == "script"
            and self._current_attrs.get("id") == "bijux-report-manifest"
            and self._current_attrs.get("type") == "application/json"
        ):
            self.inline_manifest = json.loads(text)


def _parse_html_contract(path: Path) -> _HtmlContract:
    parser = _ContractHtmlParser()
    parser.feed(path.read_text(encoding="utf-8"))
    parser.close()
    return _HtmlContract(
        title=parser.title,
        headings=parser.headings,
        local_links=sorted(set(parser.local_links)),
        inline_manifest=parser.inline_manifest,
    )


def _normalize_inline_report_manifest(payload: Any) -> Any:
    if isinstance(payload, dict):
        normalized: dict[str, Any] = {}
        for key, value in payload.items():
            if key == "input_checksums" and isinstance(value, dict):
                normalized[key] = {
                    _normalize_manifest_path(path_text): checksum
                    for path_text, checksum in value.items()
                }
                continue
            if key == "input_paths" and isinstance(value, list):
                normalized[key] = [
                    _normalize_manifest_path(item) if isinstance(item, str) else item
                    for item in value
                ]
                continue
            if key == "path" and isinstance(value, str):
                normalized[key] = _normalize_manifest_path(value)
                continue
            normalized[key] = _normalize_inline_report_manifest(value)
        return normalized
    if isinstance(payload, list):
        return [_normalize_inline_report_manifest(item) for item in payload]
    return payload


def _normalize_manifest_path(path_text: str) -> str:
    normalized = path_text.replace("\\", "/")
    marker = "/tests/fixtures/"
    if marker in normalized:
        suffix = normalized.split(marker, 1)[1]
        return f"tests/fixtures/{suffix}"
    return normalized


def _is_local_link(value: str) -> bool:
    return not (value.startswith(("#", "mailto:")) or "://" in value)
