from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.simulation import validate_geiger_sim_char_reference_examples


def load_recovery_summary(path: Path) -> dict[str, int | str]:
    rows = _read_tsv_rows(path)
    if len(rows) != 1:
        raise ValueError(f"expected one workflow summary row in {path}")
    row = rows[0]
    payload: dict[str, int | str] = {}
    for key, value in row.items():
        if key == "dataset_id":
            payload[key] = value
            continue
        payload[key] = int(value)
    return payload


def load_sim_char_summary() -> dict[str, object]:
    report = validate_geiger_sim_char_reference_examples()
    return {
        "case_count": report.case_count,
        "all_passed": report.all_passed,
    }


def load_large_tree_benchmark_summary() -> dict[str, int]:
    rows = _read_tsv_rows(
        repository_root()
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "resources"
        / "benchmarks"
        / "large_tree_model_fitting"
        / "expected"
        / "summary.tsv"
    )
    if len(rows) != 1:
        raise ValueError("expected one large-tree benchmark summary row")
    row = rows[0]
    return {
        "case_count": int(row["case_count"]),
        "geiger_match_case_count": int(row["geiger_match_case_count"]),
        "threshold_pass_case_count": int(row["threshold_pass_case_count"]),
        "too_slow_case_count": int(row["too_slow_case_count"]),
        "unstable_case_count": int(row["unstable_case_count"]),
    }


def load_real_dataset_benchmark_summary() -> dict[str, int]:
    expected_root = (
        repository_root()
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "resources"
        / "benchmarks"
        / "real_dataset_macroevolution"
        / "expected"
    )
    summary_rows = _read_tsv_rows(expected_root / "benchmark-summary.tsv")
    model_rows = _read_tsv_rows(expected_root / "model-table.tsv")
    alignment_review_rows = _read_tsv_rows(expected_root / "alignment-review.tsv")
    parity_rows = _read_tsv_rows(expected_root / "geiger-parity.tsv")
    native_summary_rows = [
        row for row in summary_rows if row["review_scope"] == "native-model-table"
    ]
    return {
        "summary_row_count": len(summary_rows),
        "model_row_count": len(model_rows),
        "alignment_review_row_count": len(alignment_review_rows),
        "parity_row_count": len(parity_rows),
        "selection_match_count": sum(
            1
            for row in native_summary_rows
            if _tsv_bool(row["selection_matches_geiger"])
        ),
        "unstable_review_count": sum(
            1
            for row in native_summary_rows
            if not _tsv_bool(row["stable_conclusion_supported"])
        ),
    }


def repository_root() -> Path:
    return Path(__file__).resolve().parents[7]


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _tsv_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise ValueError(f"expected boolean TSV cell, found {value!r}")
