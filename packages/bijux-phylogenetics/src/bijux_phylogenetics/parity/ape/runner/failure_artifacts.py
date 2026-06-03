from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil

from ..registry import ApeParityCase


def _copy_if_exists(source: Path, destination: Path) -> None:
    if source.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def _write_rows_table(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def _persist_failure_bundle(
    *,
    failure_root: Path,
    case: ApeParityCase,
    case_file: Path,
    execution_root: Path,
    execution_payload: dict[str, object] | None,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
    reference_error: dict[str, object] | None,
    bijux_error: dict[str, object] | None,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
    bijux_normalized_text: str | None,
    mismatch_reason: str,
) -> Path:
    artifact_root = failure_root / case.case_id
    artifact_root.mkdir(parents=True, exist_ok=True)
    _copy_if_exists(case_file, artifact_root / "case.json")
    _copy_if_exists(
        execution_root.parent / "bijux-reference-input.nwk",
        artifact_root / "bijux-reference-input.nwk",
    )
    _copy_if_exists(
        execution_root / "reference-execution.json",
        artifact_root / "reference-execution.json",
    )
    if execution_payload is not None:
        outputs = execution_payload.get("outputs")
        if isinstance(outputs, dict):
            for path_string in outputs.values():
                if isinstance(path_string, str):
                    source = Path(path_string)
                    _copy_if_exists(source, artifact_root / f"reference-{source.name}")
    if execution_payload is not None:
        _write_json(
            artifact_root / "reference-execution.observed.json", execution_payload
        )
    if reference_summary is not None:
        _write_json(
            artifact_root / "reference-summary.observed.json", reference_summary
        )
    if reference_rows is not None:
        _write_json(artifact_root / "reference-rows.observed.json", reference_rows)
        _write_rows_table(artifact_root / "reference-rows.observed.tsv", reference_rows)
    if bijux_summary is not None:
        _write_json(artifact_root / "bijux-summary.json", bijux_summary)
    if reference_error is not None:
        _write_json(artifact_root / "reference-error.observed.json", reference_error)
    if bijux_error is not None:
        _write_json(artifact_root / "bijux-error.json", bijux_error)
    if bijux_rows is not None:
        _write_json(artifact_root / "bijux-rows.json", bijux_rows)
        _write_rows_table(artifact_root / "bijux-rows.tsv", bijux_rows)
    if bijux_normalized_text is not None:
        (artifact_root / "bijux-normalized.txt").write_text(
            f"{bijux_normalized_text}\n",
            encoding="utf-8",
        )
    _write_json(
        artifact_root / "comparison.json",
        {
            "case_id": case.case_id,
            "function_name": case.function_name,
            "mismatch_reason": mismatch_reason,
        },
    )
    return artifact_root
