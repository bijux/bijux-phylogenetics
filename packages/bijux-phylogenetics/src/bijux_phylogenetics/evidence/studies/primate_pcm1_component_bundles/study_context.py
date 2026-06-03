from __future__ import annotations

import csv
import json
from pathlib import Path

from .definitions import REFERENCE_BUNDLE_ID, REFERENCE_SCRIPT_PATH


def study_root(repo_root: Path) -> Path:
    return Path(repo_root) / "evidence-book" / "studies" / "primate-longevity-signal"


def reference_bundle_root(repo_root: Path) -> Path:
    return study_root(repo_root) / REFERENCE_BUNDLE_ID


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def reference_script_locators(line_specs: list[str]) -> list[str]:
    locators: list[str] = []
    for spec in line_specs:
        for part in spec.split(","):
            normalized = part.strip()
            if not normalized:
                continue
            if "-" in normalized:
                start_text, end_text = normalized.split("-", maxsplit=1)
                locators.append(
                    f"{REFERENCE_SCRIPT_PATH}#L{int(start_text)}-L{int(end_text)}"
                )
            else:
                locators.append(f"{REFERENCE_SCRIPT_PATH}#L{int(normalized)}")
    return locators


def load_reference_context(repo_root: Path) -> dict[str, object]:
    bundle_root = reference_bundle_root(repo_root)
    return {
        "r_results": read_json(bundle_root / "results" / "r_reference_results.json"),
        "bijux_results": read_json(
            bundle_root / "results" / "bijux_reference_results.json"
        ),
        "reference_rows": read_csv_rows(
            study_root(repo_root) / "datasets" / "reference_primate.csv"
        ),
        "block_payloads": {
            path.stem: read_json(path)
            for path in sorted(
                (bundle_root / "results" / "block-payloads").glob("*.json")
            )
        },
    }


def missing_counts(rows: list[dict[str, str]]) -> dict[str, int]:
    counts = dict.fromkeys(rows[0].keys(), 0)
    for row in rows:
        for column, value in row.items():
            if value == "":
                counts[column] += 1
    return counts
