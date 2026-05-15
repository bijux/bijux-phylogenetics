from __future__ import annotations

from pathlib import Path

from tests.support.fake_external_engines import write_executable


def fake_phytools_rscript(
    path: Path,
    *,
    phytools_available: bool = True,
    summary_overrides: dict[str, dict[str, object]] | None = None,
) -> Path:
    summary_payload = repr(summary_overrides or {})
    return write_executable(
        path,
        """#!/usr/bin/env python3
import csv
import json
import sys
from pathlib import Path


SUMMARIES = {
    "phylosig-lambda-non-ultrametric-strong-signal-twenty-four-taxa": {
        "taxon_count": 24,
        "trait_name": "signal_strong",
        "lambda_value": 1.0,
        "log_likelihood": -27.527955447046294,
        "null_log_likelihood": -38.15175747253939,
        "p_value": 1.0,
    },
    "phylosig-k-strong-signal-twenty-four-taxa": {
        "taxon_count": 24,
        "trait_name": "signal_strong",
        "k": 0.881409701146833,
        "p_value": 0.593,
        "permutation_count": 1000,
        "simulated_k_minimum": 0.15,
        "simulated_k_mean": 0.64,
        "simulated_k_maximum": 1.37,
    },
    "phylosig-lambda-primate-longevity": {
        "taxon_count": 90,
        "trait_name": "longevity",
        "lambda_value": 0.99957,
        "log_likelihood": -170.97952,
        "null_log_likelihood": -193.33681,
        "p_value": 0.0,
    },
}

SUMMARIES.update(__SUMMARY_OVERRIDES__)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_summary_table(path: Path, summary: dict[str, object]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\\t")
        writer.writeheader()
        for key, value in summary.items():
            writer.writerow({"metric": key, "value": value})


case_payload = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
output_root = Path(sys.argv[3])
output_root.mkdir(parents=True, exist_ok=True)
execution_path = output_root / "reference-execution.json"
summary_path = output_root / "reference-summary.json"
summary_table_path = output_root / "reference-summary.tsv"

if not __PHYTOOLS_AVAILABLE__:
    write_json(
        execution_path,
        {
            "status": "unavailable",
            "mismatch_reason": "phytools_package_unavailable",
            "r_version": "4.6.0",
            "phytools_version": None,
        },
    )
    raise SystemExit(0)

case_id = case_payload["case_id"]
if case_id not in SUMMARIES:
    write_json(
        execution_path,
        {
            "status": "failed",
            "mismatch_reason": "unsupported_operation",
            "error_type": "ValueError",
            "message": "unsupported phytools parity case: " + case_id,
            "r_version": "4.6.0",
            "phytools_version": "2.5.2",
        },
    )
    raise SystemExit(0)

summary = SUMMARIES[case_id]
write_json(
    execution_path,
    {
        "status": "ok",
        "r_version": "4.6.0",
        "phytools_version": "2.5.2",
    },
)
write_json(summary_path, summary)
write_summary_table(summary_table_path, summary)
"""
        .replace("__SUMMARY_OVERRIDES__", summary_payload)
        .replace("__PHYTOOLS_AVAILABLE__", repr(phytools_available)),
    )
