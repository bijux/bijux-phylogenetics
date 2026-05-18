from __future__ import annotations

from pathlib import Path

from tests.support.fake_external_engines import write_executable


def fake_geiger_rscript(
    path: Path,
    *,
    geiger_available: bool = True,
    summary_overrides: dict[str, dict[str, object]] | None = None,
    reference_payloads: dict[str, dict[str, object]] | None = None,
) -> Path:
    summary_payload = repr(summary_overrides or {})
    reference_payload_repr = repr(reference_payloads or {})
    script = (
        """#!/usr/bin/env python3
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

SUMMARY_OVERRIDES = __SUMMARY_OVERRIDES__
REFERENCE_PAYLOADS = __REFERENCE_PAYLOADS__


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_rows_table(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["parameter", "value"], delimiter="\\t")
        writer.writeheader()
        writer.writerows(rows)


def find_package_root() -> Path:
    for candidate in (Path.cwd(), *Path.cwd().parents):
        package_root = candidate / "packages" / "bijux-phylogenetics"
        if package_root.is_dir():
            return package_root
    raise ValueError("could not locate bijux-phylogenetics package root")


def find_repo_python(package_root: Path) -> Path:
    candidates = (
        package_root / ".venv" / "bin" / "python",
        package_root.parent / ".venv" / "bin" / "python",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise ValueError("could not locate repo python for fake geiger parity payload")


def parameter_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    rows = []
    for parameter in ("root_state", "rate", "log_likelihood", "aic", "aicc"):
        value = summary.get(parameter)
        if value is None:
            continue
        rows.append({"parameter": parameter, "value": value})
    parameter_name = summary.get("parameter_name")
    parameter_value = summary.get("parameter_value")
    if parameter_name and parameter_value is not None:
        rows.append({"parameter": parameter_name, "value": parameter_value})
    return rows


def standard_error_policy() -> str:
    return "tip-standard-errors-not-supported"


def missing_value_policy() -> str:
    return "prune-tree-tip-overlap-with-missing-or-nonnumeric-trait-values"


def build_reference_payload(case_payload: dict[str, object]) -> tuple[dict[str, object], list[dict[str, object]]]:
    if case_payload["case_id"] in REFERENCE_PAYLOADS:
        summary = dict(REFERENCE_PAYLOADS[case_payload["case_id"]]["summary"])
        rows = [dict(row) for row in REFERENCE_PAYLOADS[case_payload["case_id"]]["rows"]]
        return summary, rows
    package_root = find_package_root()
    repo_python = find_repo_python(package_root)
    payload_path = package_root / "artifacts" / "fake-geiger-case.json"
    payload_path.parent.mkdir(parents=True, exist_ok=True)
    payload_path.write_text(json.dumps(case_payload), encoding="utf-8")
    inline_script = '''
import json
from pathlib import Path
from bijux_phylogenetics.comparative.common import summarize_numeric_trait_readiness
from bijux_phylogenetics.comparative.evolutionary_modes import fit_continuous_evolutionary_mode

case_payload = json.loads(Path(__PAYLOAD_PATH__).read_text(encoding="utf-8"))
mode_lookup = {
    "BM": "brownian",
    "OU": "ornstein-uhlenbeck",
    "EB": "early-burst",
}
tree_path, traits_path = case_payload["input_fixtures"]
readiness = summarize_numeric_trait_readiness(
    Path(tree_path),
    Path(traits_path),
    trait=case_payload["trait_name"],
    taxon_column=case_payload["taxon_column"],
)
report = fit_continuous_evolutionary_mode(
    Path(tree_path),
    Path(traits_path),
    trait=case_payload["trait_name"],
    mode=mode_lookup[case_payload["model_name"]],
    taxon_column=case_payload["taxon_column"],
    ou_bounds=(0.0, 10.0),
    early_burst_bounds=(0.0, 10.0),
)
excluded_taxa = sorted(
    {
        *readiness.missing_from_traits,
        *readiness.pruned_missing_value_taxa,
        *readiness.pruned_non_numeric_taxa,
    }
)
summary = {
    "taxon_count": report.taxon_count,
    "trait_name": report.trait,
    "model_name": case_payload["model_name"],
    "excluded_taxon_count": len(excluded_taxa),
    "excluded_taxa": excluded_taxa,
    "missing_value_taxa": list(readiness.pruned_missing_value_taxa),
    "non_numeric_taxa": list(readiness.pruned_non_numeric_taxa),
    "missing_from_traits": list(readiness.missing_from_traits),
    "extra_trait_taxa": list(readiness.extra_trait_taxa),
    "missing_value_policy": "__MISSING_VALUE_POLICY__",
    "standard_error_policy": "__STANDARD_ERROR_POLICY__",
    "root_state": report.root_state,
    "rate": report.rate,
    "log_likelihood": report.log_likelihood,
    "aic": report.aic,
    "aicc": report.aicc,
    "parameter_name": report.parameter_name,
    "parameter_value": report.parameter_value,
    "optimizer_settings": case_payload["optimizer_settings"],
    "optimizer_result": {
        "convergence_code": 0,
        "message": "fake geiger parity runner",
    },
}
print(json.dumps(summary))
'''
    command = [
        str(repo_python),
        "-c",
        inline_script
        .replace("__PAYLOAD_PATH__", repr(str(payload_path)))
        .replace("__MISSING_VALUE_POLICY__", missing_value_policy())
        .replace("__STANDARD_ERROR_POLICY__", standard_error_policy()),
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_root / "src")
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(package_root),
        env=env,
    )
    summary = json.loads(result.stdout)
    if case_payload["case_id"] in SUMMARY_OVERRIDES:
        summary.update(SUMMARY_OVERRIDES[case_payload["case_id"]])
    return summary, parameter_rows(summary)


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit(2)
    case_path = Path(sys.argv[2])
    output_root = Path(sys.argv[3])
    output_root.mkdir(parents=True, exist_ok=True)
    execution_path = output_root / "reference-execution.json"
    summary_path = output_root / "reference-summary.json"
    rows_path = output_root / "reference-parameters.tsv"
    case_payload = json.loads(case_path.read_text(encoding="utf-8"))
    if not __GEIGER_AVAILABLE__:
        write_json(
            execution_path,
            {
                "status": "unavailable",
                "mismatch_reason": "geiger_package_unavailable",
                "r_version": "4.4.0",
                "geiger_version": None,
            },
        )
        return 0
    summary, rows = build_reference_payload(case_payload)
    write_json(summary_path, summary)
    write_rows_table(rows_path, rows)
    write_json(
        execution_path,
        {
            "status": "ok",
            "r_version": "4.4.0",
            "geiger_version": "2.0.11",
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
"""
        .replace("__GEIGER_AVAILABLE__", "True" if geiger_available else "False")
        .replace("__SUMMARY_OVERRIDES__", summary_payload)
        .replace("__REFERENCE_PAYLOADS__", reference_payload_repr)
    )
    return write_executable(
        path,
        script,
    )
