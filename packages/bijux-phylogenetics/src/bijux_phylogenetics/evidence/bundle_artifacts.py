"""Render governed local artifact surfaces for each evidence bundle."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .bundle_contracts import (
    LOCAL_ARTIFACT_PURPOSES,
    REQUIRED_BUNDLE_LOCAL_ARTIFACTS,
    build_bundle_artifact_contract,
)
from .bundle_results import build_bundle_result_artifacts


def render_reference_r(contract: dict[str, Any]) -> str:
    reference_scripts = contract["reference_r_scripts"]
    source_locators = contract["source_basis_locators"]
    execution_mode = (
        "study_reference_wrapper" if reference_scripts else "not_applicable"
    )
    lines = [
        "#!/usr/bin/env Rscript",
        "",
        "suppressPackageStartupMessages({",
        "  library(jsonlite)",
        "})",
        "",
        'script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])',
        "bundle_root <- dirname(normalizePath(script_path, mustWork = TRUE))",
        'results_root <- file.path(bundle_root, "results")',
        "dir.create(results_root, recursive = TRUE, showWarnings = FALSE)",
        "args <- commandArgs(trailingOnly = TRUE)",
        "out_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = FALSE) else results_root",
        "dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)",
        "",
        "payload <- list(",
        f'  study_id = "{contract["study_id"]}",',
        f'  evidence_id = "{contract["evidence_id"]}",',
        f'  evidence_title = "{contract["evidence_title"]}",',
        f'  comparison_mode = "{contract["comparison_mode"]}",',
        f'  execution_mode = "{execution_mode}",',
        f'  source_intake_policy = "{contract["source_intake_policy"]}",',
        "  source_basis_locators = c(",
    ]
    for locator in source_locators:
        lines.append(f'    "{locator}",')
    lines.extend(
        [
            "  ),",
            "  reference_scripts = c(",
        ]
    )
    for locator in reference_scripts:
        lines.append(f'    "{locator}",')
    lines.extend(
        [
            "  )",
            ")",
            "",
            'write_json(payload, file.path(out_dir, "reference-contract.json"), auto_unbox = TRUE, pretty = TRUE)',
            "cat(toJSON(payload, auto_unbox = TRUE, pretty = TRUE))",
            'cat("\\n")',
            "",
        ]
    )
    if reference_scripts:
        lines.extend(
            [
                "# To run the governed study reference, invoke the study-owned R script",
                "# listed in `reference_scripts` with the source context required by that",
                "# study. This bundle-local wrapper stays focused on one evidence unit and",
                "# records the exact comparison contract without editing the untouched",
                "# Lund materials in place.",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "# This evidence unit is a Bijux-native or fixture-backed surface with no",
                "# canonical external R comparison program. The wrapper therefore records",
                "# the bundle-local contract instead of pretending an R parity run exists.",
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def render_analysis_py(contract: dict[str, Any]) -> str:
    build_script = contract["build_script_path"]
    output_paths = contract["primary_output_paths"]
    lines = [
        "from __future__ import annotations",
        "",
        "import json",
        "from pathlib import Path",
        "import subprocess",
        "import sys",
        "",
        "REPO_ROOT = Path(__file__).resolve().parents[4]",
        "BUNDLE_ROOT = Path(__file__).resolve().parent",
        "RESULTS_ROOT = BUNDLE_ROOT / 'results'",
        f'STUDY_ID = "{contract["study_id"]}"',
        f'EVIDENCE_ID = "{contract["evidence_id"]}"',
        f'COMPARISON_MODE = "{contract["comparison_mode"]}"',
        f"PRIMARY_OUTPUTS = {output_paths!r}",
        f"BUILD_SCRIPT = {build_script!r}",
        "",
        "def main() -> None:",
        "    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)",
        "    payload = {",
        '        "study_id": STUDY_ID,',
        '        "evidence_id": EVIDENCE_ID,',
        '        "comparison_mode": COMPARISON_MODE,',
        '        "build_script": BUILD_SCRIPT,',
        "    }",
        "    if BUILD_SCRIPT is not None:",
        "        subprocess.run(",
        "            [sys.executable, str(REPO_ROOT / BUILD_SCRIPT)],",
        "            cwd=str(REPO_ROOT),",
        "            check=True,",
        "        )",
        '        payload["execution_mode"] = "study_build_wrapper"',
        "    else:",
        '        payload["execution_mode"] = "bundle_contract_only"',
        '    payload["primary_outputs"] = [',
        "        path for path in PRIMARY_OUTPUTS if (REPO_ROOT / path).is_file()",
        "    ]",
        "    output_path = RESULTS_ROOT / 'analysis-run.json'",
        "    output_path.write_text(",
        "        json.dumps(payload, indent=2, sort_keys=True) + '\\n',",
        "        encoding='utf-8',",
        "    )",
        "    print(json.dumps(payload, indent=2, sort_keys=True))",
        "",
        "if __name__ == '__main__':",
        "    main()",
        "",
    ]
    return "\n".join(lines)


def render_checks_json(contract: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": 2,
        "study_id": contract["study_id"],
        "evidence_id": contract["evidence_id"],
        "comparison_mode": contract["comparison_mode"],
        "expected_verdict_status": contract["verdict"].get("status"),
        "claim_ids": contract["claim_ids"],
        "required_local_artifacts": list(REQUIRED_BUNDLE_LOCAL_ARTIFACTS),
        "required_result_artifacts": list(contract["required_result_artifacts"]),
        "local_artifact_purposes": dict(LOCAL_ARTIFACT_PURPOSES),
        "primary_output_paths": contract["primary_output_paths"],
        "source_basis_locators": contract["source_basis_locators"],
        "validation_rules": [
            {
                "rule_id": "bundle-local-artifacts-present",
                "success_condition": "reference, analysis, checks, report, and provenance files are all present in the evidence directory",
            },
            {
                "rule_id": "results-directory-governed",
                "success_condition": "the evidence-owned results directory contains a manifest and README that inventory local execution products and governed outputs",
            },
            {
                "rule_id": "primary-outputs-present",
                "success_condition": "the bundle retains its governed machine outputs alongside the authored local artifact surfaces",
            },
            {
                "rule_id": "verdict-still-explicit",
                "success_condition": "the bundle manifest continues to expose the verdict status and summary without hiding mismatches or boundaries",
            },
        ],
    }


def render_report_md(contract: dict[str, Any]) -> str:
    lines = [
        f"# {contract['evidence_title']}",
        "",
        contract["summary"],
        "",
        f"- study: `{contract['study_id']}`",
        f"- evidence: `{contract['evidence_id']}`",
        f"- comparison mode: `{contract['comparison_mode']}`",
        f"- expected verdict: `{contract['verdict'].get('status', 'unknown')}`",
        "",
        "## Local Artifacts",
        "",
    ]
    for filename, purpose in LOCAL_ARTIFACT_PURPOSES.items():
        lines.append(f"- `{filename}`: {purpose}")
    lines.extend(
        [
            "",
            "## Claims",
            "",
        ]
    )
    for claim_id in contract["claim_ids"]:
        lines.append(f"- `{claim_id}`")
    lines.extend(
        [
            "",
            "## Source Basis",
            "",
        ]
    )
    for locator in contract["source_basis_locators"]:
        lines.append(f"- `{locator}`")
    lines.extend(
        [
            "",
            "## Governed Primary Outputs",
            "",
        ]
    )
    for path in contract["primary_output_paths"]:
        lines.append(f"- `{path}`")
    lines.extend(
        [
            "",
            "## Results Directory",
            "",
            f"- `{contract['results_relative_path']}/README.md`",
            f"- `{contract['results_relative_path']}/manifest.json`",
        ]
    )
    if contract["limitations"]:
        lines.extend(
            [
                "",
                "## Limits",
                "",
            ]
        )
        for item in contract["limitations"]:
            lines.append(f"- {item}")
    lines.append("")
    return "\n".join(lines)


def render_provenance_json(contract: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": 2,
        "study_id": contract["study_id"],
        "evidence_id": contract["evidence_id"],
        "source_intake_policy": contract["source_intake_policy"],
        "study_provenance_locator": contract["study_provenance_locator"],
        "study_dataset_registry_locator": contract["study_dataset_registry_locator"],
        "bundle_provenance_paths": contract["bundle_provenance_paths"],
        "source_basis_locators": contract["source_basis_locators"],
        "authored_local_artifacts": [
            {
                "artifact_path": f"{contract['bundle_relative_path']}/{filename}",
                "artifact_kind": LOCAL_ARTIFACT_PURPOSES[filename],
                "generation_mode": "governed-bundle-template",
            }
            for filename in REQUIRED_BUNDLE_LOCAL_ARTIFACTS
        ],
        "bundle_result_surfaces": [
            {
                "artifact_path": path,
                "artifact_kind": "governed-result-surface",
            }
            for path in contract["required_result_artifacts"]
        ],
        "study_support_scripts": {
            "build_script_path": contract["build_script_path"],
            "reference_r_scripts": contract["reference_r_scripts"],
            "reference_python_scripts": contract["reference_python_scripts"],
        },
    }


def build_bundle_local_artifacts(
    repo_root: Path,
    bundle_root: Path,
) -> dict[str, str | dict[str, Any]]:
    contract = build_bundle_artifact_contract(repo_root, bundle_root)
    return {
        "reference.R": render_reference_r(contract),
        "analysis.py": render_analysis_py(contract),
        "checks.json": render_checks_json(contract),
        "report.md": render_report_md(contract),
        "provenance.json": render_provenance_json(contract),
    }


def build_bundle_governed_artifacts(
    repo_root: Path,
    bundle_root: Path,
) -> dict[str, str | dict[str, Any]]:
    contract = build_bundle_artifact_contract(repo_root, bundle_root)
    artifacts = build_bundle_local_artifacts(repo_root, bundle_root)
    artifacts.update(build_bundle_result_artifacts(repo_root, bundle_root, contract))
    return artifacts
