"""Render governed local artifact surfaces for each evidence bundle."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

REQUIRED_BUNDLE_LOCAL_ARTIFACTS = (
    "reference.R",
    "analysis.py",
    "checks.json",
    "report.md",
    "provenance.json",
)
ARTIFACT_JSON_FILENAMES = {"checks.json", "provenance.json"}
PRIMARY_OUTPUT_EXCLUDED_FILENAMES = {
    "README.md",
    "manifest.json",
    "claims.json",
    "claim_verdicts.json",
    "reviewer-summary.json",
    "reviewer-summary.md",
    "inputs.manifest.json",
    *REQUIRED_BUNDLE_LOCAL_ARTIFACTS,
}


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object at {path}")
    return payload


def iter_bundle_roots(repo_root: Path) -> list[Path]:
    studies_root = repo_root / "evidence-book" / "studies"
    return sorted(path for path in studies_root.glob("*/evidence-*") if path.is_dir())


def _relative(path: Path, repo_root: Path) -> str:
    return path.relative_to(repo_root).as_posix()


def _primary_output_paths(bundle_root: Path, repo_root: Path) -> list[str]:
    output_paths: list[str] = []
    for child in sorted(bundle_root.iterdir()):
        if not child.is_file() or child.name in PRIMARY_OUTPUT_EXCLUDED_FILENAMES:
            continue
        output_paths.append(_relative(child, repo_root))
    return output_paths


def _study_manifest(study_root: Path) -> dict[str, Any]:
    return _load_json(study_root / "study.json")


def _bundle_manifest(bundle_root: Path) -> dict[str, Any]:
    return _load_json(bundle_root / "manifest.json")


def _study_support_scripts(study_root: Path, repo_root: Path) -> dict[str, Any]:
    reference_root = study_root / "reference"
    reference_r_scripts = [
        _relative(path, repo_root)
        for path in sorted(reference_root.glob("*.R"))
        if path.is_file()
    ]
    reference_python_scripts = [
        _relative(path, repo_root)
        for path in sorted(reference_root.glob("*.py"))
        if path.is_file()
    ]
    build_script_path = study_root / "build_evidence.py"
    return {
        "build_script_path": None
        if not build_script_path.is_file()
        else _relative(build_script_path, repo_root),
        "reference_r_scripts": reference_r_scripts,
        "reference_python_scripts": reference_python_scripts,
    }


def _source_basis_locators(bundle_manifest: dict[str, Any]) -> list[str]:
    locators: list[str] = []
    for entry in bundle_manifest.get("source_basis", []):
        if isinstance(entry, dict):
            locator = entry.get("locator")
            if isinstance(locator, str) and locator:
                locators.append(locator)
    return locators


def build_bundle_artifact_contract(repo_root: Path, bundle_root: Path) -> dict[str, Any]:
    repo_root = repo_root.resolve()
    bundle_root = bundle_root.resolve()
    study_root = bundle_root.parent
    study_manifest = _study_manifest(study_root)
    bundle_manifest = _bundle_manifest(bundle_root)
    support_scripts = _study_support_scripts(study_root, repo_root)
    provenance_path = study_root / "provenance"
    provenance_locators = [
        _relative(path, repo_root)
        for path in sorted(provenance_path.glob("*.json"))
        if path.is_file()
    ]
    return {
        "study_id": str(bundle_manifest["study_id"]),
        "evidence_id": str(bundle_manifest["evidence_id"]),
        "study_title": str(study_manifest["study_title"]),
        "evidence_title": str(bundle_manifest["evidence_title"]),
        "summary": str(bundle_manifest["summary"]),
        "comparison_mode": str(bundle_manifest["comparison_mode"]),
        "claim_ids": list(bundle_manifest.get("claim_ids", [])),
        "claim_tags": list(bundle_manifest.get("claim_tags", [])),
        "verdict": dict(bundle_manifest.get("verdict", {})),
        "limitations": list(bundle_manifest.get("limitations", [])),
        "source_intake_policy": str(study_manifest.get("source_intake_policy", "")),
        "study_provenance_locator": str(
            study_manifest.get("provenance_descriptor_locator", "")
        ),
        "study_dataset_registry_locator": str(
            study_manifest.get("dataset_registry_locator", "")
        ),
        "source_basis_locators": _source_basis_locators(bundle_manifest),
        "primary_output_paths": _primary_output_paths(bundle_root, repo_root),
        "bundle_relative_path": _relative(bundle_root, repo_root),
        "build_script_path": support_scripts["build_script_path"],
        "reference_r_scripts": support_scripts["reference_r_scripts"],
        "reference_python_scripts": support_scripts["reference_python_scripts"],
        "bundle_provenance_paths": provenance_locators,
    }


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
        '  library(jsonlite)',
        "})",
        "",
        'script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])',
        "bundle_root <- dirname(normalizePath(script_path, mustWork = TRUE))",
        'args <- commandArgs(trailingOnly = TRUE)',
        'out_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = FALSE) else tempdir()',
        'dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)',
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
            'cat(toJSON(payload, auto_unbox = TRUE, pretty = TRUE))',
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
        'REPO_ROOT = Path(__file__).resolve().parents[4]',
        'BUNDLE_ROOT = Path(__file__).resolve().parent',
        f'STUDY_ID = "{contract["study_id"]}"',
        f'EVIDENCE_ID = "{contract["evidence_id"]}"',
        f'COMPARISON_MODE = "{contract["comparison_mode"]}"',
        f"PRIMARY_OUTPUTS = {output_paths!r}",
        f"BUILD_SCRIPT = {build_script!r}",
        "",
        "def main() -> None:",
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
        "    print(json.dumps(payload, indent=2, sort_keys=True))",
        "",
        'if __name__ == "__main__":',
        "    main()",
        "",
    ]
    return "\n".join(lines)


def render_checks_json(contract: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "study_id": contract["study_id"],
        "evidence_id": contract["evidence_id"],
        "comparison_mode": contract["comparison_mode"],
        "expected_verdict_status": contract["verdict"].get("status"),
        "claim_ids": contract["claim_ids"],
        "required_local_artifacts": list(REQUIRED_BUNDLE_LOCAL_ARTIFACTS),
        "primary_output_paths": contract["primary_output_paths"],
        "source_basis_locators": contract["source_basis_locators"],
        "validation_rules": [
            {
                "rule_id": "bundle-local-artifacts-present",
                "pass_when": "reference, analysis, checks, report, and provenance files are all present in the evidence directory",
            },
            {
                "rule_id": "primary-outputs-present",
                "pass_when": "the bundle retains its governed machine outputs alongside the authored local artifact surfaces",
            },
            {
                "rule_id": "verdict-still-explicit",
                "pass_when": "the bundle manifest continues to expose the verdict status and summary without hiding mismatches or boundaries",
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
        "## Claims",
        "",
    ]
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
            "## Primary Outputs",
            "",
        ]
    )
    for path in contract["primary_output_paths"]:
        lines.append(f"- `{path}`")
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
        "schema_version": 1,
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
                "artifact_kind": filename,
                "generation_mode": "governed-bundle-template",
            }
            for filename in REQUIRED_BUNDLE_LOCAL_ARTIFACTS
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
