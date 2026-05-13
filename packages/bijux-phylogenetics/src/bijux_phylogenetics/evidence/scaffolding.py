from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from .bundle_contracts import LOCAL_ARTIFACT_PURPOSES

RESULTS_DIRNAME = "results"
COMPARISON_OBSERVATION_FILENAME = "comparison-observation.json"


@dataclass(frozen=True, slots=True)
class EvidenceBundleTemplateSpec:
    study_id: str
    evidence_id: str
    evidence_title: str
    summary: str
    owner_package: str
    claim_id: str
    claim_title: str
    claim_summary: str
    analytical_surfaces: tuple[str, ...]
    claim_tags: tuple[str, ...]
    source_basis_locators: tuple[str, ...]
    governed_code_paths: tuple[str, ...]
    limitation: str
    freshness_date: str


def build_evidence_bundle_template(
    spec: EvidenceBundleTemplateSpec,
) -> dict[str, object]:
    manifest = {
        "schema_version": 1,
        "study_id": spec.study_id,
        "evidence_id": spec.evidence_id,
        "evidence_title": spec.evidence_title,
        "summary": spec.summary,
        "owner_package": spec.owner_package,
        "claim_ids": [spec.claim_id],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": f"source-{index + 1}",
                "locator": locator,
            }
            for index, locator in enumerate(spec.source_basis_locators)
        ],
        "freshness": {
            "last_generated_on": spec.freshness_date,
            "governed_code_paths": list(spec.governed_code_paths),
            "source_basis_locators": list(spec.source_basis_locators),
        },
        "ownership": {
            "owner_package": spec.owner_package,
            "analytical_surfaces": list(spec.analytical_surfaces),
        },
        "claim_tags": list(spec.claim_tags),
        "verdict": {
            "status": "not_comparable",
            "summary": "Template placeholder until governed evidence is generated.",
        },
        "limitations": [spec.limitation],
    }
    claims = {
        "schema_version": 1,
        "study_id": spec.study_id,
        "evidence_id": spec.evidence_id,
        "claim_count": 1,
        "claims": [
            {
                "claim_id": spec.claim_id,
                "claim_title": spec.claim_title,
                "summary": spec.claim_summary,
                "verdict": "not_comparable",
                "evidence_ids": [spec.evidence_id],
                "source_fragments": [],
            }
        ],
    }
    readme = "\n".join(
        [
            f"# {spec.evidence_title}",
            "",
            spec.summary,
            "",
            "Template checklist:",
            "",
            "- replace the placeholder report payload with governed evidence",
            "- replace the placeholder verdict once parity or boundary evidence exists",
            "- keep source locators repository-relative or explicitly externalized",
            "",
        ]
    )
    report_payload = {
        "schema_version": 1,
        "study_id": spec.study_id,
        "evidence_id": spec.evidence_id,
        "status": "template",
        "next_action": "replace this placeholder payload with governed evidence output",
    }
    reference_r = "\n".join(
        [
            "#!/usr/bin/env Rscript",
            "",
            "payload <- list(",
            f'  study_id = "{spec.study_id}",',
            f'  evidence_id = "{spec.evidence_id}",',
            '  execution_mode = "template-contract"',
            ")",
            "cat(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE))",
            'cat("\\n")',
            "",
        ]
    )
    analysis_py = "\n".join(
        [
            "from __future__ import annotations",
            "",
            "import json",
            "",
            "",
            "def main() -> None:",
            "    payload = {",
            f'        "study_id": "{spec.study_id}",',
            f'        "evidence_id": "{spec.evidence_id}",',
            '        "execution_mode": "template-contract",',
            "    }",
            "    print(json.dumps(payload, indent=2, sort_keys=True))",
            "",
            "",
            "if __name__ == '__main__':",
            "    main()",
            "",
        ]
    )
    checks = {
        "schema_version": 2,
        "study_id": spec.study_id,
        "evidence_id": spec.evidence_id,
        "expected_verdict_status": "not_comparable",
        "required_local_artifacts": list(LOCAL_ARTIFACT_PURPOSES),
        "required_result_artifacts": [
            f"{RESULTS_DIRNAME}/README.md",
            f"{RESULTS_DIRNAME}/manifest.json",
        ],
    }
    human_report = "\n".join(
        [
            f"# {spec.evidence_title}",
            "",
            spec.summary,
            "",
            "## Template Status",
            "",
            "- verdict: `not_comparable`",
            "- next action: replace placeholder local scripts and result observation with governed evidence",
            "",
        ]
    )
    provenance = {
        "schema_version": 2,
        "study_id": spec.study_id,
        "evidence_id": spec.evidence_id,
        "source_basis_locators": list(spec.source_basis_locators),
        "authored_local_artifacts": [
            {
                "artifact_path": filename,
                "artifact_kind": purpose,
            }
            for filename, purpose in LOCAL_ARTIFACT_PURPOSES.items()
        ],
        "bundle_result_surfaces": [
            {
                "artifact_path": f"{RESULTS_DIRNAME}/README.md",
                "artifact_kind": "governed-result-surface",
            },
            {
                "artifact_path": f"{RESULTS_DIRNAME}/manifest.json",
                "artifact_kind": "governed-result-surface",
            },
        ],
    }
    results_readme = "\n".join(
        [
            "# Results",
            "",
            "This directory stores evidence-local execution products and governed",
            "output inventories for this evidence bundle.",
            "",
        ]
    )
    results_manifest = {
        "schema_version": 1,
        "study_id": spec.study_id,
        "evidence_id": spec.evidence_id,
        "results_directory": RESULTS_DIRNAME,
        "governed_primary_output_count": 1,
        "governed_primary_outputs": [
            f"{spec.evidence_id}/{RESULTS_DIRNAME}/{COMPARISON_OBSERVATION_FILENAME}"
        ],
    }
    return {
        "manifest": manifest,
        "claims": claims,
        "readme": readme,
        "reference_r": reference_r,
        "analysis_py": analysis_py,
        "checks": checks,
        "human_report": human_report,
        "provenance": provenance,
        "results_readme": results_readme,
        "results_manifest": results_manifest,
        "report_filename": f"{RESULTS_DIRNAME}/{COMPARISON_OBSERVATION_FILENAME}",
        "report_payload": report_payload,
    }


def write_evidence_bundle_template(
    output_root: Path, spec: EvidenceBundleTemplateSpec
) -> dict[str, Path]:
    output_root.mkdir(parents=True, exist_ok=True)
    template = build_evidence_bundle_template(spec)
    manifest_path = output_root / "manifest.json"
    claims_path = output_root / "claims.json"
    readme_path = output_root / "README.md"
    reference_r_path = output_root / "reference.R"
    analysis_py_path = output_root / "analysis.py"
    checks_path = output_root / "checks.json"
    human_report_path = output_root / "report.md"
    provenance_path = output_root / "provenance.json"
    results_readme_path = output_root / RESULTS_DIRNAME / "README.md"
    results_manifest_path = output_root / RESULTS_DIRNAME / "manifest.json"
    report_path = output_root / str(template["report_filename"])
    report_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(template["manifest"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    claims_path.write_text(
        json.dumps(template["claims"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    readme_path.write_text(str(template["readme"]), encoding="utf-8")
    reference_r_path.write_text(str(template["reference_r"]), encoding="utf-8")
    analysis_py_path.write_text(str(template["analysis_py"]), encoding="utf-8")
    checks_path.write_text(
        json.dumps(template["checks"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    human_report_path.write_text(str(template["human_report"]), encoding="utf-8")
    provenance_path.write_text(
        json.dumps(template["provenance"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    results_readme_path.write_text(str(template["results_readme"]), encoding="utf-8")
    results_manifest_path.write_text(
        json.dumps(template["results_manifest"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    report_path.write_text(
        json.dumps(template["report_payload"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "manifest": manifest_path,
        "claims": claims_path,
        "readme": readme_path,
        "reference_r": reference_r_path,
        "analysis_py": analysis_py_path,
        "checks": checks_path,
        "human_report": human_report_path,
        "provenance": provenance_path,
        "results_readme": results_readme_path,
        "results_manifest": results_manifest_path,
        "report": report_path,
    }
