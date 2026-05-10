from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path


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


def build_evidence_bundle_template(spec: EvidenceBundleTemplateSpec) -> dict[str, object]:
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
    return {
        "manifest": manifest,
        "claims": claims,
        "readme": readme,
        "report_filename": "template-report.json",
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
    report_path = output_root / str(template["report_filename"])
    manifest_path.write_text(
        json.dumps(template["manifest"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    claims_path.write_text(
        json.dumps(template["claims"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    readme_path.write_text(str(template["readme"]), encoding="utf-8")
    report_path.write_text(
        json.dumps(template["report_payload"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "manifest": manifest_path,
        "claims": claims_path,
        "readme": readme_path,
        "report": report_path,
    }
