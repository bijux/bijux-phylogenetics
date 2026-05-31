from __future__ import annotations

from datetime import date
from pathlib import Path

from .definitions import COMPONENT_BUNDLE_DEFINITIONS, REFERENCE_BUNDLE_ROOT, STUDY_ID
from .report_payloads import build_component_report_payload
from .study_context import load_reference_context, reference_script_locators


def build_primate_pcm1_component_bundles(
    repo_root: Path,
) -> dict[str, dict[str, object]]:
    context = load_reference_context(repo_root)
    bundles: dict[str, dict[str, object]] = {}
    for spec in COMPONENT_BUNDLE_DEFINITIONS:
        report_payload = build_component_report_payload(spec, context)
        manifest = build_bundle_manifest(
            spec,
            report_filename=spec["report_filename"],
            report_payload=report_payload,
        )
        bundles[spec["evidence_id"]] = {
            "manifest": manifest,
            "claims": build_claims_payload(spec),
            "report_filename": spec["report_filename"],
            "report_payload": report_payload,
            "readme": render_bundle_readme(
                spec,
                report_filename=spec["report_filename"],
                manifest=manifest,
            ),
        }
    return bundles


def build_claims_payload(spec: dict[str, object]) -> dict[str, object]:
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": spec["evidence_id"],
        "claim_count": 1,
        "claims": [
            {
                "claim_id": spec["claim_id"],
                "claim_title": spec["claim_title"],
                "summary": spec["claim_summary"],
                "verdict": spec["verdict"],
                "evidence_ids": [spec["evidence_id"]],
                "source_fragments": spec["source_fragments"],
            }
        ],
    }


def build_bundle_manifest(
    spec: dict[str, object], *, report_filename: str, report_payload: dict[str, object]
) -> dict[str, object]:
    source_basis = list(spec["source_basis"]) + [
        {
            "kind": "repository-reference",
            "label": f"{spec['title']} report payload",
            "locator": (
                f"evidence-book/studies/primate-longevity-signal/{spec['evidence_id']}/"
                f"{report_filename}"
            ),
        }
    ]
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": spec["evidence_id"],
        "evidence_title": spec["title"],
        "summary": spec["summary"],
        "owner_package": "bijux-phylogenetics",
        "claim_ids": [spec["claim_id"]],
        "source_basis": source_basis,
        "freshness": {
            "last_generated_on": date.today().isoformat(),
            "governed_code_paths": [
                "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/studies/primate_pcm1_component_bundles/definitions.py",
                "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/studies/primate_pcm1_component_bundles/study_context.py",
                "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/studies/primate_pcm1_component_bundles/report_payloads.py",
                "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/studies/primate_pcm1_component_bundles/bundle_builder.py",
            ],
            "source_basis_locators": [entry["locator"] for entry in source_basis],
        },
        "ownership": {
            "owner_package": "bijux-phylogenetics",
            "analytical_surfaces": spec["analytical_surfaces"],
        },
        "claim_tags": spec["claim_tags"],
        "comparison_mode": spec["comparison_mode"],
        "verdict": {
            "status": spec["verdict"],
            "summary": spec["claim_summary"],
        },
        "limitations": spec["limitations"],
        "source_fragments": spec["source_fragments"],
        "reference_script_locators": reference_script_locators(
            spec["reference_line_specs"]
        ),
        "reference_bundle_locator": REFERENCE_BUNDLE_ROOT.as_posix(),
        "supporting_report_locator": (
            f"evidence-book/studies/primate-longevity-signal/{spec['evidence_id']}/"
            f"{report_filename}"
        ),
        "report_keys": sorted(report_payload.keys()),
    }


def render_bundle_readme(
    spec: dict[str, object], *, report_filename: str, manifest: dict[str, object]
) -> str:
    lines = [
        f"# {spec['title']}",
        "",
        spec["summary"],
        "",
        f"- evidence id: `{spec['evidence_id']}`",
        f"- claim id: `{spec['claim_id']}`",
        f"- verdict: `{spec['verdict']}`",
        f"- source fragments: `{', '.join(spec['source_fragments'])}`",
        "",
        "Governed files:",
        "",
        f"- `{report_filename}`",
        "- `claims.json`",
        "- `manifest.json`",
        "",
        "Reference script locators:",
        "",
    ]
    lines.extend(f"- `{locator}`" for locator in manifest["reference_script_locators"])
    lines.extend(
        [
            "",
            "Limitations:",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in spec["limitations"])
    lines.append("")
    return "\n".join(lines)
