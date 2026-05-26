from __future__ import annotations

from datetime import date
from pathlib import Path

from .definitions import (
    BUNDLE_DEFINITIONS,
    CLAIM_DEFINITIONS,
    PCM2_REFERENCE_SCRIPT_PATH,
    PCM2_SOURCE_LOCATOR,
    STUDY_ID,
    STUDY_ONE_REFERENCE_ROOT,
    SUMMARY_EVIDENCE_ID,
)
from .parity import (
    build_primate_pgls_signal_scalar_parity_table,
    render_primate_pgls_signal_scalar_parity_table_markdown,
)
from .runtime import (
    ancestral_reconstruction_payloads,
    ancestral_reconstruction_reports,
    baseline_gls_payload,
    baseline_gls_report,
    continuous_mode_fit_payloads,
    continuous_mode_fit_reports,
    coverage_boundary_payload,
    estimated_lambda_pgls_payload,
    estimated_lambda_pgls_report,
    fixed_reference_lambda_pgls_payload,
    fixed_reference_lambda_pgls_report,
    likelihood_ratio_test_payloads,
    load_r_reference_results,
    mode_comparison_report,
    ordered_trait_values,
    signal_report,
    signal_test_payload,
    source_reference_paths,
    transformed_tree_reports_for_repo,
    tree_rescaling_payloads,
)

__all__ = [
    "build_primate_pgls_signal_bundle",
    "build_primate_pgls_signal_bundles",
]


def _report_payload_for_bundle(repo_root: Path, evidence_id: str) -> dict[str, object]:
    r_results = load_r_reference_results(repo_root)
    if evidence_id == "evidence-001":
        scalar_table = build_primate_pgls_signal_scalar_parity_table(repo_root)
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "object_names": r_results["source_contract"]["object_names"],
            "row_count": r_results["source_contract"]["row_count"],
            "tip_count": r_results["source_contract"]["tip_count"],
            "species_tip_match": r_results["source_contract"]["species_tip_match"],
            "governed_reload_inputs": [
                (STUDY_ONE_REFERENCE_ROOT / "reference_primate.csv").as_posix(),
                (
                    STUDY_ONE_REFERENCE_ROOT / "reference_trimmed_primatetree.nwk"
                ).as_posix(),
            ],
            "scalar_row_count": scalar_table["row_count"],
            "verdict_counts": scalar_table["verdict_counts"],
        }
    if evidence_id == "evidence-002":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_baseline": r_results["baseline_gls"],
            "bijux_baseline": baseline_gls_payload(baseline_gls_report(repo_root)),
            "r_fixed_lambda_equivalence": r_results[
                "fixed_lambda_gls_matches_baseline"
            ],
        }
    if evidence_id == "evidence-003":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_fixed_reference_lambda": r_results["fixed_reference_lambda_pgls"],
            "bijux_fixed_reference_lambda": fixed_reference_lambda_pgls_payload(
                fixed_reference_lambda_pgls_report(repo_root)
            ),
            "r_estimated_lambda": r_results["estimated_lambda_pgls"],
            "bijux_estimated_lambda": estimated_lambda_pgls_payload(
                estimated_lambda_pgls_report(repo_root)
            ),
        }
    if evidence_id == "evidence-004":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_signal_test": r_results["signal_test"],
            "bijux_signal_test": signal_test_payload(signal_report(repo_root)),
        }
    if evidence_id == "evidence-005":
        baseline = baseline_gls_payload(baseline_gls_report(repo_root))
        estimated = estimated_lambda_pgls_payload(
            estimated_lambda_pgls_report(repo_root)
        )
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "baseline_diagnostics": {
                "r": r_results["baseline_gls"]["diagnostics"],
                "bijux": baseline["diagnostics"],
            },
            "estimated_lambda_diagnostics": {
                "r": r_results["estimated_lambda_pgls"]["diagnostics"],
                "bijux": estimated["diagnostics"],
            },
        }
    if evidence_id == "evidence-006":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_tree_rescaling": r_results["tree_rescaling"],
            "bijux_tree_rescaling": tree_rescaling_payloads(
                transformed_tree_reports_for_repo(repo_root)
            ),
        }
    if evidence_id == "evidence-007":
        brownian_fit, ou_fit, early_burst_fit = continuous_mode_fit_reports(repo_root)
        tip_values = ordered_trait_values(
            source_reference_paths(repo_root)[1],
            brownian_fit.taxa,
            trait="longevity",
            taxon_column="species",
        )
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_continuous_mode_fits": r_results["continuous_mode_fits"],
            "bijux_continuous_mode_fits": continuous_mode_fit_payloads(
                brownian_fit,
                ou_fit,
                early_burst_fit,
                tip_values=tip_values,
            ),
        }
    if evidence_id == "evidence-008":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_likelihood_ratio_tests": r_results["likelihood_ratio_tests"],
            "bijux_likelihood_ratio_tests": likelihood_ratio_test_payloads(
                mode_comparison_report(repo_root)
            ),
        }
    if evidence_id == "evidence-009":
        brownian_ancestral, early_burst_ancestral = ancestral_reconstruction_reports(
            repo_root
        )
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_ancestral_reconstruction": r_results["ancestral_reconstruction"],
            "bijux_ancestral_reconstruction": ancestral_reconstruction_payloads(
                brownian_ancestral,
                early_burst_ancestral,
            ),
        }
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": evidence_id,
        "coverage_boundaries": coverage_boundary_payload(),
    }


def build_primate_pgls_signal_bundle(
    repo_root: Path,
    evidence_id: str,
) -> dict[str, object]:
    definition = next(
        (
            definition
            for definition in BUNDLE_DEFINITIONS
            if definition["evidence_id"] == evidence_id
        ),
        None,
    )
    if definition is None:
        raise KeyError(evidence_id)
    report_payload = _report_payload_for_bundle(repo_root, evidence_id)
    manifest = _manifest_for_bundle(repo_root, definition, report_payload)
    bundle = {
        "manifest": manifest,
        "claims": _claims_payload(definition),
        "report_payload": report_payload,
        "report_filename": definition["report_filename"],
        "readme": _readme_for_bundle(definition),
    }
    if evidence_id == SUMMARY_EVIDENCE_ID:
        scalar_table = build_primate_pgls_signal_scalar_parity_table(repo_root)
        bundle["scalar_parity_table"] = scalar_table
        bundle["scalar_parity_markdown"] = (
            render_primate_pgls_signal_scalar_parity_table_markdown(scalar_table)
        )
    return bundle


def _manifest_for_bundle(
    repo_root: Path,
    definition: dict[str, object],
    report_payload: dict[str, object],
) -> dict[str, object]:
    evidence_id = str(definition["evidence_id"])
    source_basis = [
        {
            "kind": "external-source-descriptor",
            "label": "Lund source descriptors",
            "locator": f"evidence-book/studies/{STUDY_ID}/provenance/lund-course-sources.json",
        },
        {
            "kind": "repository-reference",
            "label": "governed primate reference table",
            "locator": (STUDY_ONE_REFERENCE_ROOT / "reference_primate.csv").as_posix(),
        },
        {
            "kind": "repository-reference",
            "label": "governed primate trimmed tree",
            "locator": (
                STUDY_ONE_REFERENCE_ROOT / "reference_trimmed_primatetree.nwk"
            ).as_posix(),
        },
        {
            "kind": "repository-reference",
            "label": "R reference results for the PGLS and signal study",
            "locator": f"evidence-book/studies/{STUDY_ID}/reference/reference_results.json",
        },
        {
            "kind": "repository-reference",
            "label": f"{definition['title']} report payload",
            "locator": f"evidence-book/studies/{STUDY_ID}/{evidence_id}/{definition['report_filename']}",
        },
    ]
    if evidence_id == SUMMARY_EVIDENCE_ID:
        source_basis.append(
            {
                "kind": "repository-reference",
                "label": "scalar parity table",
                "locator": f"evidence-book/studies/{STUDY_ID}/{SUMMARY_EVIDENCE_ID}/scalar-parity-table.json",
            }
        )
    claim = CLAIM_DEFINITIONS[str(definition["claim_id"])]
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": evidence_id,
        "evidence_title": definition["title"],
        "summary": definition["summary"],
        "owner_package": "bijux-phylogenetics",
        "claim_ids": [definition["claim_id"]],
        "source_basis": source_basis,
        "freshness": {
            "last_generated_on": date.today().isoformat(),
            "governed_code_paths": [
                "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/studies/primate_pgls_and_signal"
            ],
            "source_basis_locators": [entry["locator"] for entry in source_basis],
        },
        "ownership": {
            "owner_package": "bijux-phylogenetics",
            "analytical_surfaces": definition["analytical_surfaces"],
        },
        "claim_tags": definition["claim_tags"],
        "comparison_mode": definition["comparison_mode"],
        "verdict": {
            "status": claim["verdict"],
            "summary": claim["summary"],
        },
        "limitations": definition["limitations"],
        "source_fragments": definition["source_fragments"],
        "reference_script_locators": [f"{PCM2_REFERENCE_SCRIPT_PATH}#L1-L200"],
        "supporting_report_locator": (
            f"evidence-book/studies/{STUDY_ID}/{evidence_id}/{definition['report_filename']}"
        ),
        "report_keys": sorted(report_payload.keys()),
    }


def _claims_payload(definition: dict[str, object]) -> dict[str, object]:
    claim = CLAIM_DEFINITIONS[str(definition["claim_id"])]
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": definition["evidence_id"],
        "claim_count": 1,
        "claims": [
            {
                "claim_id": definition["claim_id"],
                "claim_title": claim["claim_title"],
                "summary": claim["summary"],
                "verdict": claim["verdict"],
                "evidence_ids": [definition["evidence_id"]],
                "source_fragments": definition["source_fragments"],
            }
        ],
    }


def _readme_for_bundle(definition: dict[str, object]) -> str:
    lines = [
        f"# {definition['title']}",
        "",
        definition["summary"],
        "",
        f"- evidence id: `{definition['evidence_id']}`",
        f"- source fragments: {', '.join(f'`{fragment}`' for fragment in definition['source_fragments'])}",
        "",
        "## Limitations",
        "",
    ]
    for limitation in definition["limitations"]:
        lines.append(f"- {limitation}")
    lines.extend(
        [
            "",
            "## Source Locators",
            "",
            f"- `{PCM2_SOURCE_LOCATOR}`",
            f"- `{PCM2_REFERENCE_SCRIPT_PATH}`",
            "",
        ]
    )
    return "\n".join(lines)


def build_primate_pgls_signal_bundles(repo_root: Path) -> dict[str, dict[str, object]]:
    bundles: dict[str, dict[str, object]] = {}
    for definition in BUNDLE_DEFINITIONS:
        bundles[definition["evidence_id"]] = build_primate_pgls_signal_bundle(
            repo_root,
            definition["evidence_id"],
        )
    return bundles
