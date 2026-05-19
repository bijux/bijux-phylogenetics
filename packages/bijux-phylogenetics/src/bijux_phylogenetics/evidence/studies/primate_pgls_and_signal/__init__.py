from __future__ import annotations

from datetime import date
from functools import lru_cache
from pathlib import Path

from .definitions import (
    BUNDLE_DEFINITIONS,
    CLAIM_DEFINITIONS,
    FAMILY_DEFINITIONS,
    FRAGMENT_DEFINITIONS,
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
    ancestral_reconstruction_payloads as _ancestral_reconstruction_payloads,
    ancestral_reconstruction_reports as _ancestral_reconstruction_reports,
    baseline_gls_payload as _baseline_gls_payload,
    baseline_gls_report as _baseline_gls_report,
    continuous_mode_fit_payloads as _continuous_mode_fit_payloads,
    continuous_mode_fit_reports as _continuous_mode_fit_reports,
    coverage_boundary_payload as _coverage_boundary_payload,
    estimated_lambda_pgls_payload as _estimated_lambda_pgls_payload,
    estimated_lambda_pgls_report as _estimated_lambda_pgls_report,
    likelihood_ratio_test_payloads as _likelihood_ratio_test_payloads,
    load_r_reference_results as _load_r_reference_results,
    mode_comparison_report as _mode_comparison_report,
    ordered_trait_values as _ordered_trait_values,
    signal_report as _signal_report,
    signal_test_payload as _signal_test_payload,
    source_reference_paths as _source_reference_paths,
    transformed_tree_reports_for_repo as _transformed_tree_reports,
    tree_rescaling_payloads as _tree_rescaling_payloads,
)



def _line_spec_to_locators(spec: str) -> list[str]:
    locators: list[str] = []
    for part in spec.split(","):
        normalized = part.strip()
        if not normalized:
            continue
        if "-" in normalized:
            start_text, end_text = normalized.split("-", maxsplit=1)
            locators.append(
                f"{PCM2_SOURCE_LOCATOR}#L{int(start_text)}-L{int(end_text)}"
            )
        else:
            locators.append(f"{PCM2_SOURCE_LOCATOR}#L{int(normalized)}")
    return locators


def _line_spec_to_spans(spec: str) -> list[dict[str, int]]:
    spans: list[dict[str, int]] = []
    for part in spec.split(","):
        normalized = part.strip()
        if not normalized:
            continue
        if "-" in normalized:
            start_text, end_text = normalized.split("-", maxsplit=1)
            spans.append({"start_line": int(start_text), "end_line": int(end_text)})
        else:
            line = int(normalized)
            spans.append({"start_line": line, "end_line": line})
    return spans


@lru_cache(maxsize=1)
def build_primate_pgls_signal_external_sources() -> dict[str, object]:
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "intake_policy": "read-only-external-source",
        "source_count": 3,
        "sources": [
            {
                "source_id": "lund-pcm2-script",
                "kind": "external-course-script",
                "label": "Lund PCM2 modes and PGLS lecture script",
                "locator": PCM2_SOURCE_LOCATOR,
                "path_hint": "PCM2_modes_pgls/Lecture/R/Scripts/PCM2_modes_pgls.R",
                "concept_tags": [
                    "pgls",
                    "phylogenetic-signal",
                    "gls",
                    "evolutionary-models",
                ],
            },
            {
                "source_id": "lund-primate-rdata",
                "kind": "external-course-workspace",
                "label": "Lund primate workspace RData",
                "locator": "external:lund/pcm2-modes-pgls/data/primate.RData",
                "provides": ["primate", "primatetree"],
            },
            {
                "source_id": "governed-primate-reference-artifacts",
                "kind": "repository-reference",
                "label": "Governed primate CSV and trimmed tree from the earlier evidence study",
                "locator": "evidence-book/studies/primate-longevity-signal/evidence-001",
                "provides": [
                    "reference_primate.csv",
                    "reference_trimmed_primatetree.nwk",
                ],
            },
        ],
    }


@lru_cache(maxsize=1)
def build_primate_pgls_signal_source_fragment_map() -> dict[str, object]:
    fragments = []
    for definition in FRAGMENT_DEFINITIONS:
        fragments.append(
            {
                "fragment_id": definition["fragment_id"],
                "fragment_title": definition["fragment_title"],
                "concept_family": definition["family_id"],
                "claim_ids": definition["claim_ids"],
                "evidence_id": definition["evidence_id"],
                "supporting_evidence_ids": definition["supporting_evidence_ids"],
                "script_line_spec": definition["script_line_spec"],
                "script_line_spans": _line_spec_to_spans(
                    definition["script_line_spec"]
                ),
                "script_locators": _line_spec_to_locators(
                    definition["script_line_spec"]
                ),
                "parity_expectation": definition["parity_expectation"],
                "comparison_kind": definition["comparison_kind"],
                "block_status": definition["block_status"],
                "review_note": definition["review_note"],
                "scope": definition["scope"],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": SUMMARY_EVIDENCE_ID,
        "fragment_count": len(fragments),
        "fragments": fragments,
    }


@lru_cache(maxsize=1)
def build_primate_pgls_signal_parity_policy() -> dict[str, object]:
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "policy_count": 10,
        "evidence_id": SUMMARY_EVIDENCE_ID,
        "policies": [
            {
                "family_id": "workflow-contracts",
                "family_title": FAMILY_DEFINITIONS["workflow-contracts"]["title"],
                "parity_expectation": "exact",
                "comparison_kind": "exact_answer",
                "metric_tolerances": [
                    {"metric_kind": "default", "tolerance_abs_diff": 0.0}
                ],
                "rule": "Object names, row counts, tip counts, and repository locators must match exactly.",
                "source_fragments": ["workspace-reload-contract"],
            },
            {
                "family_id": "transformed-tree-workflows",
                "family_title": FAMILY_DEFINITIONS["transformed-tree-workflows"][
                    "title"
                ],
                "parity_expectation": "exact",
                "comparison_kind": "exact_answer",
                "metric_tolerances": [
                    {"metric_kind": "branch_count", "tolerance_abs_diff": 0.0},
                    {"metric_kind": "total_branch_length", "tolerance_abs_diff": 0.0},
                ],
                "rule": "Rounded transformed branch counts and total branch lengths must match exactly for the governed tree-rescaling checkpoints.",
                "source_fragments": ["transformed-tree-workflows"],
            },
            {
                "family_id": "continuous-model-fitting",
                "family_title": FAMILY_DEFINITIONS["continuous-model-fitting"]["title"],
                "parity_expectation": "statistical_tolerance",
                "comparison_kind": "tolerance_or_equivalence",
                "metric_tolerances": [
                    {"metric_kind": "parameter_value", "tolerance_abs_diff": 0.05},
                    {"metric_kind": "root_state", "tolerance_abs_diff": 0.5},
                    {"metric_kind": "rate", "tolerance_abs_diff": 25000.0},
                    {"metric_kind": "log_likelihood", "tolerance_abs_diff": 0.25},
                    {"metric_kind": "aic", "tolerance_abs_diff": 0.5},
                ],
                "rule": "Brownian, OU, and early-burst intercept fits may drift numerically, but parameter ranking and fit-quality conclusions must remain aligned.",
                "source_fragments": ["continuous-model-comparison"],
            },
            {
                "family_id": "likelihood-ratio-tests",
                "family_title": FAMILY_DEFINITIONS["likelihood-ratio-tests"]["title"],
                "parity_expectation": "statistical_tolerance",
                "comparison_kind": "tolerance_or_equivalence",
                "metric_tolerances": [
                    {"metric_kind": "statistic", "tolerance_abs_diff": 0.05},
                    {"metric_kind": "p_value", "tolerance_abs_diff": 0.001},
                ],
                "rule": "Likelihood-ratio statistics may drift slightly, but the same model-comparison decisions must hold.",
                "source_fragments": [
                    "continuous-model-comparison",
                    "evolutionary-mode-likelihood-ratios",
                ],
            },
            {
                "family_id": "baseline-regression",
                "family_title": FAMILY_DEFINITIONS["baseline-regression"]["title"],
                "parity_expectation": "exact",
                "comparison_kind": "exact_answer",
                "metric_tolerances": [
                    {"metric_kind": "coefficient", "tolerance_abs_diff": 1e-06},
                    {"metric_kind": "log_likelihood", "tolerance_abs_diff": 1e-06},
                    {"metric_kind": "r_squared", "tolerance_abs_diff": 1e-06},
                ],
                "rule": "The baseline regression is expected to agree numerically to near machine precision.",
                "source_fragments": ["baseline-gls-fit"],
            },
            {
                "family_id": "phylogenetic-regression",
                "family_title": FAMILY_DEFINITIONS["phylogenetic-regression"]["title"],
                "parity_expectation": "statistical_tolerance",
                "comparison_kind": "tolerance_or_equivalence",
                "metric_tolerances": [
                    {"metric_kind": "lambda_value", "tolerance_abs_diff": 0.05},
                    {"metric_kind": "coefficient", "tolerance_abs_diff": 1.0},
                    {"metric_kind": "slope", "tolerance_abs_diff": 0.1},
                    {"metric_kind": "log_likelihood", "tolerance_abs_diff": 0.25},
                ],
                "rule": "Likelihood and parameter estimates may drift modestly, but sign, direction, and coefficient significance decisions must stay aligned.",
                "source_fragments": ["pagel-lambda-regression"],
            },
            {
                "family_id": "phylogenetic-signal",
                "family_title": FAMILY_DEFINITIONS["phylogenetic-signal"]["title"],
                "parity_expectation": "statistical_tolerance",
                "comparison_kind": "tolerance_or_equivalence",
                "metric_tolerances": [
                    {"metric_kind": "lambda_value", "tolerance_abs_diff": 0.05},
                    {"metric_kind": "likelihood_ratio", "tolerance_abs_diff": 0.05},
                    {"metric_kind": "p_value", "tolerance_abs_diff": 0.001},
                ],
                "rule": "Signal testing must preserve the same lambda-zero rejection decision while allowing bounded likelihood drift.",
                "source_fragments": ["phylogenetic-signal-test"],
            },
            {
                "family_id": "diagnostics",
                "family_title": FAMILY_DEFINITIONS["diagnostics"]["title"],
                "parity_expectation": "scientific_equivalence",
                "comparison_kind": "scientific_equivalence",
                "metric_tolerances": [
                    {"metric_kind": "qq_correlation", "tolerance_abs_diff": 0.02},
                    {
                        "metric_kind": "abs_residual_fitted_correlation",
                        "tolerance_abs_diff": 0.05,
                    },
                    {"metric_kind": "outlier_count", "tolerance_abs_diff": 1.0},
                ],
                "rule": "Diagnostics may differ numerically, but they must sustain the same high-level review conclusion about skew, heteroscedasticity, and outlier severity.",
                "source_fragments": [
                    "baseline-gls-diagnostics",
                    "estimated-lambda-diagnostics",
                ],
            },
            {
                "family_id": "ancestral-reconstruction",
                "family_title": FAMILY_DEFINITIONS["ancestral-reconstruction"]["title"],
                "parity_expectation": "statistical_tolerance",
                "comparison_kind": "tolerance_or_equivalence",
                "metric_tolerances": [
                    {"metric_kind": "node_count", "tolerance_abs_diff": 0.0},
                    {"metric_kind": "estimate", "tolerance_abs_diff": 0.5},
                ],
                "rule": "Brownian and early-burst ancestral estimates may drift slightly, but the same node-level trajectory and teaching conclusion must remain intact.",
                "source_fragments": ["ancestral-mode-comparison"],
            },
            {
                "family_id": "coverage-boundaries",
                "family_title": FAMILY_DEFINITIONS["coverage-boundaries"]["title"],
                "parity_expectation": "not_comparable",
                "comparison_kind": "not_comparable",
                "metric_tolerances": [],
                "rule": "These fragments are tracked explicitly as open trust boundaries and therefore cannot be promoted to parity claims yet.",
                "source_fragments": ["mode-linked-intercept-models"],
            },
        ],
    }


def build_primate_pgls_signal_claim_registry(repo_root: Path) -> dict[str, object]:
    bundles = build_primate_pgls_signal_bundles(repo_root)
    claims = []
    for definition in BUNDLE_DEFINITIONS:
        claim = CLAIM_DEFINITIONS[definition["claim_id"]]
        claims.append(
            {
                "claim_id": definition["claim_id"],
                "claim_title": claim["claim_title"],
                "summary": claim["summary"],
                "verdict": claim["verdict"],
                "evidence_ids": [definition["evidence_id"]],
                "source_fragments": definition["source_fragments"],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "claim_count": len(claims),
        "claims": claims,
        "bundle_count": len(bundles),
    }


def build_primate_pgls_signal_family_index(repo_root: Path) -> dict[str, object]:
    fragments = build_primate_pgls_signal_source_fragment_map()["fragments"]
    families = []
    for family_id, family in FAMILY_DEFINITIONS.items():
        matching_fragments = [
            fragment
            for fragment in fragments
            if fragment["concept_family"] == family_id
        ]
        claim_ids = sorted(
            {
                claim_id
                for fragment in matching_fragments
                for claim_id in fragment["claim_ids"]
            }
        )
        evidence_ids = sorted(
            {
                fragment["evidence_id"]
                for fragment in matching_fragments
                if fragment["evidence_id"]
            }
        )
        claims = {
            claim_id: CLAIM_DEFINITIONS[claim_id]["claim_title"]
            for claim_id in claim_ids
        }
        verdicts = {CLAIM_DEFINITIONS[claim_id]["verdict"] for claim_id in claim_ids}
        if verdicts == {"matched"}:
            family_verdict = "matched"
        elif "matched_with_tolerance" in verdicts and not (
            "mismatch_unexplained" in verdicts or "not_comparable" in verdicts
        ):
            family_verdict = "matched_with_tolerance"
        elif verdicts == {"not_comparable"}:
            family_verdict = "not_comparable"
        else:
            family_verdict = (
                "not_comparable"
                if "not_comparable" in verdicts
                else "mismatch_unexplained"
            )
        families.append(
            {
                "family_id": family_id,
                "family_title": family["title"],
                "summary": family["summary"],
                "fragment_count": len(matching_fragments),
                "fragment_ids": [
                    fragment["fragment_id"] for fragment in matching_fragments
                ],
                "claim_ids": claim_ids,
                "claim_titles": claims,
                "evidence_ids": evidence_ids,
                "family_verdict": family_verdict,
                "coverage_status": (
                    "coverage-gap" if family_id == "coverage-boundaries" else "covered"
                ),
                "known_gaps": []
                if family_id != "coverage-boundaries"
                else [
                    "The lecture corBlomberg intercept-mode likelihood sweep is still an explicit coverage boundary.",
                ],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": SUMMARY_EVIDENCE_ID,
        "source_family_id": "pcm2-modes-pgls",
        "source_family_title": "PCM2 modes and PGLS evidence family",
        "family_count": len(families),
        "families": families,
    }


def _report_payload_for_bundle(repo_root: Path, evidence_id: str) -> dict[str, object]:
    r_results = _load_r_reference_results(repo_root)
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
            "bijux_baseline": _baseline_gls_payload(_baseline_gls_report(repo_root)),
            "r_fixed_lambda_equivalence": r_results[
                "fixed_lambda_gls_matches_baseline"
            ],
        }
    if evidence_id == "evidence-003":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_estimated_lambda": r_results["estimated_lambda_pgls"],
            "bijux_estimated_lambda": _estimated_lambda_pgls_payload(
                _estimated_lambda_pgls_report(repo_root)
            ),
        }
    if evidence_id == "evidence-004":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_signal_test": r_results["signal_test"],
            "bijux_signal_test": _signal_test_payload(_signal_report(repo_root)),
        }
    if evidence_id == "evidence-005":
        baseline = _baseline_gls_payload(_baseline_gls_report(repo_root))
        estimated = _estimated_lambda_pgls_payload(
            _estimated_lambda_pgls_report(repo_root)
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
            "bijux_tree_rescaling": _tree_rescaling_payloads(
                _transformed_tree_reports(repo_root)
            ),
        }
    if evidence_id == "evidence-007":
        brownian_fit, ou_fit, early_burst_fit = _continuous_mode_fit_reports(repo_root)
        tip_values = _ordered_trait_values(
            _source_reference_paths(repo_root)[1],
            brownian_fit.taxa,
            trait="longevity",
            taxon_column="species",
        )
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_continuous_mode_fits": r_results["continuous_mode_fits"],
            "bijux_continuous_mode_fits": _continuous_mode_fit_payloads(
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
            "bijux_likelihood_ratio_tests": _likelihood_ratio_test_payloads(
                _mode_comparison_report(repo_root)
            ),
        }
    if evidence_id == "evidence-009":
        brownian_ancestral, early_burst_ancestral = _ancestral_reconstruction_reports(
            repo_root
        )
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_ancestral_reconstruction": r_results["ancestral_reconstruction"],
            "bijux_ancestral_reconstruction": _ancestral_reconstruction_payloads(
                brownian_ancestral,
                early_burst_ancestral,
            ),
        }
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": evidence_id,
        "coverage_boundaries": _coverage_boundary_payload(),
    }


def build_primate_pgls_signal_bundle(
    repo_root: Path, evidence_id: str
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


def build_primate_pgls_signal_evidence_registry(
    repo_root: Path,
) -> dict[str, object]:
    bundles = build_primate_pgls_signal_bundles(repo_root)
    evidences = []
    for definition in BUNDLE_DEFINITIONS:
        claim = CLAIM_DEFINITIONS[str(definition["claim_id"])]
        evidences.append(
            {
                "evidence_id": definition["evidence_id"],
                "title": definition["title"],
                "coverage_status": (
                    "coverage-gap"
                    if claim["verdict"] == "not_comparable"
                    else "covered"
                ),
                "claim_id": definition["claim_id"],
                "verdict": claim["verdict"],
                "analytical_surfaces": definition["analytical_surfaces"],
                "source_fragments": definition["source_fragments"],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "bundle_count": len(bundles),
        "evidence_count": len(evidences),
        "coverage_boundary_evidence_ids": [
            entry["evidence_id"]
            for entry in evidences
            if entry["coverage_status"] == "coverage-gap"
        ],
        "evidences": evidences,
    }


def render_primate_pgls_signal_study_manifest(repo_root: Path) -> dict[str, object]:
    registry = build_primate_pgls_signal_evidence_registry(repo_root)
    return {
        "study_id": STUDY_ID,
        "study_title": "Primate PGLS and signal evidence study",
        "summary": "Governed parity study for the regression, transformed-tree, evolutionary-mode fit, and ancestral sections of the Lund primate comparative lecture, with an explicit remaining intercept-mode boundary.",
        "owner_package": "bijux-phylogenetics",
        "study_categories": ["teaching-study", "migration-study"],
        "confidence_posture": "governed-parity-in-progress",
        "coverage_boundary_evidence_ids": registry["coverage_boundary_evidence_ids"],
        "evidence_registry_locator": (
            f"evidence-book/studies/{STUDY_ID}/evidence-registry.json"
        ),
        "study_scope": {
            "coverage_focus": [
                "rdata-reload",
                "baseline-regression",
                "phylogenetic-regression",
                "phylogenetic-signal",
                "transformed-tree-workflows",
                "continuous-model-fitting",
                "likelihood-ratio-tests",
                "ancestral-reconstruction",
            ],
            "untouched_source_locators": [
                PCM2_SOURCE_LOCATOR,
                "external:lund/pcm2-modes-pgls/data/primate.RData",
            ],
        },
    }


def render_primate_pgls_signal_study_readme(repo_root: Path) -> str:
    registry = build_primate_pgls_signal_evidence_registry(repo_root)
    lines = [
        "# Primate PGLS And Signal",
        "",
        "This study turns the regression, transformed-tree, evolutionary-mode fit,",
        "likelihood-ratio, and ancestral sections of the Lund primate comparative",
        "lecture into governed Evidence IDs backed by checked-in R reference outputs",
        "and canonical `bijux-phylogenetics` reproductions.",
        "",
        "It is intentionally strict about confidence posture:",
        "",
        "- baseline GLS, Pagel-lambda PGLS, signal testing, transformed-tree",
        "  workflows, fitContinuous-style mode comparisons, and ancestral-mode",
        "  reconstructions are backed by governed parity bundles",
        "- the lecture corBlomberg intercept sweep remains visible as an explicit",
        "  coverage boundary instead of being implied as validated",
        "",
        "Current bundles:",
        "",
    ]
    for entry in registry["evidences"]:
        title = str(entry["title"]).removeprefix("Primate ").removesuffix(" bundle")
        lines.append(f"- `{entry['evidence_id']}` {title}")
    lines.append("")
    return "\n".join(lines)


def _manifest_for_bundle(
    repo_root: Path, definition: dict[str, object], report_payload: dict[str, object]
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
                "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/studies/primate_pgls_and_signal.py"
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
            repo_root, definition["evidence_id"]
        )
    return bundles
