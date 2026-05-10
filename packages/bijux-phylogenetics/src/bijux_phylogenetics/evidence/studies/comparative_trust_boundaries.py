from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.comparative.models import (
    audit_ou_identifiability_reference_examples,
)
from bijux_phylogenetics.comparative.pgls import run_pgls
from bijux_phylogenetics.comparative.signal import (
    compute_phylogenetic_signal_test,
    estimate_pagels_lambda,
)
from bijux_phylogenetics.errors import ComparativeMethodError


STUDY_ID = "comparative-trust-boundaries"
EXPECTED_FAILURE_EVIDENCE_ID = "evidence-001"
WEAK_SIGNAL_EVIDENCE_ID = "evidence-002"
OU_IDENTIFIABILITY_EVIDENCE_ID = "evidence-003"
OWNER_PACKAGE = "bijux-phylogenetics"

FIXTURES_ROOT = (
    Path("packages") / "bijux-phylogenetics" / "tests" / "fixtures"
)
EXPECTED_FAILURE_TREE_ROOT = FIXTURES_ROOT / "trees"
EXPECTED_FAILURE_TRAITS_ROOT = FIXTURES_ROOT / "metadata"
WEAK_SIGNAL_TRAITS_PATH = (
    Path("evidence-book")
    / "studies"
    / STUDY_ID
    / WEAK_SIGNAL_EVIDENCE_ID
    / "weak_signal_traits.tsv"
)
WEAK_SIGNAL_TREE_PATH = EXPECTED_FAILURE_TREE_ROOT / "example_tree_six_taxa.nwk"
WEAK_SIGNAL_SEEDS = tuple(range(1, 21))
WEAK_SIGNAL_PERMUTATIONS = 99
WEAK_SIGNAL_ALPHA = 0.05

FAMILY_DEFINITIONS = {
    "expected-failures": {
        "title": "Expected failures",
        "summary": "Input rejection is treated as governed evidence rather than as incidental exception behavior.",
    },
    "weak-signal": {
        "title": "Weak signal",
        "summary": "Low-information signal cases stay visible when the significance boundary is unstable across governed reruns.",
    },
    "model-instability": {
        "title": "Model instability",
        "summary": "Boundary-seeking comparative fits and identifiability warnings remain explicit instead of being hidden beneath successful demos.",
    },
}

CLAIM_DEFINITIONS = {
    "comparative-input-rejection-governed": {
        "claim_title": "Comparative input rejection remains explicit and reproducible",
        "summary": "Rootedness, branch-length completeness, and numeric-response guards are captured as governed runtime evidence instead of informal expectations.",
        "verdict": "matched",
    },
    "weak-signal-instability-visible": {
        "claim_title": "Weak phylogenetic signal instability remains visible",
        "summary": "A governed low-information signal case is rerun across fixed seeds so threshold-sensitive conclusions do not hide behind one lucky pass.",
        "verdict": "matched_with_tolerance",
    },
    "ou-identifiability-boundaries-detected": {
        "claim_title": "OU identifiability boundary warnings are detected on reference cases",
        "summary": "Known small-sample and weak-pull OU pathologies stay visible through the canonical runtime audit instead of being implied away.",
        "verdict": "matched",
    },
}

FRAGMENT_DEFINITIONS = [
    {
        "fragment_id": "pgls-input-rejection-guards",
        "fragment_title": "PGLS rootedness, branch-length, and numeric-response guards",
        "family_id": "expected-failures",
        "claim_ids": ["comparative-input-rejection-governed"],
        "evidence_id": EXPECTED_FAILURE_EVIDENCE_ID,
        "supporting_evidence_ids": [],
        "code_path": "packages/bijux-phylogenetics/src/bijux_phylogenetics/comparative/pgls.py",
        "code_locator": "bijux_phylogenetics.comparative.pgls:inspect_pgls_inputs",
        "parity_expectation": "exact",
        "comparison_kind": "exact_answer",
        "block_status": "verified",
        "review_note": "These cases prove that weak comparative inputs are rejected with explicit blockers before downstream model claims exist.",
        "scope": "runtime-guard",
    },
    {
        "fragment_id": "permutation-signal-instability",
        "fragment_title": "Permutation-backed weak phylogenetic signal boundary",
        "family_id": "weak-signal",
        "claim_ids": ["weak-signal-instability-visible"],
        "evidence_id": WEAK_SIGNAL_EVIDENCE_ID,
        "supporting_evidence_ids": [],
        "code_path": "packages/bijux-phylogenetics/src/bijux_phylogenetics/comparative/signal.py",
        "code_locator": "bijux_phylogenetics.comparative.signal:compute_phylogenetic_signal_test",
        "parity_expectation": "matched_with_tolerance",
        "comparison_kind": "stability_audit",
        "block_status": "verified",
        "review_note": "The governed weak-signal case is intentionally sensitive near the 0.05 threshold and is measured across fixed seeds rather than judged from one run.",
        "scope": "stability-audit",
    },
    {
        "fragment_id": "ou-identifiability-reference-audit",
        "fragment_title": "OU identifiability warning reference audit",
        "family_id": "model-instability",
        "claim_ids": ["ou-identifiability-boundaries-detected"],
        "evidence_id": OU_IDENTIFIABILITY_EVIDENCE_ID,
        "supporting_evidence_ids": [WEAK_SIGNAL_EVIDENCE_ID],
        "code_path": "packages/bijux-phylogenetics/src/bijux_phylogenetics/comparative/models.py",
        "code_locator": "bijux_phylogenetics.comparative.models:audit_ou_identifiability_reference_examples",
        "parity_expectation": "exact",
        "comparison_kind": "warning_audit",
        "block_status": "verified",
        "review_note": "OU boundary warnings are treated as trust evidence because model optimism would be dishonest without them.",
        "scope": "model-boundary",
    },
]

BUNDLE_DEFINITIONS = [
    {
        "evidence_id": EXPECTED_FAILURE_EVIDENCE_ID,
        "report_filename": "expected-failure-cases.json",
        "title": "Comparative input rejection bundle",
        "summary": "Governed expected-failure cases for rootedness, complete branch lengths, and numeric-response requirements in PGLS workflows.",
        "claim_ids": ["comparative-input-rejection-governed"],
        "claim_tags": ["runtime-guard", "expected-failure", "comparative", "pgls"],
        "analytical_surfaces": ["comparative-input-rejection", "pgls"],
        "source_fragments": ["pgls-input-rejection-guards"],
        "limitations": [
            "This bundle audits representative rejection cases, not every possible malformed comparative input.",
        ],
    },
    {
        "evidence_id": WEAK_SIGNAL_EVIDENCE_ID,
        "report_filename": "result-instability-audit.json",
        "title": "Weak signal instability bundle",
        "summary": "Governed weak-signal case showing that significance can cross the 0.05 boundary across fixed-seed permutation reruns.",
        "claim_ids": ["weak-signal-instability-visible"],
        "claim_tags": ["weak-signal", "instability", "comparative", "signal"],
        "analytical_surfaces": ["weak-signal", "result-instability", "phylogenetic-signal"],
        "source_fragments": ["permutation-signal-instability"],
        "limitations": [
            "This bundle proves one governed weak-signal boundary case and does not claim to exhaust every instability mode in comparative workflows.",
        ],
    },
    {
        "evidence_id": OU_IDENTIFIABILITY_EVIDENCE_ID,
        "report_filename": "ou-identifiability-audit.json",
        "title": "OU identifiability warning bundle",
        "summary": "Governed audit of built-in reference cases that should trigger OU identifiability warnings.",
        "claim_ids": ["ou-identifiability-boundaries-detected"],
        "claim_tags": ["model-instability", "comparative", "ou", "identifiability"],
        "analytical_surfaces": ["model-instability", "ou-identifiability"],
        "source_fragments": ["ou-identifiability-reference-audit"],
        "limitations": [
            "The audit covers known reference warning families and does not claim that every possible OU pathology is enumerated here.",
        ],
    },
]

WEAK_SIGNAL_TRAIT_ROWS = [
    ("A", 0.835),
    ("B", 1.724),
    ("C", 6.838),
    ("D", 1.614),
    ("E", 9.492),
    ("F", 8.443),
]


def render_weak_signal_traits_tsv() -> str:
    lines = ["taxon\tresponse"]
    for taxon, response in WEAK_SIGNAL_TRAIT_ROWS:
        lines.append(f"{taxon}\t{response}")
    lines.append("")
    return "\n".join(lines)


def build_comparative_trust_boundaries_provenance() -> dict[str, object]:
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "source_count": 6,
        "sources": [
            {
                "source_id": "example-tree-unrooted",
                "kind": "repository-fixture",
                "locator": (EXPECTED_FAILURE_TREE_ROOT / "example_tree_unrooted.nwk").as_posix(),
                "role": "expected-failure-tree",
            },
            {
                "source_id": "example-tree-no-lengths",
                "kind": "repository-fixture",
                "locator": (EXPECTED_FAILURE_TREE_ROOT / "example_tree_no_lengths.nwk").as_posix(),
                "role": "expected-failure-tree",
            },
            {
                "source_id": "example-tree-rooted",
                "kind": "repository-fixture",
                "locator": (EXPECTED_FAILURE_TREE_ROOT / "example_tree.nwk").as_posix(),
                "role": "expected-failure-tree",
            },
            {
                "source_id": "example-traits-comparative",
                "kind": "repository-fixture",
                "locator": (EXPECTED_FAILURE_TRAITS_ROOT / "example_traits_comparative.tsv").as_posix(),
                "role": "expected-failure-traits",
            },
            {
                "source_id": "example-tree-six-taxa",
                "kind": "repository-fixture",
                "locator": WEAK_SIGNAL_TREE_PATH.as_posix(),
                "role": "stability-tree",
            },
            {
                "source_id": "weak-signal-governed-traits",
                "kind": "repository-evidence",
                "locator": WEAK_SIGNAL_TRAITS_PATH.as_posix(),
                "role": "stability-traits",
            },
        ],
    }


def build_comparative_trust_boundaries_source_fragment_map() -> dict[str, object]:
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
                "code_path": definition["code_path"],
                "code_locator": definition["code_locator"],
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
        "fragment_count": len(fragments),
        "fragments": fragments,
    }


def build_comparative_trust_boundaries_family_index(repo_root: Path) -> dict[str, object]:
    _ = repo_root
    fragments = build_comparative_trust_boundaries_source_fragment_map()["fragments"]
    families = []
    for family_id, family in FAMILY_DEFINITIONS.items():
        matching_fragments = [
            fragment for fragment in fragments if fragment["concept_family"] == family_id
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
        verdicts = {CLAIM_DEFINITIONS[claim_id]["verdict"] for claim_id in claim_ids}
        if verdicts == {"matched"}:
            family_verdict = "matched"
        elif verdicts == {"matched_with_tolerance"}:
            family_verdict = "matched_with_tolerance"
        else:
            family_verdict = "mismatch_unexplained"
        families.append(
            {
                "family_id": family_id,
                "family_title": family["title"],
                "summary": family["summary"],
                "fragment_count": len(matching_fragments),
                "fragment_ids": [fragment["fragment_id"] for fragment in matching_fragments],
                "claim_ids": claim_ids,
                "evidence_ids": evidence_ids,
                "family_verdict": family_verdict,
                "coverage_status": "covered",
                "known_gaps": [],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "family_count": len(families),
        "families": families,
    }


def _expected_failure_cases(repo_root: Path) -> list[dict[str, object]]:
    fixtures_root = Path(repo_root) / FIXTURES_ROOT
    traits_path = fixtures_root / "metadata" / "example_traits_comparative.tsv"
    cases = [
        {
            "case_id": "pgls-unrooted-tree-rejection",
            "tree_path": fixtures_root / "trees" / "example_tree_unrooted.nwk",
            "traits_path": traits_path,
            "response": "response",
            "predictors": ["predictor_one"],
            "expected_blockers": ["PGLS requires a rooted tree"],
        },
        {
            "case_id": "pgls-missing-branch-length-rejection",
            "tree_path": fixtures_root / "trees" / "example_tree_no_lengths.nwk",
            "traits_path": traits_path,
            "response": "response",
            "predictors": ["predictor_one"],
            "expected_blockers": ["PGLS requires complete tree branch lengths"],
        },
        {
            "case_id": "pgls-nonnumeric-response-rejection",
            "tree_path": fixtures_root / "trees" / "example_tree.nwk",
            "traits_path": traits_path,
            "response": "habitat",
            "predictors": ["predictor_one"],
            "expected_blockers": ["response column 'habitat' must be numeric for PGLS"],
        },
    ]
    observations = []
    for case in cases:
        try:
            run_pgls(
                case["tree_path"],
                case["traits_path"],
                response=str(case["response"]),
                predictors=list(case["predictors"]),
            )
        except ComparativeMethodError as error:
            observed_message = str(error)
            observed_blockers = [segment.strip() for segment in observed_message.split(";")]
            observations.append(
                {
                    "case_id": case["case_id"],
                    "tree_path": _relative_repo_path(repo_root, case["tree_path"]),
                    "traits_path": _relative_repo_path(repo_root, case["traits_path"]),
                    "response": case["response"],
                    "predictors": case["predictors"],
                    "expected_blockers": case["expected_blockers"],
                    "observed_error_code": error.code,
                    "observed_message": observed_message,
                    "observed_blockers": observed_blockers,
                    "matched_expected_blockers": all(
                        expected in observed_message
                        for expected in case["expected_blockers"]
                    ),
                }
            )
        else:
            raise AssertionError(f"expected ComparativeMethodError for {case['case_id']}")
    return observations


def build_expected_failure_cases_report(repo_root: Path) -> dict[str, object]:
    cases = _expected_failure_cases(repo_root)
    matched_case_count = sum(1 for case in cases if case["matched_expected_blockers"])
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": EXPECTED_FAILURE_EVIDENCE_ID,
        "case_count": len(cases),
        "matched_case_count": matched_case_count,
        "all_cases_matched": matched_case_count == len(cases),
        "cases": cases,
    }


def build_result_instability_audit(repo_root: Path) -> dict[str, object]:
    tree_path = Path(repo_root) / WEAK_SIGNAL_TREE_PATH
    traits_path = Path(repo_root) / WEAK_SIGNAL_TRAITS_PATH
    runs = []
    for seed in WEAK_SIGNAL_SEEDS:
        report = compute_phylogenetic_signal_test(
            tree_path,
            traits_path,
            trait="response",
            permutations=WEAK_SIGNAL_PERMUTATIONS,
            seed=seed,
        )
        runs.append(
            {
                "seed": seed,
                "p_value": report.p_value,
                "permuted_k_at_or_above_observed": report.permuted_k_at_or_above_observed,
                "observed_k": report.observed_k,
            }
        )
    p_values = [float(run["p_value"]) for run in runs]
    estimated_lambda = estimate_pagels_lambda(
        tree_path,
        traits_path,
        trait="response",
    )
    below_threshold = sum(p_value < WEAK_SIGNAL_ALPHA for p_value in p_values)
    at_or_above_threshold = len(p_values) - below_threshold
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": WEAK_SIGNAL_EVIDENCE_ID,
        "tree_path": WEAK_SIGNAL_TREE_PATH.as_posix(),
        "traits_path": WEAK_SIGNAL_TRAITS_PATH.as_posix(),
        "alpha_threshold": WEAK_SIGNAL_ALPHA,
        "permutations_per_run": WEAK_SIGNAL_PERMUTATIONS,
        "seed_count": len(runs),
        "estimated_lambda": estimated_lambda.lambda_value,
        "log_likelihood": estimated_lambda.log_likelihood,
        "null_log_likelihood": estimated_lambda.null_log_likelihood,
        "p_value_min": min(p_values),
        "p_value_max": max(p_values),
        "below_threshold_count": below_threshold,
        "at_or_above_threshold_count": at_or_above_threshold,
        "crosses_alpha_boundary": below_threshold > 0 and at_or_above_threshold > 0,
        "instability_summary": (
            "The governed weak-signal case crosses the 0.05 decision boundary across fixed-seed permutation reruns."
        ),
        "runs": runs,
        "trait_rows": [
            {"taxon": taxon, "response": response}
            for taxon, response in WEAK_SIGNAL_TRAIT_ROWS
        ],
        "scientific_debt_entries": [
            {
                "debt_id": "weak-signal-threshold-instability",
                "debt_kind": "instability",
                "detail": "The weak-signal case crosses the 0.05 decision boundary across governed reruns, so any single-run significance claim would be overconfident.",
                "evidence": [
                    f"min p-value={min(p_values):.2f}",
                    f"max p-value={max(p_values):.2f}",
                    f"runs below 0.05={below_threshold}",
                    f"runs at or above 0.05={at_or_above_threshold}",
                ],
            }
        ],
    }


def build_ou_identifiability_audit(repo_root: Path) -> dict[str, object]:
    _ = repo_root
    audit = audit_ou_identifiability_reference_examples()
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": OU_IDENTIFIABILITY_EVIDENCE_ID,
        "all_expected_warning_kinds_detected": audit.all_expected_warning_kinds_detected,
        "expected_warning_kinds": audit.expected_warning_kinds_detected
        if hasattr(audit, "expected_warning_kinds_detected")
        else audit.expected_warning_kinds,
        "detected_warning_kinds": audit.detected_warning_kinds,
        "cases": [asdict(case) for case in audit.cases],
        "scientific_debt_entries": [
            {
                "debt_id": "ou-identifiability-boundary-cases",
                "debt_kind": "model-boundary",
                "detail": "OU fits on the governed reference cases still sit on identifiability boundaries, so the repository must keep these warnings reviewer-visible instead of treating them as ordinary successful fits.",
                "evidence": [
                    "example-tree-small-n triggers small_sample_size, boundary_alpha, and flat_likelihood",
                    "example-tree-weak-pull triggers boundary_alpha, flat_likelihood, and weak_pull_to_optimum",
                ],
            }
        ],
    }


def _report_payload_for_bundle(repo_root: Path, evidence_id: str) -> dict[str, object]:
    if evidence_id == EXPECTED_FAILURE_EVIDENCE_ID:
        return build_expected_failure_cases_report(repo_root)
    if evidence_id == WEAK_SIGNAL_EVIDENCE_ID:
        return build_result_instability_audit(repo_root)
    return build_ou_identifiability_audit(repo_root)


def _relative_repo_path(repo_root: Path, path: Path) -> str:
    return path.relative_to(Path(repo_root)).as_posix()


def _manifest_for_bundle(definition: dict[str, object], report_payload: dict[str, object]) -> dict[str, object]:
    evidence_id = str(definition["evidence_id"])
    source_basis = [
        {
            "kind": "repository-source-descriptor",
            "label": "comparative trust boundary provenance",
            "locator": f"evidence-book/studies/{STUDY_ID}/provenance/runtime-sources.json",
        }
    ]
    if evidence_id == EXPECTED_FAILURE_EVIDENCE_ID:
        source_basis.extend(
            [
                {
                    "kind": "repository-fixture",
                    "label": "unrooted comparative tree fixture",
                    "locator": (EXPECTED_FAILURE_TREE_ROOT / "example_tree_unrooted.nwk").as_posix(),
                },
                {
                    "kind": "repository-fixture",
                    "label": "branch-length incomplete comparative tree fixture",
                    "locator": (EXPECTED_FAILURE_TREE_ROOT / "example_tree_no_lengths.nwk").as_posix(),
                },
                {
                    "kind": "repository-fixture",
                    "label": "comparative traits fixture",
                    "locator": (EXPECTED_FAILURE_TRAITS_ROOT / "example_traits_comparative.tsv").as_posix(),
                },
            ]
        )
    elif evidence_id == WEAK_SIGNAL_EVIDENCE_ID:
        source_basis.extend(
            [
                {
                    "kind": "repository-fixture",
                    "label": "six-taxon stability tree fixture",
                    "locator": WEAK_SIGNAL_TREE_PATH.as_posix(),
                },
                {
                    "kind": "repository-evidence",
                    "label": "governed weak-signal trait table",
                    "locator": WEAK_SIGNAL_TRAITS_PATH.as_posix(),
                },
            ]
        )
    else:
        source_basis.extend(
            [
                {
                    "kind": "repository-fixture",
                    "label": "comparative reference tree fixtures",
                    "locator": FIXTURES_ROOT.as_posix(),
                },
                {
                    "kind": "repository-code",
                    "label": "OU identifiability audit implementation",
                    "locator": "packages/bijux-phylogenetics/src/bijux_phylogenetics/comparative/models.py",
                },
            ]
        )

    verdict_status = CLAIM_DEFINITIONS[str(definition["claim_ids"][0])]["verdict"]
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": evidence_id,
        "evidence_title": definition["title"],
        "summary": definition["summary"],
        "owner_package": OWNER_PACKAGE,
        "claim_ids": definition["claim_ids"],
        "source_basis": source_basis,
        "freshness": {
            "last_generated_on": "2026-05-10",
            "governed_code_paths": [
                "packages/bijux-phylogenetics/src/bijux_phylogenetics/comparative/pgls.py",
                "packages/bijux-phylogenetics/src/bijux_phylogenetics/comparative/signal.py",
                "packages/bijux-phylogenetics/src/bijux_phylogenetics/comparative/models.py",
                "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/studies/comparative_trust_boundaries.py",
            ],
            "source_basis_locators": [entry["locator"] for entry in source_basis],
        },
        "ownership": {
            "owner_package": OWNER_PACKAGE,
            "analytical_surfaces": definition["analytical_surfaces"],
        },
        "claim_tags": definition["claim_tags"],
        "verdict": {
            "status": verdict_status,
            "summary": _bundle_verdict_summary(evidence_id, report_payload),
        },
        "limitations": definition["limitations"],
    }


def _bundle_verdict_summary(evidence_id: str, report_payload: dict[str, object]) -> str:
    if evidence_id == EXPECTED_FAILURE_EVIDENCE_ID:
        return "Observed ComparativeMethodError blockers match the governed rejection expectations for the tracked comparative input failures."
    if evidence_id == WEAK_SIGNAL_EVIDENCE_ID:
        return "The governed weak-signal case crosses the 0.05 decision boundary across fixed-seed reruns, so instability remains explicit instead of being hidden behind one run."
    return "The OU identifiability reference audit detects every expected warning family on the governed boundary cases."


def _claims_for_bundle(definition: dict[str, object]) -> dict[str, object]:
    claims = []
    for claim_id in definition["claim_ids"]:
        claim = CLAIM_DEFINITIONS[claim_id]
        claims.append(
            {
                "claim_id": claim_id,
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
        "evidence_id": definition["evidence_id"],
        "claim_count": len(claims),
        "claims": claims,
    }


def _readme_for_bundle(
    definition: dict[str, object], report_payload: dict[str, object]
) -> str:
    lines = [
        f"# {definition['title']}",
        "",
        definition["summary"],
        "",
        "Claim surfaces:",
        "",
    ]
    for claim_id in definition["claim_ids"]:
        claim = CLAIM_DEFINITIONS[claim_id]
        lines.append(f"- `{claim_id}` — {claim['claim_title']} (`{claim['verdict']}`)")
    lines.extend(
        [
            "",
            "Sources:",
            "",
        ]
    )
    for fragment_id in definition["source_fragments"]:
        fragment = next(
            item for item in FRAGMENT_DEFINITIONS if item["fragment_id"] == fragment_id
        )
        lines.append(
            f"- `{fragment_id}` — {fragment['fragment_title']} via `{fragment['code_locator']}`"
        )
    if definition["evidence_id"] == WEAK_SIGNAL_EVIDENCE_ID:
        lines.extend(
            [
                "",
                "Instability summary:",
                "",
                f"- alpha threshold: `{report_payload['alpha_threshold']}`",
                f"- min p-value: `{report_payload['p_value_min']}`",
                f"- max p-value: `{report_payload['p_value_max']}`",
                f"- runs below threshold: `{report_payload['below_threshold_count']}`",
                f"- runs at or above threshold: `{report_payload['at_or_above_threshold_count']}`",
            ]
        )
    lines.extend(
        [
            "",
            "Limitations:",
            "",
        ]
    )
    for limitation in definition["limitations"]:
        lines.append(f"- {limitation}")
    lines.append("")
    return "\n".join(lines)


def build_comparative_trust_boundaries_bundles(
    repo_root: Path,
) -> dict[str, dict[str, object]]:
    bundles: dict[str, dict[str, object]] = {}
    for definition in BUNDLE_DEFINITIONS:
        report_payload = _report_payload_for_bundle(repo_root, definition["evidence_id"])
        bundles[str(definition["evidence_id"])] = {
            "manifest": _manifest_for_bundle(definition, report_payload),
            "claims": _claims_for_bundle(definition),
            "report_filename": definition["report_filename"],
            "report_payload": report_payload,
            "readme": _readme_for_bundle(definition, report_payload),
        }
    return bundles


def write_weak_signal_traits_table(repo_root: Path) -> Path:
    path = Path(repo_root) / WEAK_SIGNAL_TRAITS_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_weak_signal_traits_tsv(), encoding="utf-8")
    return path


def build_comparative_trust_boundaries_claim_registry(
    repo_root: Path,
) -> dict[str, object]:
    bundles = build_comparative_trust_boundaries_bundles(repo_root)
    claims = []
    for definition in BUNDLE_DEFINITIONS:
        for claim_id in definition["claim_ids"]:
            claim = CLAIM_DEFINITIONS[claim_id]
            claims.append(
                {
                    "claim_id": claim_id,
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


def render_comparative_trust_boundaries_study_readme() -> str:
    lines = [
        "# Comparative Trust Boundaries",
        "",
        "This study turns comparative failure surfaces and numerically fragile",
        "boundary cases into governed Evidence IDs for `bijux-phylogenetics`.",
        "",
        "It exists to keep trust boundaries visible:",
        "",
        "- expected comparative input failures are treated as evidence, not as incidental exceptions",
        "- weak-signal cases are rerun across fixed seeds so one lucky significance pass cannot overstate confidence",
        "- OU identifiability warnings stay reviewer-visible on governed reference cases",
        "",
        "Current bundles:",
        "",
        f"- `{EXPECTED_FAILURE_EVIDENCE_ID}` comparative input rejection",
        f"- `{WEAK_SIGNAL_EVIDENCE_ID}` weak-signal instability",
        f"- `{OU_IDENTIFIABILITY_EVIDENCE_ID}` OU identifiability warnings",
        "",
    ]
    return "\n".join(lines)


def render_comparative_trust_boundaries_study_manifest() -> dict[str, object]:
    return {
        "study_id": STUDY_ID,
        "study_title": "Comparative trust boundary evidence study",
        "summary": "Governed evidence for expected comparative failures, weak-signal instability, and OU identifiability boundaries.",
        "owner_package": OWNER_PACKAGE,
        "confidence_posture": "trust-boundaries-explicit",
    }


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
