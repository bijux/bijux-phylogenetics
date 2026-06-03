from __future__ import annotations

from pathlib import Path

from .definitions import EVIDENCE_ID, FAMILY_DEFINITIONS, STUDY_ID
from .registry import build_primate_source_fragment_map


def build_primate_parity_policy(repo_root: Path) -> dict[str, object]:
    fragment_map = build_primate_source_fragment_map(repo_root)
    policies: list[dict[str, object]] = []
    for family_id in [
        "workflow-contracts",
        "data-preparation",
        "tree-operations",
        "visual-surfaces",
        "simulation-inputs",
        "comparative-signal",
        "ancestral-reconstruction",
        "artifact-provenance",
    ]:
        family_fragments = [
            fragment
            for fragment in fragment_map["fragments"]
            if fragment["concept_family"] == family_id
        ]
        if family_id in {
            "workflow-contracts",
            "visual-surfaces",
            "artifact-provenance",
        }:
            expectation = "not_comparable"
            metric_tolerances = []
            rule = "Tracked for reviewer transparency without a numerical parity claim."
        elif family_id in {"data-preparation", "tree-operations", "simulation-inputs"}:
            expectation = "exact"
            metric_tolerances = [{"metric_kind": "default", "tolerance_abs_diff": 0.0}]
            rule = (
                "Scalar counts, boolean checks, and frozen inputs must match exactly."
            )
        elif family_id == "comparative-signal":
            expectation = "statistical_tolerance"
            metric_tolerances = [
                {"metric_kind": "lambda_value", "tolerance_abs_diff": 0.001},
                {"metric_kind": "log_likelihood", "tolerance_abs_diff": 0.001},
                {"metric_kind": "likelihood_ratio", "tolerance_abs_diff": 0.0001},
                {"metric_kind": "p_value", "tolerance_abs_diff": 1e-12},
                {"metric_kind": "covariance_entry", "tolerance_abs_diff": 1e-8},
            ]
            rule = "Likelihood-based metrics must stay within governed floating-point tolerance while preserving the same scientific conclusion."
        else:
            expectation = "near_exact"
            metric_tolerances = [
                {"metric_kind": "point_estimate", "tolerance_abs_diff": 1e-9},
                {"metric_kind": "interval_bound", "tolerance_abs_diff": 1e-9},
                {"metric_kind": "count", "tolerance_abs_diff": 0.0},
            ]
            rule = "Internal-node estimates may differ only by floating-point noise; derived counts must still match exactly."

        policies.append(
            {
                "family_id": family_id,
                "family_title": FAMILY_DEFINITIONS[family_id]["title"],
                "parity_expectation": expectation,
                "metric_tolerances": metric_tolerances,
                "source_fragments": [
                    fragment["fragment_id"] for fragment in family_fragments
                ],
                "rule": rule,
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": EVIDENCE_ID,
        "policy_count": len(policies),
        "policies": policies,
    }
