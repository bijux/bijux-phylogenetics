from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.ancestral.discrete.review.ordering import (
    summarize_ordered_discrete_reconstruction,
    summarize_ordered_discrete_report,
)
from bijux_phylogenetics.ancestral.discrete.review.transition_constraints import (
    summarize_irreversible_discrete_reconstruction,
    summarize_irreversible_discrete_report,
)
from bijux_phylogenetics.ancestral.sensitivity import (
    summarize_ancestral_root_sensitivity,
    summarize_ancestral_root_sensitivity_report,
)
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_trait_table_fixture,
    get_shared_tree_fixture,
)


@dataclass(slots=True)
class DiscreteAncestralReferenceProbabilityRow:
    """One state-probability comparison for a governed discrete ancestral case."""

    node: str
    state: str
    expected_probability: float
    observed_probability: float
    absolute_delta: float


@dataclass(slots=True)
class DiscreteAncestralReferenceObservation:
    """One governed reference or policy observation for discrete ancestral review."""

    case_id: str
    category: str
    source: str
    model: str
    tolerance: float
    expected_metrics: dict[str, object]
    observed_metrics: dict[str, object]
    probability_rows: list[DiscreteAncestralReferenceProbabilityRow]
    passed: bool
    notes: list[str]


@dataclass(slots=True)
class DiscreteAncestralReferenceValidationReport:
    """Integrated governed validation over discrete ancestral reference surfaces."""

    case_count: int
    external_case_count: int
    all_passed: bool
    observations: list[DiscreteAncestralReferenceObservation]


def validate_discrete_ancestral_reference_examples(
    *,
    tolerance: float | None = None,
) -> DiscreteAncestralReferenceValidationReport:
    """Validate discrete ancestral reconstruction against governed reference examples."""
    observations = _validate_external_probability_cases(
        override_tolerance=tolerance,
    )
    observations.append(_validate_root_prior_behavior())
    observations.append(_validate_ambiguous_state_behavior())
    observations.append(_validate_ordered_constraint_behavior())
    observations.append(_validate_irreversible_constraint_behavior())
    external_case_count = sum(
        observation.category == "external-marginal-probability"
        for observation in observations
    )
    return DiscreteAncestralReferenceValidationReport(
        case_count=len(observations),
        external_case_count=external_case_count,
        all_passed=all(observation.passed for observation in observations),
        observations=observations,
    )


def _validate_external_probability_cases(
    *,
    override_tolerance: float | None,
) -> list[DiscreteAncestralReferenceObservation]:
    fixture = _load_ace_reference_fixture()
    repository_root = _repository_root()
    observations: list[DiscreteAncestralReferenceObservation] = []
    for case in fixture["cases"]:
        report = reconstruct_discrete_ancestral_states(
            repository_root / case["tree"],
            repository_root / case["table"],
            trait=case["trait"],
            taxon_column=case["taxon_column"],
            model=case["model"],
        )
        tolerance = (
            float(case["abs_tolerance"])
            if override_tolerance is None
            else override_tolerance
        )
        observed_rows = {
            estimate.node: estimate
            for estimate in report.estimates
            if not estimate.is_tip
        }
        probability_rows: list[DiscreteAncestralReferenceProbabilityRow] = []
        for row in case["rows"]:
            observed_probabilities = observed_rows[row["node"]].state_probabilities
            for state, expected_probability in row["probabilities"].items():
                observed_probability = observed_probabilities[state]
                probability_rows.append(
                    DiscreteAncestralReferenceProbabilityRow(
                        node=str(row["node"]),
                        state=str(state),
                        expected_probability=float(expected_probability),
                        observed_probability=observed_probability,
                        absolute_delta=abs(observed_probability - expected_probability),
                    )
                )
        root_observed = _root_estimate(report)
        root_expected_probabilities = {
            str(state): float(probability)
            for state, probability in case["rows"][0]["probabilities"].items()
        }
        max_probability_delta = max(
            (row.absolute_delta for row in probability_rows),
            default=0.0,
        )
        passed = max_probability_delta <= tolerance
        observations.append(
            DiscreteAncestralReferenceObservation(
                case_id=str(case["case_id"]),
                category="external-marginal-probability",
                source=str(case.get("source", "ape::ace")),
                model=str(case["model"]),
                tolerance=tolerance,
                expected_metrics={
                    "node_count": len(case["rows"]),
                    "root_state": _expected_most_likely_state(
                        root_expected_probabilities
                    ),
                    "root_ambiguous": _expected_ambiguous(root_expected_probabilities),
                },
                observed_metrics={
                    "node_count": len(observed_rows),
                    "root_state": root_observed.most_likely_state,
                    "root_ambiguous": _expected_ambiguous(
                        root_observed.state_probabilities
                    ),
                    "max_probability_delta": max_probability_delta,
                },
                probability_rows=probability_rows,
                passed=passed,
                notes=[]
                if passed
                else ["one or more node-state probabilities drifted beyond tolerance"],
            )
        )
    return observations


def _validate_root_prior_behavior() -> DiscreteAncestralReferenceObservation:
    tree_fixture = get_shared_tree_fixture("balanced_rooted_ultrametric")
    trait_fixture = get_shared_trait_table_fixture("multistate_discrete_match")
    report = summarize_ancestral_root_sensitivity(
        tree_fixture.path,
        trait_fixture.path,
        trait="region",
        model="equal-rates",
        fixed_root_state="island",
    )
    summary = summarize_ancestral_root_sensitivity_report(report)
    assumption_rows = {row.assumption_id: row for row in report.assumption_rows}
    expected_metrics = {
        "assumption_count": 3,
        "state_changed_node_count": 2,
        "support_changed_node_count": 1,
        "top_sensitive_node": "A|B|C|D",
        "equal_root_state": "north",
        "empirical_root_state": "north",
        "fixed_root_state": "island",
        "fixed_root_confidence": 1.0,
    }
    observed_metrics = {
        "assumption_count": summary.assumption_count,
        "state_changed_node_count": summary.state_changed_node_count,
        "support_changed_node_count": summary.support_changed_node_count,
        "top_sensitive_node": summary.top_sensitive_node,
        "equal_root_state": assumption_rows["equal_root_prior"].root_most_likely_state,
        "empirical_root_state": assumption_rows[
            "empirical_root_prior"
        ].root_most_likely_state,
        "fixed_root_state": assumption_rows["fixed_root_state"].root_most_likely_state,
        "fixed_root_confidence": assumption_rows["fixed_root_state"].root_confidence,
    }
    passed = (
        observed_metrics["assumption_count"] == expected_metrics["assumption_count"]
        and observed_metrics["state_changed_node_count"]
        == expected_metrics["state_changed_node_count"]
        and observed_metrics["support_changed_node_count"]
        == expected_metrics["support_changed_node_count"]
        and observed_metrics["top_sensitive_node"]
        == expected_metrics["top_sensitive_node"]
        and observed_metrics["equal_root_state"] == expected_metrics["equal_root_state"]
        and observed_metrics["empirical_root_state"]
        == expected_metrics["empirical_root_state"]
        and observed_metrics["fixed_root_state"] == expected_metrics["fixed_root_state"]
        and math.isclose(
            float(observed_metrics["fixed_root_confidence"]),
            float(expected_metrics["fixed_root_confidence"]),
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    )
    return DiscreteAncestralReferenceObservation(
        case_id="root_prior_sensitivity_policy",
        category="root-prior-policy",
        source="owned root-prior policy",
        model="equal-rates",
        tolerance=1e-12,
        expected_metrics=expected_metrics,
        observed_metrics=observed_metrics,
        probability_rows=[],
        passed=passed,
        notes=[]
        if passed
        else [
            "root-prior sensitivity behavior diverged from the governed policy surface"
        ],
    )


def _validate_ambiguous_state_behavior() -> DiscreteAncestralReferenceObservation:
    tree_fixture = get_shared_tree_fixture("balanced_rooted_ultrametric")
    trait_fixture = get_shared_trait_table_fixture("categorical_predictor_match")
    report = reconstruct_discrete_ancestral_states(
        tree_fixture.path,
        trait_fixture.path,
        trait="habitat",
        model="fitch",
    )
    internal_rows = {
        estimate.node: estimate for estimate in report.estimates if not estimate.is_tip
    }
    root = internal_rows["A|B|C|D"]
    expected_metrics = {
        "ambiguous_nodes": ["A|B|C|D"],
        "root_state_set": ["forest", "tundra"],
        "root_probabilities": {"forest": 0.5, "tundra": 0.5},
        "weak_support_nodes": ["A|B|C|D"],
    }
    observed_metrics = {
        "ambiguous_nodes": sorted(
            estimate.node for estimate in internal_rows.values() if estimate.ambiguous
        ),
        "root_state_set": list(root.state_set),
        "root_probabilities": dict(root.state_probabilities),
        "weak_support_nodes": list(report.weak_support_nodes),
    }
    probability_rows = [
        DiscreteAncestralReferenceProbabilityRow(
            node="A|B|C|D",
            state=state,
            expected_probability=expected_probability,
            observed_probability=root.state_probabilities[state],
            absolute_delta=abs(root.state_probabilities[state] - expected_probability),
        )
        for state, expected_probability in expected_metrics[
            "root_probabilities"
        ].items()
    ]
    max_probability_delta = max(row.absolute_delta for row in probability_rows)
    passed = (
        observed_metrics["ambiguous_nodes"] == expected_metrics["ambiguous_nodes"]
        and observed_metrics["root_state_set"] == expected_metrics["root_state_set"]
        and observed_metrics["weak_support_nodes"]
        == expected_metrics["weak_support_nodes"]
        and max_probability_delta <= 1e-12
    )
    return DiscreteAncestralReferenceObservation(
        case_id="ambiguous_internal_states",
        category="ambiguity-policy",
        source="owned ambiguity policy",
        model="fitch",
        tolerance=1e-12,
        expected_metrics=expected_metrics,
        observed_metrics=observed_metrics,
        probability_rows=probability_rows,
        passed=passed,
        notes=[]
        if passed
        else [
            "ambiguous internal states were not reported with the governed state sets or probabilities"
        ],
    )


def _validate_ordered_constraint_behavior() -> DiscreteAncestralReferenceObservation:
    tree_fixture = get_shared_tree_fixture("balanced_rooted_ultrametric")
    trait_fixture = get_shared_trait_table_fixture("multistate_discrete_match")
    report = summarize_ordered_discrete_reconstruction(
        tree_fixture.path,
        trait_fixture.path,
        trait="region",
        model="equal-rates",
        ordered_states=["north", "south", "island"],
    )
    summary = summarize_ordered_discrete_report(report)
    north_to_island = next(
        row
        for row in report.transition_rows
        if row.source_state == "north" and row.target_state == "island"
    )
    expected_metrics = {
        "restricted_transition_count": 2,
        "preferred_ordering": "ordered",
        "north_to_island_allowed": False,
        "north_to_island_rate": 0.0,
    }
    observed_metrics = {
        "restricted_transition_count": summary.restricted_transition_count,
        "preferred_ordering": summary.preferred_ordering,
        "north_to_island_allowed": north_to_island.ordered_transition_allowed,
        "north_to_island_rate": north_to_island.ordered_rate,
    }
    passed = (
        observed_metrics["restricted_transition_count"]
        == expected_metrics["restricted_transition_count"]
        and observed_metrics["preferred_ordering"]
        == expected_metrics["preferred_ordering"]
        and observed_metrics["north_to_island_allowed"]
        == expected_metrics["north_to_island_allowed"]
        and math.isclose(
            float(observed_metrics["north_to_island_rate"]),
            float(expected_metrics["north_to_island_rate"]),
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    )
    return DiscreteAncestralReferenceObservation(
        case_id="ordered_transition_constraints",
        category="ordered-constraint-policy",
        source="owned ordered-state policy",
        model="equal-rates",
        tolerance=1e-12,
        expected_metrics=expected_metrics,
        observed_metrics=observed_metrics,
        probability_rows=[],
        passed=passed,
        notes=[]
        if passed
        else [
            "ordered-state constraints no longer produce the governed transition restriction surface"
        ],
    )


def _validate_irreversible_constraint_behavior() -> (
    DiscreteAncestralReferenceObservation
):
    fixtures_root = _fixtures_root()
    report = summarize_irreversible_discrete_reconstruction(
        fixtures_root / "trees/example_tree_six_taxa.nwk",
        fixtures_root / "metadata/example_traits_irreversible_loss.tsv",
        trait="state",
        model="all-rates-different",
        allowed_transition_pairs=[("present", "absent")],
    )
    summary = summarize_irreversible_discrete_report(report)
    absent_to_present = next(
        row
        for row in report.transition_rows
        if row.source_state == "absent" and row.target_state == "present"
    )
    expected_metrics = {
        "preferred_constraint": "constrained",
        "forbidden_transition_count": 1,
        "absent_to_present_allowed": False,
        "absent_to_present_rate": 0.0,
    }
    observed_metrics = {
        "preferred_constraint": summary.preferred_constraint,
        "forbidden_transition_count": summary.forbidden_transition_count,
        "absent_to_present_allowed": absent_to_present.constrained_transition_allowed,
        "absent_to_present_rate": absent_to_present.constrained_rate,
    }
    passed = (
        observed_metrics["preferred_constraint"]
        == expected_metrics["preferred_constraint"]
        and observed_metrics["forbidden_transition_count"]
        == expected_metrics["forbidden_transition_count"]
        and observed_metrics["absent_to_present_allowed"]
        == expected_metrics["absent_to_present_allowed"]
        and math.isclose(
            float(observed_metrics["absent_to_present_rate"]),
            float(expected_metrics["absent_to_present_rate"]),
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    )
    return DiscreteAncestralReferenceObservation(
        case_id="irreversible_loss_constraints",
        category="irreversible-constraint-policy",
        source="owned irreversible-transition policy",
        model="all-rates-different",
        tolerance=1e-12,
        expected_metrics=expected_metrics,
        observed_metrics=observed_metrics,
        probability_rows=[],
        passed=passed,
        notes=[]
        if passed
        else [
            "irreversible-transition constraints no longer preserve the governed one-way loss behavior"
        ],
    )


def _load_ace_reference_fixture() -> dict[str, object]:
    return json.loads(
        (_fixtures_root() / "expected/discrete_ancestral_ace_reference.json").read_text(
            encoding="utf-8"
        )
    )


def _fixtures_root() -> Path:
    return (
        _repository_root() / "packages" / "bijux-phylogenetics" / "tests" / "fixtures"
    )


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[7]


def _root_estimate(report):
    internal_rows = [estimate for estimate in report.estimates if not estimate.is_tip]
    if not internal_rows:
        raise ValueError("discrete ancestral report did not contain internal nodes")
    return max(internal_rows, key=lambda estimate: len(estimate.descendant_taxa))


def _expected_most_likely_state(state_probabilities: dict[str, float]) -> str:
    best_probability = max(state_probabilities.values())
    tied_states = [
        state
        for state, probability in state_probabilities.items()
        if math.isclose(probability, best_probability, rel_tol=1e-9, abs_tol=1e-12)
    ]
    return sorted(tied_states)[0]


def _expected_ambiguous(state_probabilities: dict[str, float]) -> bool:
    best_probability = max(state_probabilities.values())
    return (
        sum(
            math.isclose(probability, best_probability, rel_tol=1e-9, abs_tol=1e-12)
            for probability in state_probabilities.values()
        )
        > 1
    )
