from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.continuous import fit_ornstein_uhlenbeck_model


@dataclass(slots=True)
class ComparativeOUIdentifiabilityCase:
    """One fixture-level OU identifiability audit case."""

    case: str
    alpha: float
    warning_kinds: list[str]


@dataclass(slots=True)
class ComparativeOUIdentifiabilityAudit:
    """Whether the expected OU failure modes are detected on reference fixtures."""

    cases: list[ComparativeOUIdentifiabilityCase]
    detected_warning_kinds: list[str]
    expected_warning_kinds: list[str]
    all_expected_warning_kinds_detected: bool


def audit_ou_identifiability_reference_examples() -> ComparativeOUIdentifiabilityAudit:
    """Verify that all expected OU warning modes are triggered on built-in fixtures."""
    root = Path(__file__).resolve().parents[4] / "tests/fixtures"
    cases = [
        (
            "example-tree-small-n",
            root / "trees/example_tree.nwk",
            root / "metadata/example_traits_comparative.tsv",
            "response",
        ),
        (
            "example-tree-weak-pull",
            root / "trees/example_tree_six_taxa.nwk",
            root / "metadata/example_traits_comparative_multiple.tsv",
            "response_growth",
        ),
    ]
    observations: list[ComparativeOUIdentifiabilityCase] = []
    detected: set[str] = set()
    for case, tree, traits, trait in cases:
        report = fit_ornstein_uhlenbeck_model(tree, traits, trait=trait)
        warning_kinds = [warning.kind for warning in report.identifiability_warnings]
        detected.update(warning_kinds)
        observations.append(
            ComparativeOUIdentifiabilityCase(
                case=case,
                alpha=report.alpha,
                warning_kinds=warning_kinds,
            )
        )
    expected = [
        "small_sample_size",
        "boundary_alpha",
        "flat_likelihood",
        "weak_pull_to_optimum",
    ]
    return ComparativeOUIdentifiabilityAudit(
        cases=observations,
        detected_warning_kinds=sorted(detected),
        expected_warning_kinds=expected,
        all_expected_warning_kinds_detected=all(kind in detected for kind in expected),
    )
