from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    estimate_pagels_lambda,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


@pytest.mark.parametrize(
    ("fixture_name", "case", "tree_name", "traits_name", "trait", "taxon_column"),
    [
        (
            "reference_parity_core.json",
            "blomberg-k-example-tree",
            "example_tree.nwk",
            "example_traits_comparative.tsv",
            "response",
            None,
        ),
        (
            "reference_parity_extended_comparative.json",
            "blomberg-k-strong-signal-twenty-four-taxa",
            "example_tree_phytools_ultrametric_twenty_four_taxa.nwk",
            "example_traits_phytools_signal_twenty_four_taxa.tsv",
            "signal_strong",
            "taxon",
        ),
        (
            "reference_parity_extended_comparative.json",
            "blomberg-k-weak-signal-twenty-four-taxa",
            "example_tree_phytools_ultrametric_twenty_four_taxa.nwk",
            "example_traits_phytools_signal_twenty_four_taxa.tsv",
            "signal_weak",
            "taxon",
        ),
    ],
)
def test_blombergs_k_matches_governed_phytools_reference(
    fixture_name: str,
    case: str,
    tree_name: str,
    traits_name: str,
    trait: str,
    taxon_column: str | None,
) -> None:
    observation = _reference_observation(fixture(fixture_name), case=case)
    report = compute_blombergs_k(
        fixture(tree_name),
        fixture(traits_name),
        trait=trait,
        taxon_column=taxon_column,
    )
    tolerance = float(observation["tolerance"])
    assert math.isclose(
        report.k,
        float(observation["expected_output"]["k"]),
        rel_tol=tolerance,
        abs_tol=tolerance,
    )


@pytest.mark.parametrize(
    ("fixture_name", "case", "tree_name", "traits_name", "trait", "taxon_column"),
    [
        (
            "reference_parity_core.json",
            "pagel-lambda-example-tree",
            "example_tree.nwk",
            "example_traits_comparative.tsv",
            "response",
            None,
        ),
        (
            "reference_parity_extended_comparative.json",
            "pagel-lambda-non-ultrametric-strong-signal-twenty-four-taxa",
            "example_tree_phytools_non_ultrametric_twenty_four_taxa.nwk",
            "example_traits_phytools_signal_non_ultrametric_twenty_four_taxa.tsv",
            "signal_strong",
            "taxon",
        ),
        (
            "reference_parity_extended_comparative.json",
            "pagel-lambda-weak-signal-twenty-four-taxa",
            "example_tree_phytools_ultrametric_twenty_four_taxa.nwk",
            "example_traits_phytools_signal_twenty_four_taxa.tsv",
            "signal_weak",
            "taxon",
        ),
    ],
)
def test_pagels_lambda_matches_governed_phytools_reference(
    fixture_name: str,
    case: str,
    tree_name: str,
    traits_name: str,
    trait: str,
    taxon_column: str | None,
) -> None:
    observation = _reference_observation(fixture(fixture_name), case=case)
    report = estimate_pagels_lambda(
        fixture(tree_name),
        fixture(traits_name),
        trait=trait,
        taxon_column=taxon_column,
    )
    tolerance = float(observation["tolerance"])
    assert math.isclose(
        report.lambda_value,
        float(observation["expected_output"]["lambda_value"]),
        rel_tol=tolerance,
        abs_tol=tolerance,
    )
    assert math.isclose(
        report.log_likelihood,
        float(observation["expected_output"]["log_likelihood"]),
        rel_tol=tolerance,
        abs_tol=tolerance,
    )


@pytest.mark.slow
def test_pagels_lambda_matches_governed_primate_reference() -> None:
    observation = _reference_observation(
        fixture("reference_parity_extended_comparative.json"),
        case="pagel-lambda-primate-longevity",
    )
    repository_root = Path(__file__).resolve().parents[3]
    report = estimate_pagels_lambda(
        repository_root
        / "evidence-book/studies/primate-longevity-signal/datasets/reference_trimmed_primatetree.nwk",
        repository_root
        / "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
        trait="longevity",
        taxon_column="species",
    )
    tolerance = float(observation["tolerance"])
    assert math.isclose(
        report.lambda_value,
        float(observation["expected_output"]["lambda_value"]),
        rel_tol=tolerance,
        abs_tol=tolerance,
    )
    assert math.isclose(
        report.log_likelihood,
        float(observation["expected_output"]["log_likelihood"]),
        rel_tol=tolerance,
        abs_tol=tolerance,
    )


def _reference_observation(path: Path, *, case: str) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return next(row for row in payload["observations"] if row["case"] == case)
