from __future__ import annotations

import json
import math
from pathlib import Path

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


def test_blombergs_k_matches_core_phytools_reference() -> None:
    observation = _reference_observation(
        fixture("reference_parity_core.json"),
        case="blomberg-k-example-tree",
    )
    report = compute_blombergs_k(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    tolerance = float(observation["tolerance"])
    assert math.isclose(
        report.k,
        float(observation["expected_output"]["k"]),
        rel_tol=tolerance,
        abs_tol=tolerance,
    )


def test_pagels_lambda_matches_core_phytools_reference() -> None:
    observation = _reference_observation(
        fixture("reference_parity_core.json"),
        case="pagel-lambda-example-tree",
    )
    report = estimate_pagels_lambda(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
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
