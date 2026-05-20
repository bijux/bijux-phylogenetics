from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

from bijux_phylogenetics.comparative.continuous import (
    fit_brownian_motion_model,
    fit_ornstein_uhlenbeck_model,
)
from bijux_phylogenetics.comparative.pgls import run_pgls
from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    compute_phylogenetic_independent_contrasts,
    estimate_pagels_lambda,
)


@dataclass(slots=True)
class ComparativeReferenceObservation:
    """One deterministic comparative-model validation example."""

    case: str
    model: str
    trait: str
    source: str
    expected_parameters: dict[str, float]
    observed_parameters: dict[str, float]
    tolerance: float
    passed: bool


@dataclass(slots=True)
class ComparativeReferenceValidationReport:
    """Validation of standalone comparative models against durable reference examples."""

    observations: list[ComparativeReferenceObservation]
    all_passed: bool


def validate_comparative_reference_examples() -> ComparativeReferenceValidationReport:
    """Validate comparative-model outputs against checked-in external reference expectations."""
    fixtures = _load_comparative_reference_fixture()
    root = Path(__file__).resolve().parents[4] / "tests/fixtures"
    tree = root / "trees/example_tree.nwk"
    traits = root / "metadata/example_traits_comparative.tsv"

    contrasts = compute_phylogenetic_independent_contrasts(
        tree, traits, trait="response"
    )
    contrast_lookup = {row.node: row.contrast for row in contrasts.contrasts}
    pgls = run_pgls(
        tree,
        traits,
        response="response",
        predictors=["predictor_one"],
        lambda_value=1.0,
    )
    brownian = fit_brownian_motion_model(tree, traits, trait="response")
    ou = fit_ornstein_uhlenbeck_model(tree, traits, trait="response")
    blomberg = compute_blombergs_k(tree, traits, trait="response")
    pagel_lambda = estimate_pagels_lambda(tree, traits, trait="response")

    observed_by_case = {
        "brownian-example-tree": {
            "root_state": brownian.root_state,
            "rate": brownian.rate,
            "log_likelihood": brownian.log_likelihood,
        },
        "ou-example-tree-grid": {
            "alpha": ou.alpha,
            "theta": ou.theta,
            "log_likelihood": ou.log_likelihood,
        },
        "pic-example-tree": {
            node: contrast_lookup[node] for node in ("A|B", "C|D", "A|B|C|D")
        },
        "pgls-example-tree-brownian": {
            "intercept": pgls.coefficients[0].estimate,
            "predictor_one": pgls.coefficients[1].estimate,
            "log_likelihood": pgls.log_likelihood,
            **{
                f"residual_{taxon}": residual
                for taxon, residual in zip(pgls.taxa, pgls.residuals, strict=True)
            },
        },
        "blomberg-k-example-tree": {
            "k": blomberg.k,
        },
        "pagel-lambda-example-tree": {
            "lambda_value": pagel_lambda.lambda_value,
            "log_likelihood": pagel_lambda.log_likelihood,
        },
    }

    observations: list[ComparativeReferenceObservation] = []
    for example in fixtures["observations"]:
        observed = observed_by_case[example["case"]]
        tolerance = float(example["tolerance"])
        passed = all(
            math.isclose(
                observed[name], expected_value, rel_tol=tolerance, abs_tol=tolerance
            )
            for name, expected_value in example["expected_parameters"].items()
        )
        observations.append(
            ComparativeReferenceObservation(
                case=example["case"],
                model=example["model"],
                trait=example["trait"],
                source=example["source"],
                expected_parameters={
                    name: float(value)
                    for name, value in example["expected_parameters"].items()
                },
                observed_parameters={
                    name: float(observed[name])
                    for name in example["expected_parameters"]
                },
                tolerance=tolerance,
                passed=passed,
            )
        )
    return ComparativeReferenceValidationReport(
        observations=observations,
        all_passed=all(observation.passed for observation in observations),
    )


def _load_comparative_reference_fixture() -> dict[str, object]:
    fixture_path = (
        Path(__file__).resolve().parents[4]
        / "tests/fixtures/expected/comparative_reference_validation.json"
    )
    return json.loads(fixture_path.read_text(encoding="utf-8"))
