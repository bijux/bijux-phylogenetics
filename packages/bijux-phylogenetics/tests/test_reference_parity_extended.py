from __future__ import annotations

import pytest

from bijux_phylogenetics.parity import validate_reference_parity_examples

pytestmark = [pytest.mark.scientific_validation, pytest.mark.slow]


def test_validate_reference_parity_examples_extended_passes() -> None:
    report = validate_reference_parity_examples(include_extended=True)
    assert report.all_passed is True
    assert report.case_count == 25
    assert report.failed_case_count == 0
    assert report.covered_methods == [
        "blombergs-k",
        "branch-score-distance",
        "brownian-trait-model",
        "consensus-tree-generation",
        "ornstein-uhlenbeck-trait-model",
        "pagels-lambda",
        "pgls",
        "phylogenetic-independent-contrasts",
        "posterior-clade-frequencies",
        "robinson-foulds-distance",
    ]
    posterior_row = next(
        row
        for row in report.summary_rows
        if row.method == "posterior-clade-frequencies"
    )
    pgls_row = next(row for row in report.summary_rows if row.method == "pgls")
    pagel_row = next(
        row for row in report.summary_rows if row.method == "pagels-lambda"
    )
    blomberg_row = next(
        row for row in report.summary_rows if row.method == "blombergs-k"
    )
    brownian_row = next(
        row for row in report.summary_rows if row.method == "brownian-trait-model"
    )
    ou_row = next(
        row
        for row in report.summary_rows
        if row.method == "ornstein-uhlenbeck-trait-model"
    )
    assert posterior_row.case_count == 2
    assert posterior_row.suite == "mixed"
    assert pgls_row.case_count == 5
    assert pgls_row.suite == "mixed"
    assert pagel_row.case_count == 4
    assert pagel_row.suite == "mixed"
    assert blomberg_row.case_count == 3
    assert blomberg_row.suite == "mixed"
    assert brownian_row.case_count == 2
    assert brownian_row.suite == "mixed"
    assert ou_row.case_count == 2
    assert ou_row.suite == "mixed"
    primate_estimated = next(
        row
        for row in report.observations
        if row.case == "pgls-primate-longevity-estimated-lambda"
    )
    primate_fixed = next(
        row
        for row in report.observations
        if row.case == "pgls-primate-longevity-fixed-reference-lambda"
    )
    assert primate_estimated.observed_output["aic"] <= 896.25
    assert primate_estimated.observed_output["lambda_value"] >= 0.75
    assert (
        primate_estimated.observed_output[
            "coefficient.social_group_size.standard_error"
        ]
        > 0.0
    )
    assert (
        primate_fixed.observed_output["aic"] < primate_estimated.observed_output["aic"]
    )
    assert primate_fixed.observed_output["lambda_value"] >= 0.75
    assert (
        primate_fixed.observed_output["coefficient.social_group_size.standard_error"]
        > 0.0
    )
    strong_lambda = next(
        row
        for row in report.observations
        if row.case == "pagel-lambda-non-ultrametric-strong-signal-twenty-four-taxa"
    )
    weak_lambda = next(
        row
        for row in report.observations
        if row.case == "pagel-lambda-weak-signal-twenty-four-taxa"
    )
    strong_k = next(
        row
        for row in report.observations
        if row.case == "blomberg-k-strong-signal-twenty-four-taxa"
    )
    weak_k = next(
        row
        for row in report.observations
        if row.case == "blomberg-k-weak-signal-twenty-four-taxa"
    )
    assert strong_lambda.observed_output["lambda_value"] >= 0.99
    assert weak_lambda.observed_output["lambda_value"] <= 0.001
    assert strong_k.observed_output["k"] > weak_k.observed_output["k"]
