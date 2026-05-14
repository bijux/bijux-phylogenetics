from __future__ import annotations

import pytest

from bijux_phylogenetics.reference_parity import validate_reference_parity_examples

pytestmark = [pytest.mark.scientific_validation]


def test_validate_reference_parity_examples_extended_passes() -> None:
    report = validate_reference_parity_examples(include_extended=True)
    assert report.all_passed is True
    assert report.case_count == 14
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
    assert posterior_row.case_count == 2
    assert posterior_row.suite == "mixed"
