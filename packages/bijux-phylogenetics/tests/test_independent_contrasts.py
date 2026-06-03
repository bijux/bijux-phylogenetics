from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.comparative.common import load_comparative_dataset
from bijux_phylogenetics.comparative.signal import (
    compute_phylogenetic_independent_contrasts,
    compute_phylogenetic_independent_contrasts_from_dataset,
    summarize_independent_contrast_regression,
    write_independent_contrast_regression_table,
    write_independent_contrast_table,
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


def test_independent_contrast_report_matches_reference_fixture_case() -> None:
    report = compute_phylogenetic_independent_contrasts(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    reference = json.loads(
        fixture("comparative_reference_validation.json").read_text(encoding="utf-8")
    )
    expected = next(
        case for case in reference["observations"] if case["case"] == "pic-example-tree"
    )
    observed = {row.node: row.contrast for row in report.contrasts}
    observed_node_ids = {row.node: row.node_id for row in report.contrasts}
    for node, value in expected["expected_parameters"].items():
        assert math.isclose(observed[node], value, rel_tol=1e-12, abs_tol=1e-12)
    assert observed_node_ids == {"A|B": 6, "C|D": 7, "A|B|C|D": 5}


def test_independent_contrasts_from_dataset_matches_path_surface() -> None:
    dataset = load_comparative_dataset(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        require_binary=True,
        minimum_taxa=2,
    )

    direct_report = compute_phylogenetic_independent_contrasts_from_dataset(dataset)
    path_report = compute_phylogenetic_independent_contrasts(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )

    assert direct_report.tree_path == path_report.tree_path
    assert direct_report.traits_path == path_report.traits_path
    assert direct_report.trait == path_report.trait
    assert direct_report.taxon_count == path_report.taxon_count
    assert direct_report.input_audit == path_report.input_audit
    assert math.isclose(
        direct_report.root_estimate,
        path_report.root_estimate,
        rel_tol=1e-12,
        abs_tol=1e-12,
    )
    assert direct_report.contrasts == path_report.contrasts


def test_independent_contrast_regression_supports_origin_fit() -> None:
    report = summarize_independent_contrast_regression(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response_trait="response",
        predictor_trait="predictor_one",
    )
    assert report.contrast_count == 3
    assert len(report.rows) == 3
    assert math.isclose(report.slope, 0.9576271186440678, rel_tol=1e-12, abs_tol=1e-12)
    assert math.isclose(
        report.r_squared_through_origin,
        0.7869953775038521,
        rel_tol=1e-12,
        abs_tol=1e-12,
    )
    assert report.lower_95_confidence_interval < report.slope
    assert report.upper_95_confidence_interval > report.slope


def test_write_independent_contrast_tables_write_review_rows(tmp_path: Path) -> None:
    contrast_out = tmp_path / "independent-contrasts.tsv"
    regression_out = tmp_path / "independent-contrast-regression.tsv"
    contrast_report = compute_phylogenetic_independent_contrasts(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    regression_report = summarize_independent_contrast_regression(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response_trait="response",
        predictor_trait="predictor_one",
    )
    write_independent_contrast_table(contrast_out, contrast_report)
    write_independent_contrast_regression_table(regression_out, regression_report)
    contrast_rows = contrast_out.read_text(encoding="utf-8").splitlines()
    regression_rows = regression_out.read_text(encoding="utf-8").splitlines()
    assert contrast_rows[0].startswith("trait\tnode_id\tnode\tleft_taxa\tright_taxa")
    assert regression_rows[0].startswith("response_trait\tpredictor_trait\tnode")
    assert len(contrast_rows) == 4
    assert len(regression_rows) == 4
