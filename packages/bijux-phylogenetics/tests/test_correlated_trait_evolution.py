from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.comparative.trait_dependence import (
    summarize_correlated_trait_evolution,
    write_correlated_trait_comparison_table,
    write_correlated_trait_exclusion_table,
    write_correlated_trait_observation_table,
    write_correlated_trait_summary_table,
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


def test_summarize_correlated_trait_evolution_reports_continuous_coupling() -> None:
    report = summarize_correlated_trait_evolution(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        left_trait="response",
        right_trait="predictor_one",
    )
    reference = json.loads(
        fixture("comparative_reference_validation.json").read_text(encoding="utf-8")
    )
    expected = next(
        case for case in reference["observations"] if case["case"] == "pic-example-tree"
    )
    observation_lookup = {row.label: row for row in report.observation_rows}
    assert report.analysis_kind == "continuous-brownian-contrasts"
    assert report.analyzed_taxa == ["A", "B", "C", "D"]
    assert report.association_measure_name == "evolutionary_correlation"
    assert math.isclose(report.evolutionary_covariance, 5.231481481481482)
    assert math.isclose(report.evolutionary_correlation, 0.8871275993361114)
    assert report.better_model == "correlated"
    assert report.likelihood_ratio_statistic > 0.0
    for node, value in expected["expected_parameters"].items():
        assert math.isclose(
            observation_lookup[node].left_numeric_value,
            value,
            rel_tol=1e-12,
            abs_tol=1e-12,
        )


def test_summarize_correlated_trait_evolution_reports_binary_coupling() -> None:
    report = summarize_correlated_trait_evolution(
        fixture("example_tree_eight_taxa.nwk"),
        fixture("example_traits_correlated_binary.tsv"),
        left_trait="trait_left",
        right_trait="trait_right",
    )
    assert report.analysis_kind == "binary-joint-state"
    assert report.analyzed_taxa == ["A", "B", "C", "D", "E", "F", "G", "H"]
    assert report.left_state_order == ["0", "1"]
    assert report.right_state_order == ["0", "1"]
    assert report.joint_state_counts == {"00": 3, "01": 1, "10": 1, "11": 3}
    assert report.association_measure_name == "phi_correlation"
    assert math.isclose(report.association_measure_value, 0.5)
    assert math.isclose(report.evolutionary_covariance, 0.125)
    assert math.isclose(report.evolutionary_correlation, 0.5)
    assert report.better_model in {"independent", "correlated"}
    assert len(report.comparison_rows) == 2
    assert len(report.observation_rows) == 8


def test_summarize_correlated_trait_evolution_tracks_excluded_taxa() -> None:
    report = summarize_correlated_trait_evolution(
        fixture("example_tree_eight_taxa.nwk"),
        fixture("example_traits_correlated_binary_missing.tsv"),
        left_trait="trait_left",
        right_trait="trait_right",
    )
    assert report.analyzed_taxa == ["A", "C", "D", "E", "F"]
    assert {row.taxon: row.reason for row in report.excluded_taxa} == {
        "B": "missing_trait_value",
        "G": "missing_from_trait_table",
        "H": "missing_from_trait_table",
        "I": "missing_from_tree",
    }


def test_write_correlated_trait_tables_write_review_ledgers(tmp_path: Path) -> None:
    report = summarize_correlated_trait_evolution(
        fixture("example_tree_eight_taxa.nwk"),
        fixture("example_traits_correlated_binary.tsv"),
        left_trait="trait_left",
        right_trait="trait_right",
    )
    summary_out = tmp_path / "correlated-traits-summary.tsv"
    comparison_out = tmp_path / "correlated-traits-comparison.tsv"
    observation_out = tmp_path / "correlated-traits-observations.tsv"
    exclusion_out = tmp_path / "correlated-traits-excluded.tsv"
    write_correlated_trait_summary_table(summary_out, report)
    write_correlated_trait_comparison_table(comparison_out, report)
    write_correlated_trait_observation_table(observation_out, report)
    write_correlated_trait_exclusion_table(exclusion_out, report)
    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    comparison_rows = comparison_out.read_text(encoding="utf-8").splitlines()
    observation_rows = observation_out.read_text(encoding="utf-8").splitlines()
    exclusion_rows = exclusion_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("analysis_kind\tleft_trait\tright_trait")
    assert comparison_rows[0].startswith(
        "model_kind\tmodel_description\tparameter_count"
    )
    assert observation_rows[0].startswith("row_kind\tlabel\ttaxon\tleft_taxa")
    assert exclusion_rows == ["taxon\treason\tmissing_traits"]
