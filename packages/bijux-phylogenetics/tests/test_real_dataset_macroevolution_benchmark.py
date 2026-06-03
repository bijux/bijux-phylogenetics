from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark.macroevolution import (
    benchmark_real_dataset_macroevolution,
    run_real_dataset_macroevolution_benchmark_demo,
    write_real_dataset_macroevolution_bundle,
)
from tests.support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)


@pytest.mark.slow
def test_benchmark_real_dataset_macroevolution_reports_native_and_review_surfaces() -> (
    None
):
    report = benchmark_real_dataset_macroevolution()

    assert report.dataset.dataset_id == "central_european_seashore_flora"
    assert report.provenance_doi == "10.5061/dryad.0st06f0"
    assert len(report.summary_rows) == 4
    assert len(report.model_rows) == 8
    assert len(report.alignment_review_rows) == 2
    assert len(report.parity_rows) == 10

    summary_rows = {row.surface_id: row for row in report.summary_rows}
    assert summary_rows["seed-mass-native-model-table"].bijux_selected_model == (
        "ornstein-uhlenbeck"
    )
    assert summary_rows["seed-mass-native-model-table"].selection_matches_geiger is True
    assert (
        summary_rows["seed-mass-native-model-table"].stable_conclusion_supported
        is False
    )
    assert summary_rows["lifeform-native-model-table"].bijux_selected_model == (
        "equal-rates"
    )
    assert summary_rows["lifeform-native-model-table"].selection_matches_geiger is True
    assert (
        summary_rows["lifeform-native-model-table"].stable_conclusion_supported is True
    )

    alignment_rows = {row.surface_id: row for row in report.alignment_review_rows}
    continuous_review = alignment_rows["seed-mass-alignment-review"]
    assert continuous_review.aligned_taxa_count == 40
    assert continuous_review.dropped_tree_taxa == ["Triglochin_maritimum"]
    assert continuous_review.dropped_trait_taxa == ["unmatched_review_taxon"]
    assert continuous_review.dropped_missing_value_taxa == ["Juncus_maritimus"]

    discrete_review = alignment_rows["lifeform-alignment-review"]
    assert discrete_review.aligned_taxa_count == 40
    assert discrete_review.dropped_tree_taxa == ["Triglochin_maritimum"]
    assert discrete_review.dropped_trait_taxa == ["unmatched_review_taxon"]
    assert discrete_review.dropped_missing_value_taxa == ["Juncus_gerardii"]

    model_rows = {(row.surface_id, row.model): row for row in report.model_rows}
    assert (
        model_rows[
            ("seed-mass-native-model-table", "ornstein-uhlenbeck")
        ].bijux_selected
        is True
    )
    assert (
        model_rows[("lifeform-native-model-table", "equal-rates")].bijux_selected
        is True
    )
    assert any(
        "did not converge" in note
        for note in model_rows[("lifeform-native-model-table", "symmetric")].notes
    )
    assert any(
        "did not converge" in note
        for note in model_rows[
            ("lifeform-native-model-table", "all-rates-different")
        ].notes
    )

    parity_rows = {
        (row.surface_id, row.model, row.comparison_scope): row
        for row in report.parity_rows
    }
    assert (
        parity_rows[
            ("seed-mass-native-model-table", "ornstein-uhlenbeck", "native-model-table")
        ].within_aicc_tolerance
        is True
    )
    assert (
        parity_rows[
            ("lifeform-native-model-table", "all-rates-different", "native-model-table")
        ].within_aicc_tolerance
        is False
    )
    assert (
        parity_rows[
            ("seed-mass-alignment-review", "ornstein-uhlenbeck", "alignment-review")
        ].within_parameter_tolerance
        is True
    )


def test_public_runtime_exports_include_real_dataset_macroevolution_benchmark() -> None:
    assert (
        benchmark_api.benchmark_real_dataset_macroevolution
        is benchmark_real_dataset_macroevolution
    )
    assert (
        benchmark_api.run_real_dataset_macroevolution_benchmark_demo
        is run_real_dataset_macroevolution_benchmark_demo
    )
    assert (
        benchmark_api.write_real_dataset_macroevolution_bundle
        is write_real_dataset_macroevolution_bundle
    )


@pytest.mark.slow
def test_write_real_dataset_macroevolution_bundle_matches_expected_outputs(
    tmp_path: Path,
) -> None:
    bundle = write_real_dataset_macroevolution_bundle(
        tmp_path / "real-dataset-macroevolution"
    )

    expected_root = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "bijux_phylogenetics"
        / "resources"
        / "benchmarks"
        / "real_dataset_macroevolution"
        / "expected"
    )
    generated = {
        Path("benchmark-summary.tsv"): bundle.summary_path,
        Path("model-table.tsv"): bundle.model_table_path,
        Path("alignment-review.tsv"): bundle.alignment_review_path,
        Path("geiger-parity.tsv"): bundle.parity_table_path,
        Path("geiger-reference.tsv"): bundle.geiger_reference_path,
    }
    assert_selected_scientific_outputs_equivalent(expected_root, generated)
