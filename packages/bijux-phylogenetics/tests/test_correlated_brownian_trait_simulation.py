from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.simulation import (
    simulate_correlated_brownian_trait_collection,
    simulate_correlated_brownian_traits,
    write_correlated_continuous_trait_collection_summary_table,
    write_correlated_continuous_trait_collection_table,
    write_correlated_continuous_trait_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def _row_by_label(report, *, row_kind: str, label: str):
    return next(
        row for row in report.rows if row.row_kind == row_kind and row.label == label
    )


def test_simulate_correlated_brownian_traits_accepts_correlation_matrix() -> None:
    report = simulate_correlated_brownian_traits(
        fixture("example_tree.nwk"),
        trait_names=["trait_alpha", "trait_beta"],
        evolutionary_correlation_matrix=[
            [1.0, 0.4],
            [0.4, 1.0],
        ],
        trait_standard_deviations=[2.0, 1.5],
        root_states=[1.0, -0.5],
        seed=7,
    )

    assert report.model == "multivariate-brownian-motion"
    assert report.tip_count == 4
    assert report.trait_names == ["trait_alpha", "trait_beta"]
    assert report.root_states == [1.0, -0.5]
    assert report.evolutionary_covariance_matrix == [
        [4.0, 1.2],
        [1.2, 2.25],
    ]
    assert len(report.traits) == 8
    assert report.traits[0].taxon == "A"
    assert report.traits[0].trait == "trait_alpha"


def test_simulate_correlated_brownian_traits_rejects_invalid_covariance() -> None:
    with pytest.raises(
        ValueError,
        match="positive-definite evolutionary covariance matrix",
    ):
        simulate_correlated_brownian_traits(
            fixture("example_tree.nwk"),
            trait_names=["trait_alpha", "trait_beta"],
            evolutionary_covariance_matrix=[
                [1.0, 2.0],
                [2.0, 1.0],
            ],
            seed=7,
        )


def test_correlated_brownian_collection_reports_tree_and_trait_covariance(
    tmp_path: Path,
) -> None:
    report = simulate_correlated_brownian_trait_collection(
        fixture("example_tree.nwk"),
        trait_names=["trait_alpha", "trait_beta"],
        evolutionary_covariance_matrix=[
            [0.8, 0.3],
            [0.3, 0.4],
        ],
        root_states=[0.0, 1.0],
        replicates=64,
        seed=11,
    )

    assert report.model == "multivariate-brownian-motion"
    assert report.tip_count == 4
    assert report.branch_count == 6
    assert report.trait_names == ["trait_alpha", "trait_beta"]
    assert report.root_states == [0.0, 1.0]
    assert len(report.simulations) == 64
    assert len([row for row in report.rows if row.row_kind == "root_state"]) == 2
    assert (
        len([row for row in report.rows if row.row_kind == "evolutionary_covariance"])
        == 3
    )
    assert len([row for row in report.rows if row.row_kind == "tip_distribution"]) == 8
    assert len([row for row in report.rows if row.row_kind == "tip_covariance"]) == 36

    same_taxon_cross_trait = _row_by_label(
        report,
        row_kind="tip_covariance",
        label="A|trait_alpha||A|trait_beta",
    )
    sister_same_trait = _row_by_label(
        report,
        row_kind="tip_covariance",
        label="A|trait_alpha||B|trait_alpha",
    )
    distant_same_trait = _row_by_label(
        report,
        row_kind="tip_covariance",
        label="A|trait_alpha||C|trait_alpha",
    )
    covariance_parameter = _row_by_label(
        report,
        row_kind="evolutionary_covariance",
        label="trait_alpha|trait_beta",
    )

    assert same_taxon_cross_trait.covariance is not None
    assert same_taxon_cross_trait.covariance > 0.0
    assert sister_same_trait.covariance is not None
    assert distant_same_trait.covariance is not None
    assert sister_same_trait.covariance > distant_same_trait.covariance
    assert covariance_parameter.covariance == 0.3
    assert covariance_parameter.correlation is not None
    assert covariance_parameter.correlation > 0.5

    collection_path = write_correlated_continuous_trait_collection_table(
        tmp_path / "correlated-brownian-collection.tsv",
        report,
    )
    summary_path = write_correlated_continuous_trait_collection_summary_table(
        tmp_path / "correlated-brownian-summary.tsv",
        report,
    )
    single_path = write_correlated_continuous_trait_table(
        tmp_path / "correlated-brownian-single.tsv",
        report.simulations[0],
    )

    assert collection_path.read_text(encoding="utf-8").splitlines()[0] == (
        "replicate_index\ttaxon\ttrait\tvalue"
    )
    assert summary_path.read_text(encoding="utf-8").splitlines()[0] == (
        "row_kind\tlabel\tmean_value\tstandard_deviation\tminimum\tmedian\tmaximum\tcovariance\tcorrelation"
    )
    assert single_path.read_text(encoding="utf-8").splitlines()[0] == (
        "taxon\ttrait\tvalue"
    )
