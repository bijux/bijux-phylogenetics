from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.simulation import (
    simulate_brownian_trait_collection,
    simulate_brownian_traits,
    simulate_speciational_trait_collection,
    simulate_speciational_traits,
    write_continuous_trait_collection_summary_table,
    write_continuous_trait_collection_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def _row_by_label(report, *, row_kind: str, label: str):
    return next(
        row for row in report.rows if row.row_kind == row_kind and row.label == label
    )


def test_simulate_brownian_traits_accepts_sigma_squared() -> None:
    sigma_report = simulate_brownian_traits(
        fixture("example_tree.nwk"),
        seed=7,
        root_state=1.0,
        sigma=0.5,
    )
    sigma_squared_report = simulate_brownian_traits(
        fixture("example_tree.nwk"),
        seed=7,
        root_state=1.0,
        sigma_squared=0.25,
    )

    assert sigma_squared_report.sigma_squared == 0.25
    assert sigma_squared_report.sigma == 0.5
    assert sigma_squared_report.traits == sigma_report.traits
    assert sigma_squared_report.node_values == sigma_report.node_values


def test_simulate_brownian_trait_collection_reports_tree_covariance_signal() -> None:
    report = simulate_brownian_trait_collection(
        fixture("example_tree.nwk"),
        root_state=0.0,
        sigma_squared=0.25,
        replicates=32,
        seed=7,
    )

    assert report.model == "brownian-motion"
    assert report.tip_count == 4
    assert report.branch_count == 6
    assert report.replicate_count == 32
    assert len(report.simulations) == 32
    assert len([row for row in report.rows if row.row_kind == "tip_distribution"]) == 4
    assert len([row for row in report.rows if row.row_kind == "tip_covariance"]) == 10
    assert report.simulations[0] == simulate_brownian_traits(
        fixture("example_tree.nwk"),
        root_state=0.0,
        sigma_squared=0.25,
        seed=7,
    )

    sister_covariance = _row_by_label(
        report,
        row_kind="tip_covariance",
        label="A|B",
    )
    distant_covariance = _row_by_label(
        report,
        row_kind="tip_covariance",
        label="A|C",
    )
    assert sister_covariance.covariance is not None
    assert distant_covariance.covariance is not None
    assert sister_covariance.covariance > distant_covariance.covariance


def test_brownian_collection_writers_and_cli_emit_sigma_squared(
    tmp_path: Path,
    capsys,
) -> None:
    report = simulate_brownian_trait_collection(
        fixture("example_tree.nwk"),
        root_state=1.0,
        sigma_squared=0.25,
        replicates=4,
        seed=7,
    )
    collection_path = write_continuous_trait_collection_table(
        tmp_path / "brownian-collection.tsv",
        report,
    )
    summary_path = write_continuous_trait_collection_summary_table(
        tmp_path / "brownian-summary.tsv",
        report,
    )

    assert collection_path.read_text(encoding="utf-8").splitlines()[0] == (
        "replicate_index\ttaxon\tvalue"
    )
    assert summary_path.read_text(encoding="utf-8").splitlines()[0] == (
        "row_kind\tlabel\tmean_value\tstandard_deviation\tminimum\tmedian\tmaximum\tcovariance\tcorrelation"
    )

    output = tmp_path / "brownian.tsv"
    exit_code = main(
        [
            "simulate",
            "traits-brownian",
            str(fixture("example_tree.nwk")),
            "--root-state",
            "1.0",
            "--sigma-squared",
            "0.25",
            "--seed",
            "7",
            "--out",
            str(output),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["sigma_squared"] == 0.25
    assert payload["data"]["sigma_squared"] == 0.25
    assert output.read_text(encoding="utf-8").splitlines()[0] == "taxon\tvalue"


def test_simulate_speciational_traits_equalizes_positive_branch_lengths() -> None:
    report = simulate_speciational_traits(
        fixture("example_tree_internal_long_branch.nwk"),
        root_state=1.5,
        sigma_squared=0.25,
        seed=11,
    )

    assert report.model == "speciational"
    assert report.sigma_squared == 0.25
    node_value_by_taxon = {row.taxon: row.value for row in report.traits}
    assert node_value_by_taxon["A"] != node_value_by_taxon["B"]
    assert node_value_by_taxon["C"] != node_value_by_taxon["D"]


def test_speciational_collection_uses_branch_count_not_branch_magnitude() -> None:
    brownian_report = simulate_brownian_trait_collection(
        fixture("example_tree_internal_long_branch.nwk"),
        root_state=0.0,
        sigma_squared=0.25,
        replicates=64,
        seed=11,
    )
    speciational_report = simulate_speciational_trait_collection(
        fixture("example_tree_internal_long_branch.nwk"),
        root_state=0.0,
        sigma_squared=0.25,
        replicates=64,
        seed=11,
    )

    brownian_tip_sd = {
        row.label: row.standard_deviation
        for row in brownian_report.rows
        if row.row_kind == "tip_distribution"
    }
    speciational_tip_sd = {
        row.label: row.standard_deviation
        for row in speciational_report.rows
        if row.row_kind == "tip_distribution"
    }

    assert brownian_report.model == "brownian-motion"
    assert speciational_report.model == "speciational"
    assert brownian_tip_sd["A"] is not None
    assert brownian_tip_sd["C"] is not None
    assert speciational_tip_sd["A"] is not None
    assert speciational_tip_sd["C"] is not None
    assert brownian_tip_sd["A"] > brownian_tip_sd["C"]
    assert abs(speciational_tip_sd["A"] - speciational_tip_sd["C"]) < 0.15
