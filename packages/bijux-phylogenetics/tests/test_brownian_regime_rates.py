from __future__ import annotations

import math
from pathlib import Path
import random

from bijux_phylogenetics.comparative.common import descendant_taxa, node_signature
from bijux_phylogenetics.comparative.continuous import (
    summarize_brownian_regime_rates,
    write_brownian_regime_branch_table,
    write_brownian_regime_comparison_table,
    write_brownian_regime_exclusion_table,
    write_brownian_regime_profile_table,
    write_brownian_regime_rate_table,
    write_brownian_regime_summary_table,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

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


def test_summarize_brownian_regime_rates_reports_branch_and_rate_ledgers() -> None:
    report = summarize_brownian_regime_rates(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        fixture("example_branch_regimes.tsv"),
        trait="response",
    )

    assert report.tree_taxon_count == 4
    assert report.analyzed_taxon_count == 4
    assert report.branch_id_column == "branch_id"
    assert report.regime_column == "regime"
    assert report.excluded_taxa == []
    assert [row.branch_id for row in report.branch_rows] == [
        "A",
        "A|B",
        "B",
        "C",
        "C|D",
        "D",
    ]
    assert [row.regime for row in report.regime_rows] == ["fast", "slow"]
    assert len(report.profile_rows) >= 162
    assert {row.regime for row in report.profile_rows if row.selected} == {
        "fast",
        "slow",
    }
    assert [row.model for row in report.comparison_rows] == [
        "brownian",
        "brownian-regimes",
    ]
    assert report.likelihood_ratio_degrees_of_freedom == 1
    assert 0.0 <= report.likelihood_ratio_p_value <= 1.0
    assert report.root_state_interval.name == "root_state"
    assert (
        report.root_state_interval.lower_95
        < report.root_state
        < report.root_state_interval.upper_95
    )
    assert all(
        row.lower_95 <= row.sigma_squared <= row.upper_95 for row in report.regime_rows
    )


def test_summarize_brownian_regime_rates_tracks_pruned_and_extra_taxa() -> None:
    report = summarize_brownian_regime_rates(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_brownian_missing.tsv"),
        fixture("example_branch_regimes_six_taxa.tsv"),
        trait="response_growth",
    )

    assert report.analyzed_taxa == ["A", "D", "E", "F"]
    assert {row.taxon: row.reason for row in report.excluded_taxa} == {
        "B": "missing_trait_value",
        "C": "non_numeric_trait_value",
        "G": "absent_from_tree",
    }
    branch_lookup = {row.branch_id: row for row in report.branch_rows}
    assert branch_lookup["A|B"].analyzed_descendant_taxa == ["A"]
    assert branch_lookup["C|D"].analyzed_descendant_taxa == ["D"]
    assert branch_lookup["E|F"].analyzed_descendant_taxa == ["E", "F"]
    assert (
        "one or more overlapping taxa have missing trait values and will be pruned"
        in report.warnings
    )
    assert (
        "one or more overlapping taxa have non-numeric trait values and will be pruned"
        in report.warnings
    )
    assert "trait table contains taxa absent from the tree" in report.warnings


def test_brownian_regime_rate_writers_emit_review_ledgers(tmp_path: Path) -> None:
    report = summarize_brownian_regime_rates(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_brownian_missing.tsv"),
        fixture("example_branch_regimes_six_taxa.tsv"),
        trait="response_growth",
    )
    summary_out = tmp_path / "brownian-regimes-summary.tsv"
    rates_out = tmp_path / "brownian-regimes-rates.tsv"
    profile_out = tmp_path / "brownian-regimes-profile.tsv"
    comparison_out = tmp_path / "brownian-regimes-comparison.tsv"
    branch_out = tmp_path / "brownian-regimes-branches.tsv"
    exclusion_out = tmp_path / "brownian-regimes-excluded.tsv"

    write_brownian_regime_summary_table(summary_out, report)
    write_brownian_regime_rate_table(rates_out, report)
    write_brownian_regime_profile_table(profile_out, report)
    write_brownian_regime_comparison_table(comparison_out, report)
    write_brownian_regime_branch_table(branch_out, report)
    write_brownian_regime_exclusion_table(exclusion_out, report)

    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("trait\ttaxon_column\tbranch_id_column\tregime_column")
    )
    assert (
        rates_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("regime\tbranch_count\tcontributing_branch_count")
    )
    assert (
        profile_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("regime\tsigma_squared\tlog_likelihood")
    )
    comparison_rows = comparison_out.read_text(encoding="utf-8").splitlines()
    assert comparison_rows[0].startswith("row_kind\tmodel\tcomparison_id")
    assert comparison_rows[1].startswith("model_fit\tbrownian\t")
    assert any(
        row.startswith("likelihood_ratio_test\t\tbrownian-vs-brownian-regimes\t")
        for row in comparison_rows[1:]
    )
    assert (
        branch_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("branch_id\tregime\tbranch_length")
    )
    assert exclusion_out.read_text(encoding="utf-8").splitlines() == [
        "taxon\treason",
        "B\tmissing_trait_value",
        "C\tnon_numeric_trait_value",
        "G\tabsent_from_tree",
    ]


def test_summarize_brownian_regime_rates_recovers_known_rate_order_on_simulated_data(
    tmp_path: Path,
) -> None:
    tree_path, traits_path, regime_map_path = _write_simulated_regime_dataset(tmp_path)

    report = summarize_brownian_regime_rates(
        tree_path,
        traits_path,
        regime_map_path,
        trait="trait",
    )

    rates = {row.regime: row.sigma_squared for row in report.regime_rows}
    assert report.better_model == "brownian-regimes"
    assert report.likelihood_ratio_statistic > 0.0
    assert rates["fast"] > rates["background"] > rates["slow"]
    assert any(row.regime == "fast" and row.selected for row in report.profile_rows)
    assert any(row.regime == "slow" and row.selected for row in report.profile_rows)


def _write_simulated_regime_dataset(tmp_path: Path) -> tuple[Path, Path, Path]:
    tree_path = tmp_path / "brownian-regime-tree.nwk"
    traits_path = tmp_path / "brownian-regime-traits.tsv"
    regime_map_path = tmp_path / "brownian-regimes.tsv"
    tree_path.write_text(
        (
            "((((A:0.8,B:0.4):0.7,(C:0.6,D:0.3):0.5):1.2,"
            "((E:0.7,F:0.35):0.8,(G:0.5,H:0.25):0.9):1.0):0.6,"
            "((I:0.9,J:0.45):0.6,(K:0.55,L:0.2):1.1):0.85);"
        ),
        encoding="utf-8",
    )
    tree = load_tree(tree_path)
    regime_by_branch = _simulated_regime_assignment(tree)
    regime_map_path.write_text(
        "branch_id\tregime\n"
        + "".join(
            f"{branch_id}\t{regime}\n" for branch_id, regime in regime_by_branch.items()
        ),
        encoding="utf-8",
    )
    covariance = _multirate_covariance(
        tree,
        regime_by_branch,
        regime_rates={
            "background": 1.0,
            "fast": 8.0,
            "slow": 0.1,
        },
    )
    cholesky = _cholesky(covariance)
    rng = random.Random(1)
    standard_normal = [rng.gauss(0.0, 1.0) for _ in tree.tip_names]
    values = [
        3.0
        + sum(
            cholesky[row][column] * standard_normal[column] for column in range(row + 1)
        )
        for row in range(len(tree.tip_names))
    ]
    traits_path.write_text(
        "taxon\ttrait\n"
        + "".join(
            f"{taxon}\t{value}\n"
            for taxon, value in zip(tree.tip_names, values, strict=True)
        ),
        encoding="utf-8",
    )
    return tree_path, traits_path, regime_map_path


def _simulated_regime_assignment(tree: PhyloTree) -> dict[str, str]:
    return {
        node_signature(node): (
            "fast"
            if node_signature(node) in {"A|B", "A", "B", "E|F", "E", "F"}
            else "slow"
            if node_signature(node) in {"I|J", "I", "J", "K|L", "K", "L"}
            else "background"
        )
        for node in tree.iter_nodes()
        if node is not tree.root
    }


def _multirate_covariance(
    tree: PhyloTree,
    regime_by_branch: dict[str, str],
    *,
    regime_rates: dict[str, float],
) -> list[list[float]]:
    taxon_index = {taxon: index for index, taxon in enumerate(tree.tip_names)}
    covariance = [[0.0] * len(tree.tip_names) for _ in tree.tip_names]
    for node in tree.iter_nodes():
        if node is tree.root:
            continue
        descendants = descendant_taxa(node)
        rate = regime_rates[regime_by_branch[node_signature(node)]]
        branch_length = float(node.branch_length or 0.0)
        for left_taxon in descendants:
            for right_taxon in descendants:
                covariance[taxon_index[left_taxon]][taxon_index[right_taxon]] += (
                    rate * branch_length
                )
    return covariance


def _cholesky(matrix: list[list[float]]) -> list[list[float]]:
    size = len(matrix)
    lower = [[0.0] * size for _ in range(size)]
    for row in range(size):
        for column in range(row + 1):
            partial = sum(
                lower[row][index] * lower[column][index] for index in range(column)
            )
            if row == column:
                lower[row][column] = math.sqrt(max(matrix[row][row] - partial, 1e-12))
            else:
                lower[row][column] = (matrix[row][column] - partial) / lower[column][
                    column
                ]
    return lower
