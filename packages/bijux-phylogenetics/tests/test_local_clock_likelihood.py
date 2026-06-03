from __future__ import annotations

import math
from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    fit_local_clock_likelihood_from_alignment,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_likelihood_gateway_exports_local_clock_surface() -> None:
    assert (
        likelihood_api.fit_local_clock_likelihood_from_alignment
        is fit_local_clock_likelihood_from_alignment
    )


def test_local_clock_likelihood_fits_clade_and_branch_regimes_and_beats_strict_clock() -> (
    None
):
    report = fit_local_clock_likelihood_from_alignment(
        fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
        fixture("alignments", "local_clock_likelihood_alignment_4_taxa.fasta"),
        fixture("metadata", "local_clock_regimes_4_taxa.tsv"),
    )

    assert report.model_name == "JC69 local-clock"
    assert report.taxa == ["A", "B", "C", "D"]
    assert report.site_count == 60
    assert report.pattern_count == 4
    assert report.branch_count == 6
    assert report.regime_count == 3
    assert report.parameter_count == 3
    assert report.compression_used is True
    assert report.preferred_model_by_aic == "local-clock"
    assert report.aic < report.strict_clock_aic
    assert report.aic_delta_vs_strict_clock < 0.0
    assert report.optimized_log_likelihood > report.strict_clock_log_likelihood
    assert report.function_evaluation_count > report.optimization_pass_count
    assert report.converged is True

    regime_rows_by_id = {row.regime_id: row for row in report.regime_rows}
    assert set(regime_rows_by_id) == {"background", "ab_clade", "abc_stem"}
    assert regime_rows_by_id["background"].branch_count == 2
    assert regime_rows_by_id["ab_clade"].branch_count == 3
    assert regime_rows_by_id["abc_stem"].branch_count == 1
    assert regime_rows_by_id["background"].target_kind == "background"
    assert regime_rows_by_id["ab_clade"].target_kind == "clade"
    assert regime_rows_by_id["abc_stem"].target_kind == "branch"
    assert (
        regime_rows_by_id["abc_stem"].optimized_clock_rate
        > regime_rows_by_id["ab_clade"].optimized_clock_rate
    )
    assert (
        regime_rows_by_id["ab_clade"].optimized_clock_rate
        > regime_rows_by_id["background"].optimized_clock_rate
    )

    rows_by_descendant_taxa = {
        tuple(row.descendant_taxa): row for row in report.branch_rows
    }
    assert rows_by_descendant_taxa[("A", "B", "C")].regime_id == "abc_stem"
    assert rows_by_descendant_taxa[("A", "B")].regime_id == "ab_clade"
    assert rows_by_descendant_taxa[("A",)].regime_id == "ab_clade"
    assert rows_by_descendant_taxa[("B",)].regime_id == "ab_clade"
    assert rows_by_descendant_taxa[("C",)].regime_id == "background"
    assert rows_by_descendant_taxa[("D",)].regime_id == "background"
    assert math.isclose(
        rows_by_descendant_taxa[("A", "B", "C")].optimized_clock_rate,
        regime_rows_by_id["abc_stem"].optimized_clock_rate,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_local_clock_likelihood_rejects_ambiguous_overlapping_regimes() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="ambiguous",
    ):
        fit_local_clock_likelihood_from_alignment(
            fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
            fixture("alignments", "local_clock_likelihood_alignment_4_taxa.fasta"),
            fixture("metadata", "local_clock_conflicting_regimes_4_taxa.tsv"),
        )
