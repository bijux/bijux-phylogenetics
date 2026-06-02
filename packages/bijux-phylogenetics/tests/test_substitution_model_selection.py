from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.likelihood import (
    compare_nucleotide_substitution_models_from_alignment,
    default_substitution_model_selection_candidates,
)

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_default_substitution_model_selection_candidates_cover_supported_models() -> (
    None
):
    assert default_substitution_model_selection_candidates() == (
        "JC69",
        "JC69+G",
        "JC69+I",
        "JC69+G+I",
        "K80",
        "K80+G",
        "K80+I",
        "K80+G+I",
        "F81",
        "F81+G",
        "F81+I",
        "F81+G+I",
        "HKY85",
        "HKY85+G",
        "HKY85+I",
        "HKY85+G+I",
        "GTR",
        "GTR+G",
        "GTR+I",
        "GTR+G+I",
    )


def test_substitution_model_selection_report_emits_ranked_rows_and_failed_candidates() -> (
    None
):
    report = compare_nucleotide_substitution_models_from_alignment(
        fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta"),
        candidate_models=(
            "JC69",
            "JC69+G",
            "K80+I",
            "F81",
            "HKY",
            "GTR",
            "NOT-A-MODEL",
        ),
        max_coordinate_passes=2,
    )

    assert report.taxa == ["A", "B"]
    assert report.site_count == 4
    assert report.pattern_count == 4
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert report.best_model_aic == "JC69"
    assert report.best_model_aicc == "JC69"
    assert report.best_model_bic == "JC69"
    assert len(report.rows) == 7

    rows = {row.model_name: row for row in report.rows}
    assert set(rows) == {
        "JC69",
        "JC69+G",
        "K80+I",
        "F81",
        "HKY85",
        "GTR",
        "NOT-A-MODEL",
    }
    assert "HKY" not in rows

    assert math.isclose(
        sum(row.akaike_weight or 0.0 for row in report.rows if row.fit_succeeded),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    jc69_row = rows["JC69"]
    assert jc69_row.fit_succeeded is True
    assert jc69_row.parameter_count == 0
    assert math.isclose(
        jc69_row.log_likelihood or 0.0,
        -11.105102730824235,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        jc69_row.aic or 0.0, 22.21020546164847, rel_tol=0.0, abs_tol=1e-12
    )
    assert math.isclose(
        jc69_row.aicc or 0.0, 22.21020546164847, rel_tol=0.0, abs_tol=1e-12
    )
    assert math.isclose(
        jc69_row.bic or 0.0, 22.21020546164847, rel_tol=0.0, abs_tol=1e-12
    )
    assert math.isclose(jc69_row.delta_aic or 0.0, 0.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        jc69_row.akaike_weight or 0.0,
        0.622367926223253,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert jc69_row.rank == 1
    assert jc69_row.selected_by_aic is True
    assert jc69_row.selected_by_aicc is True
    assert jc69_row.selected_by_bic is True
    assert jc69_row.warnings == [
        "JC69 has no free substitution parameters; skipping parameter search"
    ]

    jc69_gamma_row = rows["JC69+G"]
    assert jc69_gamma_row.parameter_count == 1
    assert math.isclose(
        jc69_gamma_row.delta_aic or 0.0,
        2.01553984612546,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert jc69_gamma_row.warnings == ["alpha hit upper search boundary"]

    k80_invariant_row = rows["K80+I"]
    assert k80_invariant_row.parameter_count == 2
    assert math.isclose(
        k80_invariant_row.aicc or 0.0,
        37.99276321667401,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert k80_invariant_row.rank == 3
    assert k80_invariant_row.warnings == [
        "invariant_proportion hit lower search boundary"
    ]

    f81_row = rows["F81"]
    assert f81_row.parameter_count == 3
    assert f81_row.aicc is None
    assert f81_row.comparable_on_aicc is False
    assert f81_row.warnings == [
        "sample size is too small to compute finite AICc for this parameter count"
    ]

    hky85_row = rows["HKY85"]
    assert hky85_row.parameter_count == 4
    assert hky85_row.rank == 5
    assert hky85_row.aicc is None
    assert hky85_row.warnings == [
        "sample size is too small to compute finite AICc for this parameter count"
    ]

    gtr_row = rows["GTR"]
    assert gtr_row.parameter_count == 8
    assert gtr_row.rank == 6
    assert gtr_row.aicc is None
    assert gtr_row.comparable_on_aicc is False
    assert gtr_row.warnings == [
        "AT hit lower search boundary",
        "CG hit lower search boundary",
        "CT hit lower search boundary",
        "GT hit upper search boundary",
        "sample size is too small to compute finite AICc for this parameter count",
    ]

    failed_row = rows["NOT-A-MODEL"]
    assert failed_row.fit_succeeded is False
    assert failed_row.parameter_count is None
    assert failed_row.log_likelihood is None
    assert failed_row.aic is None
    assert failed_row.aicc is None
    assert failed_row.bic is None
    assert failed_row.delta_aic is None
    assert failed_row.akaike_weight is None
    assert failed_row.rank is None
    assert failed_row.comparable_on_aic is False
    assert failed_row.comparable_on_aicc is False
    assert failed_row.comparable_on_bic is False
    assert failed_row.warnings == [
        "candidate substitution model must use one of JC69, K80, F81, HKY85, GTR"
    ]

    assert report.warnings == [
        "one or more substitution candidates failed but remain visible in the table: NOT-A-MODEL"
    ]
