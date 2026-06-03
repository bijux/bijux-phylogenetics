from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.likelihood import (
    evaluate_nucleotide_nested_likelihood_ratio_test_from_alignment,
    list_declared_nucleotide_likelihood_ratio_pairs,
    validate_declared_nucleotide_likelihood_ratio_pair,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_declared_nucleotide_likelihood_ratio_pairs_are_explicit() -> None:
    assert list_declared_nucleotide_likelihood_ratio_pairs() == [
        ("JC69", "K80"),
        ("JC69", "F81"),
        ("JC69", "HKY85"),
        ("JC69", "GTR"),
        ("K80", "HKY85"),
        ("K80", "GTR"),
        ("F81", "HKY85"),
        ("F81", "GTR"),
        ("HKY85", "GTR"),
    ]
    assert validate_declared_nucleotide_likelihood_ratio_pair("jc69", "k80") == (
        "JC69",
        "K80",
    )
    assert validate_declared_nucleotide_likelihood_ratio_pair("f81", "hky85") == (
        "F81",
        "HKY85",
    )


def test_nucleotide_nested_likelihood_ratio_report_governed_jc69_to_k80_case() -> None:
    report = evaluate_nucleotide_nested_likelihood_ratio_test_from_alignment(
        fixture("trees", "k80_kappa_optimization_tree_2_taxa.nwk"),
        fixture("alignments", "k80_kappa_optimization_alignment_2_taxa.fasta"),
        null_model_name="jc69",
        alternative_model_name="k80",
    )

    assert report.taxa == ["A", "B"]
    assert report.site_count == 13
    assert report.pattern_count == 9
    assert report.tree_newick == "(A:0.15,B:0.15);"
    assert report.null_fit.model_name == "JC69"
    assert report.alternative_fit.model_name == "K80"
    assert report.null_fit.parameter_count == 0
    assert report.alternative_fit.parameter_count == 1
    assert report.null_fit.parameter_values == {}
    assert math.isclose(
        report.alternative_fit.parameter_values["kappa"],
        9.17375541937712,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert report.likelihood_ratio_statistic > 0.0
    assert math.isclose(
        report.likelihood_ratio_statistic,
        2.0 * (report.alternative_fit.log_likelihood - report.null_fit.log_likelihood),
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.degrees_of_freedom == 1
    assert math.isclose(
        report.p_value,
        math.erfc(math.sqrt(report.likelihood_ratio_statistic / 2.0)),
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.p_value_method == "chi-square approximation"
    assert report.boundary_caveat == "no boundary parameter in declared nesting"
    assert report.null_fit.warnings == [
        "JC69 has no free substitution parameters; skipping parameter search"
    ]
    assert report.alternative_fit.warnings == []
    assert report.warnings == []


def test_nucleotide_nested_likelihood_ratio_supports_f81_to_hky85_pair() -> None:
    report = evaluate_nucleotide_nested_likelihood_ratio_test_from_alignment(
        fixture("trees", "hky85_kappa_optimization_tree_2_taxa.nwk"),
        fixture("alignments", "hky85_kappa_optimization_alignment_2_taxa.fasta"),
        null_model_name="f81",
        alternative_model_name="hky85",
    )

    assert report.null_fit.model_name == "F81"
    assert report.alternative_fit.model_name == "HKY85"
    assert report.degrees_of_freedom == 1
    assert report.likelihood_ratio_statistic > 0.0
    assert report.alternative_fit.log_likelihood > report.null_fit.log_likelihood


def test_nucleotide_nested_likelihood_ratio_rejects_nonnested_pair() -> None:
    with pytest.raises(
        ValueError,
        match="nested likelihood-ratio tests are only declared",
    ):
        evaluate_nucleotide_nested_likelihood_ratio_test_from_alignment(
            fixture("trees", "k80_kappa_optimization_tree_2_taxa.nwk"),
            fixture("alignments", "k80_kappa_optimization_alignment_2_taxa.fasta"),
            null_model_name="k80",
            alternative_model_name="f81",
        )
