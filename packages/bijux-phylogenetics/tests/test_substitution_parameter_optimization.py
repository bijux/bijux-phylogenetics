from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    optimize_nucleotide_substitution_parameters_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_nucleotide_substitution_parameter_recovery_benchmark_reports_governed_models() -> (
    None
):
    jc69_report = optimize_nucleotide_substitution_parameters_from_alignment(
        fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta"),
        model_name="jc69",
    )

    assert jc69_report.model_name == "JC69"
    assert jc69_report.parameter_count == 0
    assert jc69_report.parameter_rows == []
    assert jc69_report.fixed_parameter_values == {}
    assert jc69_report.optimization_pass_count == 0
    assert jc69_report.function_evaluation_count == 1
    assert jc69_report.converged is True
    assert jc69_report.warnings == [
        "JC69 has no free substitution parameters; skipping parameter search"
    ]
    assert math.isclose(
        jc69_report.initial_log_likelihood,
        jc69_report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        jc69_report.initial_aic,
        jc69_report.optimized_aic,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    k80_report = optimize_nucleotide_substitution_parameters_from_alignment(
        fixture("trees", "k80_kappa_optimization_tree_2_taxa.nwk"),
        fixture("alignments", "k80_kappa_optimization_alignment_2_taxa.fasta"),
        model_name="k80",
        initial_kappa=1.0,
    )
    k80_rows = {row.parameter_name: row for row in k80_report.parameter_rows}

    assert k80_report.model_name == "K80"
    assert k80_report.tree_newick == "(A:0.15,B:0.15);"
    assert k80_report.site_count == 13
    assert k80_report.pattern_count == 9
    assert k80_report.parameter_count == 1
    assert k80_report.optimization_pass_count == 1
    assert k80_report.converged is True
    assert k80_report.warnings == []
    assert math.isclose(
        k80_rows["kappa"].initial_value,
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        k80_rows["kappa"].optimized_value,
        9.17375541937712,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert math.isclose(k80_rows["kappa"].lower_bound, 0.05, rel_tol=0.0, abs_tol=0.0)
    assert math.isclose(k80_rows["kappa"].upper_bound, 20.0, rel_tol=0.0, abs_tol=0.0)
    assert k80_rows["kappa"].hit_lower_bound is False
    assert k80_rows["kappa"].hit_upper_bound is False
    assert k80_report.optimized_log_likelihood > k80_report.initial_log_likelihood

    hky85_report = optimize_nucleotide_substitution_parameters_from_alignment(
        fixture("trees", "hky85_kappa_optimization_tree_2_taxa.nwk"),
        fixture("alignments", "hky85_kappa_optimization_alignment_2_taxa.fasta"),
        model_name="hky85",
        initial_kappa=1.0,
    )
    hky85_rows = {row.parameter_name: row for row in hky85_report.parameter_rows}

    assert hky85_report.model_name == "HKY85"
    assert hky85_report.tree_newick == "(A:0.15,B:0.15);"
    assert hky85_report.site_count == 128
    assert hky85_report.pattern_count == 16
    assert hky85_report.parameter_count == 4
    assert hky85_report.base_frequency_source == "estimated"
    assert hky85_report.optimization_pass_count == 1
    assert hky85_report.converged is True
    assert hky85_report.warnings == []
    assert math.isclose(
        hky85_rows["kappa"].initial_value,
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        hky85_rows["kappa"].optimized_value,
        5.55724924899235,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert hky85_rows["kappa"].hit_lower_bound is False
    assert hky85_rows["kappa"].hit_upper_bound is False
    assert hky85_report.optimized_log_likelihood > hky85_report.initial_log_likelihood

    gtr_report = optimize_nucleotide_substitution_parameters_from_alignment(
        fixture("trees", "gtr_exchangeability_optimization_tree_2_taxa.nwk"),
        fixture("alignments", "gtr_exchangeability_optimization_alignment_2_taxa.fasta"),
        model_name="gtr",
    )
    gtr_rows = {row.parameter_name: row for row in gtr_report.parameter_rows}

    assert gtr_report.model_name == "GTR"
    assert gtr_report.tree_newick == "(A:0.15,B:0.15);"
    assert gtr_report.site_count == 256
    assert gtr_report.pattern_count == 15
    assert gtr_report.parameter_count == 8
    assert gtr_report.base_frequency_source == "estimated"
    assert gtr_report.fixed_parameter_values == {"AC": 1.0}
    assert gtr_report.optimization_pass_count == 24
    assert gtr_report.converged is True
    assert gtr_report.warnings == []
    assert math.isclose(
        gtr_rows["AG"].initial_value,
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        gtr_rows["AG"].optimized_value,
        10.436093496916191,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert gtr_rows["AG"].optimized_value > gtr_rows["CT"].optimized_value
    assert gtr_rows["AG"].optimized_value > gtr_rows["GT"].optimized_value
    assert gtr_rows["AG"].optimized_value > gtr_rows["CG"].optimized_value
    assert gtr_rows["AG"].optimized_value > gtr_rows["AT"].optimized_value
    assert gtr_report.optimized_log_likelihood > gtr_report.initial_log_likelihood
