from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
    evaluate_named_empirical_protein_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    evaluate_compiled_finite_state_likelihood_analysis,
    evaluate_compiled_finite_state_site_log_likelihood_rows,
)
from bijux_phylogenetics.phylo.likelihood.protein_finite_state import (
    build_named_empirical_protein_finite_state_surface,
    compile_named_empirical_protein_likelihood_analysis,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_named_empirical_protein_surface_and_compiled_analysis_stay_governed() -> None:
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta")
    )
    surface = build_named_empirical_protein_finite_state_surface(
        records,
        model_name="WAG",
    )
    analysis = compile_named_empirical_protein_likelihood_analysis(
        tree,
        records,
        model_name="WAG",
    )
    report = evaluate_named_empirical_protein_tree_likelihood(
        tree,
        records,
        model_name="WAG",
    )

    assert surface.model_name == "empirical protein matrix (WAG)"
    assert surface.root_prior_source == "stationary"
    assert (
        surface.normalization_policy == EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY
    )
    assert analysis.context.surface.model_name == "empirical protein matrix (WAG)"
    assert analysis.context.surface.root_prior_source == "stationary"
    assert (
        analysis.context.surface.normalization_policy
        == EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY
    )
    assert math.isclose(
        evaluate_compiled_finite_state_likelihood_analysis(analysis),
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_named_empirical_protein_compiled_site_rows_reconstruct_report() -> None:
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta")
    )
    analysis = compile_named_empirical_protein_likelihood_analysis(
        tree,
        records,
        model_name="LG",
    )
    report = evaluate_named_empirical_protein_tree_likelihood(
        tree,
        records,
        model_name="LG",
    )
    rows, weighted_total = evaluate_compiled_finite_state_site_log_likelihood_rows(
        analysis
    )

    assert len(rows) == report.site_count
    assert [row.source_site_position for row in rows] == [1, 2, 3, 4]
    assert all(row.partition_name == "full-alignment" for row in rows)
    assert all(row.partition_data_type is None for row in rows)
    assert math.isclose(
        weighted_total,
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
