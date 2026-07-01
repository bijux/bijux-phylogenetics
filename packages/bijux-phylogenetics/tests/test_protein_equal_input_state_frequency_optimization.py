from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    optimize_protein_equal_input_state_frequencies,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_protein_equal_input_state_frequency_optimization_improves_fixed_tree_likelihood() -> (
    None
):
    report = optimize_protein_equal_input_state_frequencies(
        load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk")),
        load_fasta_alignment(
            fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta")
        ),
        max_coordinate_passes=4,
    )

    optimized_values = [row.optimized_value for row in report.state_rows]

    assert report.model_name == "protein equal-input"
    assert report.state_count == 20
    assert report.state_frequency_source == "optimized"
    assert report.fit_assessment.classification in {"failed", "overparameterized"}
    assert report.function_evaluation_count > 1
    assert report.optimization_pass_count >= 1
    assert report.optimized_log_likelihood >= report.initial_log_likelihood
    assert math.isclose(
        sum(optimized_values),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert min(optimized_values) > 0.0
    assert max(optimized_values) > min(optimized_values)
    if not report.converged:
        assert (
            "protein equal-input stationary-frequency search did not converge"
            in report.warnings
        )
    assert report.parameter_diagnostics
