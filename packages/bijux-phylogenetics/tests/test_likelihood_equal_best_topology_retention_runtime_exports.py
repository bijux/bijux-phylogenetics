from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodEqualBestTreeReport,
    NucleotideLikelihoodEqualBestTreeRow,
    build_nucleotide_likelihood_equal_best_tree_report,
    record_nucleotide_likelihood_equal_best_topology,
    validate_nucleotide_likelihood_equal_best_likelihood_tolerance,
    validate_nucleotide_likelihood_equal_best_tree_cap,
)


def test_public_runtime_exports_equal_best_topology_retention_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodEqualBestTreeReport
        is NucleotideLikelihoodEqualBestTreeReport
    )
    assert (
        likelihood_api.NucleotideLikelihoodEqualBestTreeRow
        is NucleotideLikelihoodEqualBestTreeRow
    )
    assert (
        likelihood_api.build_nucleotide_likelihood_equal_best_tree_report
        is build_nucleotide_likelihood_equal_best_tree_report
    )
    assert (
        likelihood_api.record_nucleotide_likelihood_equal_best_topology
        is record_nucleotide_likelihood_equal_best_topology
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_equal_best_likelihood_tolerance
        is validate_nucleotide_likelihood_equal_best_likelihood_tolerance
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_equal_best_tree_cap
        is validate_nucleotide_likelihood_equal_best_tree_cap
    )
