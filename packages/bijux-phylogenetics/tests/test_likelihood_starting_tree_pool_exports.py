from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodStartingTreePoolReport,
    NucleotideLikelihoodStartingTreeSummary,
    build_nucleotide_likelihood_starting_tree_pool,
    build_nucleotide_likelihood_starting_tree_pool_from_alignment,
    validate_nucleotide_likelihood_random_start_tree_count,
    validate_nucleotide_likelihood_starting_tree_pool_model,
    write_nucleotide_likelihood_starting_tree_pool_artifacts,
    write_nucleotide_likelihood_starting_tree_pool_run_json,
    write_nucleotide_likelihood_starting_tree_score_table,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodStartingTreePoolReport as NucleotideLikelihoodStartingTreePoolReportImpl,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodStartingTreeSummary as NucleotideLikelihoodStartingTreeSummaryImpl,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_pool import (
    build_nucleotide_likelihood_starting_tree_pool as build_nucleotide_likelihood_starting_tree_pool_impl,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_pool import (
    build_nucleotide_likelihood_starting_tree_pool_from_alignment as build_nucleotide_likelihood_starting_tree_pool_from_alignment_impl,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_pool import (
    validate_nucleotide_likelihood_random_start_tree_count as validate_nucleotide_likelihood_random_start_tree_count_impl,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_pool import (
    validate_nucleotide_likelihood_starting_tree_pool_model as validate_nucleotide_likelihood_starting_tree_pool_model_impl,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_pool import (
    write_nucleotide_likelihood_starting_tree_pool_artifacts as write_nucleotide_likelihood_starting_tree_pool_artifacts_impl,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_pool import (
    write_nucleotide_likelihood_starting_tree_pool_run_json as write_nucleotide_likelihood_starting_tree_pool_run_json_impl,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_pool import (
    write_nucleotide_likelihood_starting_tree_score_table as write_nucleotide_likelihood_starting_tree_score_table_impl,
)


def test_likelihood_exports_starting_tree_pool_surface() -> None:
    assert (
        NucleotideLikelihoodStartingTreePoolReport
        is NucleotideLikelihoodStartingTreePoolReportImpl
    )
    assert (
        NucleotideLikelihoodStartingTreeSummary
        is NucleotideLikelihoodStartingTreeSummaryImpl
    )
    assert (
        build_nucleotide_likelihood_starting_tree_pool
        is build_nucleotide_likelihood_starting_tree_pool_impl
    )
    assert (
        build_nucleotide_likelihood_starting_tree_pool_from_alignment
        is build_nucleotide_likelihood_starting_tree_pool_from_alignment_impl
    )
    assert (
        validate_nucleotide_likelihood_random_start_tree_count
        is validate_nucleotide_likelihood_random_start_tree_count_impl
    )
    assert (
        validate_nucleotide_likelihood_starting_tree_pool_model
        is validate_nucleotide_likelihood_starting_tree_pool_model_impl
    )
    assert (
        write_nucleotide_likelihood_starting_tree_score_table
        is write_nucleotide_likelihood_starting_tree_score_table_impl
    )
    assert (
        write_nucleotide_likelihood_starting_tree_pool_run_json
        is write_nucleotide_likelihood_starting_tree_pool_run_json_impl
    )
    assert (
        write_nucleotide_likelihood_starting_tree_pool_artifacts
        is write_nucleotide_likelihood_starting_tree_pool_artifacts_impl
    )
