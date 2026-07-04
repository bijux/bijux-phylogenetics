from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    FiniteStateTransitionMatrixEvaluator,
    build_transition_matrix_evaluator,
    clear_shared_finite_state_transition_cache,
    shared_finite_state_transition_cache_entry_count,
    update_compiled_finite_state_analysis_branch_length,
    update_compiled_finite_state_analysis_topology,
    update_compiled_finite_state_context_branch_length,
    update_compiled_finite_state_context_topology,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    clear_shared_finite_state_transition_cache as clear_shared_finite_state_transition_cache_impl,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    shared_finite_state_transition_cache_entry_count as shared_finite_state_transition_cache_entry_count_impl,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    update_compiled_finite_state_analysis_branch_length as update_compiled_finite_state_analysis_branch_length_impl,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    update_compiled_finite_state_analysis_topology as update_compiled_finite_state_analysis_topology_impl,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    update_compiled_finite_state_context_branch_length as update_compiled_finite_state_context_branch_length_impl,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    update_compiled_finite_state_context_topology as update_compiled_finite_state_context_topology_impl,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    FiniteStateTransitionMatrixEvaluator as FiniteStateTransitionMatrixEvaluatorImpl,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator as build_transition_matrix_evaluator_impl,
)


def test_phylo_likelihood_exports_transition_matrix_caching_surface() -> None:
    assert (
        FiniteStateTransitionMatrixEvaluator is FiniteStateTransitionMatrixEvaluatorImpl
    )
    assert build_transition_matrix_evaluator is build_transition_matrix_evaluator_impl
    assert (
        clear_shared_finite_state_transition_cache
        is clear_shared_finite_state_transition_cache_impl
    )
    assert (
        shared_finite_state_transition_cache_entry_count
        is shared_finite_state_transition_cache_entry_count_impl
    )
    assert (
        update_compiled_finite_state_context_branch_length
        is update_compiled_finite_state_context_branch_length_impl
    )
    assert (
        update_compiled_finite_state_analysis_branch_length
        is update_compiled_finite_state_analysis_branch_length_impl
    )
    assert (
        update_compiled_finite_state_context_topology
        is update_compiled_finite_state_context_topology_impl
    )
    assert (
        update_compiled_finite_state_analysis_topology
        is update_compiled_finite_state_analysis_topology_impl
    )
