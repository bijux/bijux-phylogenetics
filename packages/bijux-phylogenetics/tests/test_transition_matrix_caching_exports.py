from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    FiniteStateTransitionMatrixEvaluator,
    build_transition_matrix_evaluator,
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
