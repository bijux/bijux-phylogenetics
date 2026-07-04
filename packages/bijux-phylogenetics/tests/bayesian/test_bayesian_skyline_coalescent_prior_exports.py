from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    COALESCENT_TREE_PRIOR_FAMILIES,
    SkylineCoalescentEpoch,
    SkylineCoalescentPriorEvaluationReport,
    SkylineCoalescentPriorModel,
    SkylineCoalescentSegmentRow,
    build_skyline_coalescent_tree_prior,
    evaluate_skyline_coalescent_tree_log_prior,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    COALESCENT_TREE_PRIOR_FAMILIES as COALESCENT_TREE_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    SkylineCoalescentEpoch as SkylineCoalescentEpochImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    SkylineCoalescentPriorEvaluationReport as SkylineCoalescentPriorEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    SkylineCoalescentPriorModel as SkylineCoalescentPriorModelImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    SkylineCoalescentSegmentRow as SkylineCoalescentSegmentRowImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    build_skyline_coalescent_tree_prior as build_skyline_coalescent_tree_prior_impl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    evaluate_skyline_coalescent_tree_log_prior as evaluate_skyline_coalescent_tree_log_prior_impl,
)


def test_bayesian_exports_skyline_coalescent_prior_surface() -> None:
    assert COALESCENT_TREE_PRIOR_FAMILIES == COALESCENT_TREE_PRIOR_FAMILIES_IMPL
    assert SkylineCoalescentEpoch is SkylineCoalescentEpochImpl
    assert SkylineCoalescentPriorModel is SkylineCoalescentPriorModelImpl
    assert SkylineCoalescentSegmentRow is SkylineCoalescentSegmentRowImpl
    assert (
        SkylineCoalescentPriorEvaluationReport
        is SkylineCoalescentPriorEvaluationReportImpl
    )
    assert (
        build_skyline_coalescent_tree_prior is build_skyline_coalescent_tree_prior_impl
    )
    assert (
        evaluate_skyline_coalescent_tree_log_prior
        is evaluate_skyline_coalescent_tree_log_prior_impl
    )
