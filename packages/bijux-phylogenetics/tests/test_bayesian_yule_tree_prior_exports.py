from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    YULE_TREE_PRIOR_FAMILIES,
    YuleTreePriorEvaluationReport,
    YuleTreePriorIntervalRow,
    YuleTreePriorModel,
    build_crown_conditioned_yule_tree_prior,
    evaluate_yule_tree_log_prior,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    YULE_TREE_PRIOR_FAMILIES as YULE_TREE_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    YuleTreePriorEvaluationReport as YuleTreePriorEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    YuleTreePriorIntervalRow as YuleTreePriorIntervalRowImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    YuleTreePriorModel as YuleTreePriorModelImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    build_crown_conditioned_yule_tree_prior as build_crown_conditioned_yule_tree_prior_impl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    evaluate_yule_tree_log_prior as evaluate_yule_tree_log_prior_impl,
)


def test_bayesian_exports_yule_tree_prior_surface() -> None:
    assert YULE_TREE_PRIOR_FAMILIES == YULE_TREE_PRIOR_FAMILIES_IMPL
    assert YuleTreePriorModel is YuleTreePriorModelImpl
    assert YuleTreePriorIntervalRow is YuleTreePriorIntervalRowImpl
    assert YuleTreePriorEvaluationReport is YuleTreePriorEvaluationReportImpl
    assert (
        build_crown_conditioned_yule_tree_prior
        is build_crown_conditioned_yule_tree_prior_impl
    )
    assert evaluate_yule_tree_log_prior is evaluate_yule_tree_log_prior_impl
