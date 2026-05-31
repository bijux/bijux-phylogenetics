from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    BIRTH_DEATH_TREE_PRIOR_FAMILIES,
    TIME_TREE_PRIOR_CONDITIONING_MODES,
    BirthDeathTreePriorBranchingRow,
    BirthDeathTreePriorEvaluationReport,
    BirthDeathTreePriorModel,
    build_crown_conditioned_birth_death_tree_prior,
    evaluate_birth_death_tree_log_prior,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    BIRTH_DEATH_TREE_PRIOR_FAMILIES as BIRTH_DEATH_TREE_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    TIME_TREE_PRIOR_CONDITIONING_MODES as TIME_TREE_PRIOR_CONDITIONING_MODES_IMPL,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    BirthDeathTreePriorBranchingRow as BirthDeathTreePriorBranchingRowImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    BirthDeathTreePriorEvaluationReport as BirthDeathTreePriorEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    BirthDeathTreePriorModel as BirthDeathTreePriorModelImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    build_crown_conditioned_birth_death_tree_prior as build_crown_conditioned_birth_death_tree_prior_impl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    evaluate_birth_death_tree_log_prior as evaluate_birth_death_tree_log_prior_impl,
)


def test_bayesian_exports_birth_death_tree_prior_surface() -> None:
    assert BIRTH_DEATH_TREE_PRIOR_FAMILIES == BIRTH_DEATH_TREE_PRIOR_FAMILIES_IMPL
    assert TIME_TREE_PRIOR_CONDITIONING_MODES == TIME_TREE_PRIOR_CONDITIONING_MODES_IMPL
    assert BirthDeathTreePriorModel is BirthDeathTreePriorModelImpl
    assert BirthDeathTreePriorBranchingRow is BirthDeathTreePriorBranchingRowImpl
    assert (
        BirthDeathTreePriorEvaluationReport is BirthDeathTreePriorEvaluationReportImpl
    )
    assert (
        build_crown_conditioned_birth_death_tree_prior
        is build_crown_conditioned_birth_death_tree_prior_impl
    )
    assert (
        evaluate_birth_death_tree_log_prior is evaluate_birth_death_tree_log_prior_impl
    )
