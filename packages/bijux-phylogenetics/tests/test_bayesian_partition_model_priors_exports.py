from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    PARTITION_MODEL_PRIOR_TARGETS,
    PARTITION_PARAMETER_LINKAGE_POLICIES,
    PARTITION_SUBSTITUTION_BASE_MODELS,
    PartitionModelPriorBundle,
    PartitionModelPriorEvaluationReport,
    PartitionModelPriorRow,
    PartitionParameterLinkagePlan,
    PartitionSubstitutionModelDefinition,
    PartitionSubstitutionParameterState,
    build_partition_model_prior_bundle,
    build_partition_parameter_linkage_plan,
    build_partition_substitution_model_definition,
    evaluate_partition_model_log_prior,
    validate_partition_substitution_model_name,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    PARTITION_MODEL_PRIOR_TARGETS as PARTITION_MODEL_PRIOR_TARGETS_IMPL,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    PARTITION_PARAMETER_LINKAGE_POLICIES as PARTITION_PARAMETER_LINKAGE_POLICIES_IMPL,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    PARTITION_SUBSTITUTION_BASE_MODELS as PARTITION_SUBSTITUTION_BASE_MODELS_IMPL,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    PartitionModelPriorBundle as PartitionModelPriorBundleImpl,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    PartitionModelPriorEvaluationReport as PartitionModelPriorEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    PartitionModelPriorRow as PartitionModelPriorRowImpl,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    PartitionParameterLinkagePlan as PartitionParameterLinkagePlanImpl,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    PartitionSubstitutionModelDefinition as PartitionSubstitutionModelDefinitionImpl,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    PartitionSubstitutionParameterState as PartitionSubstitutionParameterStateImpl,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    build_partition_model_prior_bundle as build_partition_model_prior_bundle_impl,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    build_partition_parameter_linkage_plan as build_partition_parameter_linkage_plan_impl,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    build_partition_substitution_model_definition as build_partition_substitution_model_definition_impl,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    evaluate_partition_model_log_prior as evaluate_partition_model_log_prior_impl,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    validate_partition_substitution_model_name as validate_partition_substitution_model_name_impl,
)


def test_bayesian_exports_partition_model_prior_surface() -> None:
    assert PARTITION_PARAMETER_LINKAGE_POLICIES == (
        PARTITION_PARAMETER_LINKAGE_POLICIES_IMPL
    )
    assert PARTITION_MODEL_PRIOR_TARGETS == PARTITION_MODEL_PRIOR_TARGETS_IMPL
    assert PARTITION_SUBSTITUTION_BASE_MODELS == PARTITION_SUBSTITUTION_BASE_MODELS_IMPL
    assert (
        PartitionSubstitutionModelDefinition is PartitionSubstitutionModelDefinitionImpl
    )
    assert (
        PartitionSubstitutionParameterState is PartitionSubstitutionParameterStateImpl
    )
    assert PartitionParameterLinkagePlan is PartitionParameterLinkagePlanImpl
    assert PartitionModelPriorBundle is PartitionModelPriorBundleImpl
    assert (
        PartitionModelPriorEvaluationReport is PartitionModelPriorEvaluationReportImpl
    )
    assert PartitionModelPriorRow is PartitionModelPriorRowImpl
    assert (
        build_partition_substitution_model_definition
        is build_partition_substitution_model_definition_impl
    )
    assert (
        build_partition_parameter_linkage_plan
        is build_partition_parameter_linkage_plan_impl
    )
    assert build_partition_model_prior_bundle is build_partition_model_prior_bundle_impl
    assert evaluate_partition_model_log_prior is evaluate_partition_model_log_prior_impl
    assert (
        validate_partition_substitution_model_name
        is validate_partition_substitution_model_name_impl
    )
