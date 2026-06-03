from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
import re

from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    SUBSTITUTION_PARAMETER_PRIOR_TARGETS,
    SubstitutionParameterPriorBundle,
    build_substitution_parameter_prior_bundle,
    evaluate_substitution_parameter_log_prior,
)
from bijux_phylogenetics.phylo.alignment.partitions import normalize_partition_data_type
from bijux_phylogenetics.phylo.likelihood.dna import (
    normalize_dna_exchangeabilities_by_anchor,
    validate_dna_base_frequencies,
    validate_positive_kappa,
)
from bijux_phylogenetics.phylo.likelihood.gamma import validate_discrete_gamma_alpha
from bijux_phylogenetics.phylo.likelihood.invariant import (
    validate_invariant_proportion,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

PARTITION_PARAMETER_LINKAGE_POLICIES = ("linked", "unlinked")
PARTITION_MODEL_PRIOR_TARGETS = SUBSTITUTION_PARAMETER_PRIOR_TARGETS
PARTITION_SUBSTITUTION_BASE_MODELS = ("F81", "GTR", "HKY85", "JC69", "K80")

_PARTITION_MODEL_MODIFIER = re.compile(r"^G\d+$")
_SCALAR_VALUE_TOLERANCE = 1e-12


@dataclass(frozen=True, slots=True)
class PartitionSubstitutionModelDefinition:
    """One validated substitution-model definition for one named partition."""

    partition_name: str
    model_name: str
    base_model_name: str
    gamma_enabled: bool
    invariant_enabled: bool
    empirical_base_frequencies_enabled: bool
    data_type: str | None = None

    def required_targets(self) -> tuple[str, ...]:
        """Return the substitution-parameter targets required by this partition."""
        targets: list[str] = []
        if self.base_model_name == "K80":
            targets.append("kappa")
        elif self.base_model_name == "F81":
            targets.append("base-frequencies")
        elif self.base_model_name == "HKY85":
            targets.extend(("kappa", "base-frequencies"))
        elif self.base_model_name == "GTR":
            targets.extend(("exchangeabilities", "base-frequencies"))
        if self.gamma_enabled:
            targets.append("gamma-alpha")
        if self.invariant_enabled:
            targets.append("invariant-proportion")
        return tuple(targets)


@dataclass(frozen=True, slots=True)
class PartitionSubstitutionParameterState:
    """One explicit realized substitution-parameter state for one partition."""

    partition_name: str
    kappa: float | None = None
    exchangeabilities: (
        Mapping[tuple[str, str], float] | Mapping[str, float] | Sequence[float] | None
    ) = None
    base_frequencies: Mapping[str, float] | Sequence[float] | None = None
    gamma_alpha: float | None = None
    invariant_proportion: float | None = None


@dataclass(frozen=True, slots=True)
class PartitionParameterLinkagePlan:
    """One explicit linkage plan for partition substitution parameters."""

    partition_names: tuple[str, ...]
    target_partition_groups: dict[str, dict[str, str]]

    def groups_for_target(self, target_name: str) -> dict[str, str]:
        return dict(self.target_partition_groups[target_name])


@dataclass(frozen=True, slots=True)
class PartitionModelPriorBundle:
    """One explicit prior bundle over partition substitution models."""

    partition_models: tuple[PartitionSubstitutionModelDefinition, ...]
    linkage_plan: PartitionParameterLinkagePlan
    substitution_prior_bundle: SubstitutionParameterPriorBundle


@dataclass(frozen=True, slots=True)
class PartitionModelPriorRow:
    """One grouped partition-model prior contribution."""

    target_name: str
    group_name: str
    partition_names: tuple[str, ...]
    partition_model_names: dict[str, str]
    family: str
    component_values: dict[str, float]
    hyperparameter_values: dict[str, float]
    log_prior_contribution: float


@dataclass(frozen=True, slots=True)
class PartitionModelPriorEvaluationReport:
    """One partition-model prior evaluation report."""

    partition_count: int
    parameter_count: int
    total_log_prior: float
    rows: list[PartitionModelPriorRow]


def validate_partition_substitution_model_name(model_name: str) -> str:
    """Validate and canonicalize one partition substitution-model label."""
    tokens = [token for token in model_name.strip().upper().split("+") if token]
    if not tokens:
        raise PhylogeneticsError(
            "partition substitution model name must not be empty",
            code="partition_model_name_empty",
        )
    base_model_name = tokens[0]
    if base_model_name not in PARTITION_SUBSTITUTION_BASE_MODELS:
        raise PhylogeneticsError(
            "partition substitution model must use one owned nucleotide base model",
            code="partition_model_name_invalid",
            details={
                "model_name": model_name,
                "allowed_base_models": list(PARTITION_SUBSTITUTION_BASE_MODELS),
            },
        )
    empirical_base_frequencies_enabled = False
    gamma_enabled = False
    invariant_enabled = False
    for modifier in tokens[1:]:
        if modifier == "F":
            if base_model_name not in {"F81", "HKY85", "GTR"}:
                raise PhylogeneticsError(
                    "partition substitution model uses +F on a base model without owned frequency parameters",
                    code="partition_model_frequency_modifier_invalid",
                    details={"model_name": model_name},
                )
            empirical_base_frequencies_enabled = True
            continue
        if modifier == "I":
            invariant_enabled = True
            continue
        if modifier == "G" or _PARTITION_MODEL_MODIFIER.fullmatch(modifier) is not None:
            gamma_enabled = True
            continue
        raise PhylogeneticsError(
            "partition substitution model uses an unsupported modifier",
            code="partition_model_modifier_invalid",
            details={"model_name": model_name, "modifier": modifier},
        )
    canonical_tokens = [base_model_name]
    if empirical_base_frequencies_enabled:
        canonical_tokens.append("F")
    if gamma_enabled:
        canonical_tokens.append("G")
    if invariant_enabled:
        canonical_tokens.append("I")
    return "+".join(canonical_tokens)


def build_partition_substitution_model_definition(
    *,
    partition_name: str,
    model_name: str,
    data_type: str | None = "DNA",
) -> PartitionSubstitutionModelDefinition:
    """Build one validated per-partition substitution-model definition."""
    validated_partition_name = _validate_partition_name(partition_name)
    canonical_model_name = validate_partition_substitution_model_name(model_name)
    normalized_data_type = normalize_partition_data_type(data_type)
    if normalized_data_type not in {None, "DNA"}:
        raise PhylogeneticsError(
            "partition model priors currently support DNA partitions only",
            code="partition_model_data_type_invalid",
            details={
                "partition_name": validated_partition_name,
                "data_type": normalized_data_type,
            },
        )
    model_tokens = canonical_model_name.split("+")
    base_model_name = model_tokens[0]
    modifiers = set(model_tokens[1:])
    return PartitionSubstitutionModelDefinition(
        partition_name=validated_partition_name,
        model_name=canonical_model_name,
        base_model_name=base_model_name,
        gamma_enabled="G" in modifiers,
        invariant_enabled="I" in modifiers,
        empirical_base_frequencies_enabled="F" in modifiers,
        data_type=normalized_data_type,
    )


def build_partition_parameter_linkage_plan(
    *,
    partition_names: Sequence[str],
    linkage_policies: Mapping[str, str] | None = None,
    group_assignments: Mapping[str, Mapping[str, str]] | None = None,
) -> PartitionParameterLinkagePlan:
    """Build one explicit per-target partition linkage plan."""
    validated_partition_names = _validate_unique_partition_names(partition_names)
    normalized_linkage_policies = _normalize_linkage_policies(linkage_policies)
    normalized_group_assignments = _normalize_group_assignments(
        group_assignments,
        partition_names=validated_partition_names,
    )
    target_partition_groups: dict[str, dict[str, str]] = {}
    for target_name in PARTITION_MODEL_PRIOR_TARGETS:
        explicit_groups = normalized_group_assignments.get(target_name)
        if explicit_groups is not None:
            target_partition_groups[target_name] = explicit_groups
            continue
        linkage_policy = normalized_linkage_policies.get(target_name, "linked")
        if linkage_policy == "linked":
            target_partition_groups[target_name] = dict.fromkeys(
                validated_partition_names,
                f"{target_name}-shared",
            )
            continue
        target_partition_groups[target_name] = {
            partition_name: partition_name
            for partition_name in validated_partition_names
        }
    return PartitionParameterLinkagePlan(
        partition_names=validated_partition_names,
        target_partition_groups=target_partition_groups,
    )


def build_partition_model_prior_bundle(
    *,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
    linkage_plan: PartitionParameterLinkagePlan,
    substitution_prior_bundle: SubstitutionParameterPriorBundle,
) -> PartitionModelPriorBundle:
    """Build one explicit prior bundle over partition substitution models."""
    validated_partition_models = tuple(partition_models)
    if not validated_partition_models:
        raise PhylogeneticsError(
            "partition model prior bundle requires at least one partition model",
            code="partition_model_prior_bundle_empty",
        )
    validated_partition_names = _validate_unique_partition_names(
        model.partition_name for model in validated_partition_models
    )
    if validated_partition_names != linkage_plan.partition_names:
        raise PhylogeneticsError(
            "partition model prior bundle requires linkage coverage for the exact partition set",
            code="partition_model_linkage_partition_mismatch",
            details={
                "model_partition_names": list(validated_partition_names),
                "linkage_partition_names": list(linkage_plan.partition_names),
            },
        )
    required_targets = {
        target_name
        for model in validated_partition_models
        for target_name in model.required_targets()
    }
    _validate_prior_target_usage(
        substitution_prior_bundle=substitution_prior_bundle,
        required_targets=required_targets,
    )
    return PartitionModelPriorBundle(
        partition_models=validated_partition_models,
        linkage_plan=linkage_plan,
        substitution_prior_bundle=substitution_prior_bundle,
    )


def evaluate_partition_model_log_prior(
    *,
    prior_bundle: PartitionModelPriorBundle,
    partition_parameter_states: Sequence[PartitionSubstitutionParameterState],
) -> PartitionModelPriorEvaluationReport:
    """Evaluate one partition-model prior bundle on explicit partition parameters."""
    state_by_partition_name = _build_partition_state_lookup(
        partition_parameter_states,
        expected_partition_names=tuple(
            model.partition_name for model in prior_bundle.partition_models
        ),
    )
    rows: list[PartitionModelPriorRow] = []
    parameter_count = 0
    partition_model_by_name = {
        model.partition_name: model for model in prior_bundle.partition_models
    }
    for target_name in PARTITION_MODEL_PRIOR_TARGETS:
        target_prior_bundle = _target_prior_bundle(
            substitution_prior_bundle=prior_bundle.substitution_prior_bundle,
            target_name=target_name,
        )
        if target_prior_bundle is None:
            continue
        required_partition_names = tuple(
            model.partition_name
            for model in prior_bundle.partition_models
            if target_name in model.required_targets()
        )
        if not required_partition_names:
            continue
        grouped_partition_names = _group_partition_names_for_target(
            prior_bundle.linkage_plan.groups_for_target(target_name),
            required_partition_names=required_partition_names,
        )
        for group_name, grouped_names in grouped_partition_names:
            representative_state = state_by_partition_name[grouped_names[0]]
            representative_value = _extract_partition_target_value(
                target_name=target_name,
                state=representative_state,
            )
            _validate_linked_target_values(
                target_name=target_name,
                representative_partition_name=grouped_names[0],
                representative_value=representative_value,
                grouped_partition_names=grouped_names,
                state_by_partition_name=state_by_partition_name,
            )
            target_report = _evaluate_target_group(
                target_name=target_name,
                target_prior_bundle=target_prior_bundle,
                realized_value=representative_value,
            )
            target_row = target_report.rows[0]
            rows.append(
                PartitionModelPriorRow(
                    target_name=target_name,
                    group_name=group_name,
                    partition_names=grouped_names,
                    partition_model_names={
                        partition_name: partition_model_by_name[
                            partition_name
                        ].model_name
                        for partition_name in grouped_names
                    },
                    family=target_row.family,
                    component_values=target_row.component_values,
                    hyperparameter_values=target_row.hyperparameter_values,
                    log_prior_contribution=target_row.log_prior_contribution,
                )
            )
            parameter_count += target_report.prior_count
    total_log_prior = math.fsum(row.log_prior_contribution for row in rows)
    return PartitionModelPriorEvaluationReport(
        partition_count=len(prior_bundle.partition_models),
        parameter_count=parameter_count,
        total_log_prior=float(format(total_log_prior, ".15g")),
        rows=rows,
    )


def _validate_partition_name(partition_name: str) -> str:
    normalized_partition_name = partition_name.strip()
    if not normalized_partition_name:
        raise PhylogeneticsError(
            "partition name must not be empty",
            code="partition_model_partition_name_empty",
        )
    return normalized_partition_name


def _validate_unique_partition_names(
    partition_names: Sequence[str] | Sequence[object],
) -> tuple[str, ...]:
    normalized_partition_names = tuple(
        _validate_partition_name(str(partition_name))
        for partition_name in partition_names
    )
    duplicate_partition_names = sorted(
        {
            partition_name
            for partition_name in normalized_partition_names
            if normalized_partition_names.count(partition_name) > 1
        }
    )
    if duplicate_partition_names:
        raise PhylogeneticsError(
            "partition names must be unique in one partition-model prior surface",
            code="partition_model_partition_names_duplicated",
            details={"duplicate_partition_names": duplicate_partition_names},
        )
    return normalized_partition_names


def _normalize_linkage_policies(
    linkage_policies: Mapping[str, str] | None,
) -> dict[str, str]:
    if linkage_policies is None:
        return {}
    normalized: dict[str, str] = {}
    for target_name, linkage_policy in linkage_policies.items():
        _validate_partition_target_name(target_name)
        normalized_policy = linkage_policy.strip().lower()
        if normalized_policy not in PARTITION_PARAMETER_LINKAGE_POLICIES:
            raise PhylogeneticsError(
                "partition parameter linkage policy is unsupported",
                code="partition_model_linkage_policy_invalid",
                details={
                    "target_name": target_name,
                    "linkage_policy": linkage_policy,
                    "allowed_policies": list(PARTITION_PARAMETER_LINKAGE_POLICIES),
                },
            )
        normalized[target_name] = normalized_policy
    return normalized


def _normalize_group_assignments(
    group_assignments: Mapping[str, Mapping[str, str]] | None,
    *,
    partition_names: tuple[str, ...],
) -> dict[str, dict[str, str]]:
    if group_assignments is None:
        return {}
    normalized_assignments: dict[str, dict[str, str]] = {}
    for target_name, partition_groups in group_assignments.items():
        _validate_partition_target_name(target_name)
        group_by_partition_name = {
            _validate_partition_name(partition_name): group_name.strip()
            for partition_name, group_name in partition_groups.items()
        }
        if set(group_by_partition_name) != set(partition_names):
            raise PhylogeneticsError(
                "partition group assignments must cover the exact partition set",
                code="partition_model_group_assignments_incomplete",
                details={
                    "target_name": target_name,
                    "expected_partition_names": list(partition_names),
                    "observed_partition_names": sorted(group_by_partition_name),
                },
            )
        empty_group_names = sorted(
            partition_name
            for partition_name, group_name in group_by_partition_name.items()
            if not group_name
        )
        if empty_group_names:
            raise PhylogeneticsError(
                "partition group assignments must use non-empty group names",
                code="partition_model_group_name_empty",
                details={
                    "target_name": target_name,
                    "partition_names": empty_group_names,
                },
            )
        normalized_assignments[target_name] = group_by_partition_name
    return normalized_assignments


def _validate_partition_target_name(target_name: str) -> None:
    if target_name not in PARTITION_MODEL_PRIOR_TARGETS:
        raise PhylogeneticsError(
            "partition model prior target is unsupported",
            code="partition_model_target_invalid",
            details={
                "target_name": target_name,
                "allowed_targets": list(PARTITION_MODEL_PRIOR_TARGETS),
            },
        )


def _validate_prior_target_usage(
    *,
    substitution_prior_bundle: SubstitutionParameterPriorBundle,
    required_targets: set[str],
) -> None:
    configured_target_names: set[str] = set()
    if substitution_prior_bundle.kappa_prior is not None:
        configured_target_names.add("kappa")
    if substitution_prior_bundle.exchangeability_prior is not None:
        configured_target_names.add("exchangeabilities")
    if substitution_prior_bundle.base_frequency_prior is not None:
        configured_target_names.add("base-frequencies")
    if substitution_prior_bundle.gamma_alpha_prior is not None:
        configured_target_names.add("gamma-alpha")
    if substitution_prior_bundle.invariant_proportion_prior is not None:
        configured_target_names.add("invariant-proportion")
    unused_target_names = sorted(configured_target_names - required_targets)
    if unused_target_names:
        raise PhylogeneticsError(
            "partition model prior bundle configures substitution priors that no partition model can use",
            code="partition_model_prior_targets_unused",
            details={"unused_target_names": unused_target_names},
        )


def _build_partition_state_lookup(
    partition_parameter_states: Sequence[PartitionSubstitutionParameterState],
    *,
    expected_partition_names: tuple[str, ...],
) -> dict[str, PartitionSubstitutionParameterState]:
    if not partition_parameter_states:
        raise PhylogeneticsError(
            "partition model prior evaluation requires explicit partition parameter states",
            code="partition_model_parameter_states_empty",
        )
    state_by_partition_name: dict[str, PartitionSubstitutionParameterState] = {}
    for state in partition_parameter_states:
        partition_name = _validate_partition_name(state.partition_name)
        if partition_name in state_by_partition_name:
            raise PhylogeneticsError(
                "partition model prior evaluation received more than one state for the same partition",
                code="partition_model_parameter_state_duplicated",
                details={"partition_name": partition_name},
            )
        state_by_partition_name[partition_name] = state
    if set(state_by_partition_name) != set(expected_partition_names):
        raise PhylogeneticsError(
            "partition model prior evaluation requires states for the exact partition set",
            code="partition_model_parameter_state_partition_mismatch",
            details={
                "expected_partition_names": list(expected_partition_names),
                "observed_partition_names": sorted(state_by_partition_name),
            },
        )
    return state_by_partition_name


def _group_partition_names_for_target(
    group_by_partition_name: Mapping[str, str],
    *,
    required_partition_names: tuple[str, ...],
) -> list[tuple[str, tuple[str, ...]]]:
    grouped_partition_names: dict[str, list[str]] = defaultdict(list)
    for partition_name in required_partition_names:
        group_name = group_by_partition_name[partition_name]
        grouped_partition_names[group_name].append(partition_name)
    return [
        (group_name, tuple(partition_names))
        for group_name, partition_names in sorted(grouped_partition_names.items())
    ]


def _extract_partition_target_value(
    *,
    target_name: str,
    state: PartitionSubstitutionParameterState,
) -> object:
    if target_name == "kappa":
        if state.kappa is None:
            raise PhylogeneticsError(
                "partition model prior evaluation requires kappa for one partition group",
                code="partition_model_missing_kappa",
                details={"partition_name": state.partition_name},
            )
        return state.kappa
    if target_name == "exchangeabilities":
        if state.exchangeabilities is None:
            raise PhylogeneticsError(
                "partition model prior evaluation requires exchangeabilities for one partition group",
                code="partition_model_missing_exchangeabilities",
                details={"partition_name": state.partition_name},
            )
        return state.exchangeabilities
    if target_name == "base-frequencies":
        if state.base_frequencies is None:
            raise PhylogeneticsError(
                "partition model prior evaluation requires base frequencies for one partition group",
                code="partition_model_missing_base_frequencies",
                details={"partition_name": state.partition_name},
            )
        return state.base_frequencies
    if target_name == "gamma-alpha":
        if state.gamma_alpha is None:
            raise PhylogeneticsError(
                "partition model prior evaluation requires gamma alpha for one partition group",
                code="partition_model_missing_gamma_alpha",
                details={"partition_name": state.partition_name},
            )
        return state.gamma_alpha
    if target_name == "invariant-proportion":
        if state.invariant_proportion is None:
            raise PhylogeneticsError(
                "partition model prior evaluation requires invariant proportion for one partition group",
                code="partition_model_missing_invariant_proportion",
                details={"partition_name": state.partition_name},
            )
        return state.invariant_proportion
    raise AssertionError(f"unsupported target {target_name!r}")


def _validate_linked_target_values(
    *,
    target_name: str,
    representative_partition_name: str,
    representative_value: object,
    grouped_partition_names: tuple[str, ...],
    state_by_partition_name: Mapping[str, PartitionSubstitutionParameterState],
) -> None:
    representative_signature = _target_value_signature(
        target_name=target_name,
        value=representative_value,
    )
    mismatched_partition_names: list[str] = []
    for partition_name in grouped_partition_names[1:]:
        partition_value = _extract_partition_target_value(
            target_name=target_name,
            state=state_by_partition_name[partition_name],
        )
        if (
            _target_value_signature(target_name=target_name, value=partition_value)
            != representative_signature
        ):
            mismatched_partition_names.append(partition_name)
    if mismatched_partition_names:
        raise PhylogeneticsError(
            "linked partition parameter group received mismatched realized values",
            code="partition_model_linked_parameter_values_mismatched",
            details={
                "target_name": target_name,
                "representative_partition_name": representative_partition_name,
                "mismatched_partition_names": mismatched_partition_names,
            },
        )


def _target_value_signature(*, target_name: str, value: object) -> object:
    if target_name == "kappa":
        if not isinstance(value, float):
            raise PhylogeneticsError(
                "partition prior target 'kappa' requires a float realized value",
                code="partition_model_prior_target_value_invalid",
                details={
                    "target_name": target_name,
                    "value_type": type(value).__name__,
                },
            )
        return float(
            format(validate_positive_kappa(value, model_name="partition prior"), ".15g")
        )
    if target_name == "exchangeabilities":
        normalized = normalize_dna_exchangeabilities_by_anchor(
            value,
            model_name="partition prior",
        )
        return tuple(
            float(format(component_value, ".15g")) for component_value in normalized
        )
    if target_name == "base-frequencies":
        normalized = validate_dna_base_frequencies(
            value,
            model_name="partition prior",
        )
        return tuple(
            float(format(component_value, ".15g")) for component_value in normalized
        )
    if target_name == "gamma-alpha":
        if not isinstance(value, float):
            raise PhylogeneticsError(
                "partition prior target 'gamma-alpha' requires a float realized value",
                code="partition_model_prior_target_value_invalid",
                details={
                    "target_name": target_name,
                    "value_type": type(value).__name__,
                },
            )
        return float(
            format(
                validate_discrete_gamma_alpha(value),
                ".15g",
            )
        )
    if target_name == "invariant-proportion":
        if not isinstance(value, float):
            raise PhylogeneticsError(
                "partition prior target 'invariant-proportion' requires a float realized value",
                code="partition_model_prior_target_value_invalid",
                details={
                    "target_name": target_name,
                    "value_type": type(value).__name__,
                },
            )
        return float(
            format(
                validate_invariant_proportion(
                    value,
                    model_name="partition prior",
                ),
                ".15g",
            )
        )
    raise AssertionError(f"unsupported target {target_name!r}")


def _evaluate_target_group(
    *,
    target_name: str,
    target_prior_bundle: SubstitutionParameterPriorBundle,
    realized_value: object,
) -> object:
    if target_name == "kappa":
        if not isinstance(realized_value, float):
            raise PhylogeneticsError(
                "partition prior evaluation target 'kappa' requires a float realized value",
                code="partition_model_prior_target_value_invalid",
                details={
                    "target_name": target_name,
                    "value_type": type(realized_value).__name__,
                },
            )
        return evaluate_substitution_parameter_log_prior(
            prior_bundle=target_prior_bundle,
            kappa=realized_value,
        )
    if target_name == "exchangeabilities":
        return evaluate_substitution_parameter_log_prior(
            prior_bundle=target_prior_bundle,
            exchangeabilities=realized_value,
        )
    if target_name == "base-frequencies":
        return evaluate_substitution_parameter_log_prior(
            prior_bundle=target_prior_bundle,
            base_frequencies=realized_value,
        )
    if target_name == "gamma-alpha":
        if not isinstance(realized_value, float):
            raise PhylogeneticsError(
                "partition prior evaluation target 'gamma-alpha' requires a float realized value",
                code="partition_model_prior_target_value_invalid",
                details={
                    "target_name": target_name,
                    "value_type": type(realized_value).__name__,
                },
            )
        return evaluate_substitution_parameter_log_prior(
            prior_bundle=target_prior_bundle,
            gamma_alpha=realized_value,
        )
    if target_name == "invariant-proportion":
        if not isinstance(realized_value, float):
            raise PhylogeneticsError(
                "partition prior evaluation target 'invariant-proportion' requires a float realized value",
                code="partition_model_prior_target_value_invalid",
                details={
                    "target_name": target_name,
                    "value_type": type(realized_value).__name__,
                },
            )
        return evaluate_substitution_parameter_log_prior(
            prior_bundle=target_prior_bundle,
            invariant_proportion=realized_value,
        )
    raise AssertionError(f"unsupported target {target_name!r}")


def _target_prior_bundle(
    *,
    substitution_prior_bundle: SubstitutionParameterPriorBundle,
    target_name: str,
) -> SubstitutionParameterPriorBundle | None:
    if target_name == "kappa":
        if substitution_prior_bundle.kappa_prior is None:
            return None
        return build_substitution_parameter_prior_bundle(
            kappa_prior=substitution_prior_bundle.kappa_prior
        )
    if target_name == "exchangeabilities":
        if substitution_prior_bundle.exchangeability_prior is None:
            return None
        return build_substitution_parameter_prior_bundle(
            exchangeability_prior=substitution_prior_bundle.exchangeability_prior
        )
    if target_name == "base-frequencies":
        if substitution_prior_bundle.base_frequency_prior is None:
            return None
        return build_substitution_parameter_prior_bundle(
            base_frequency_prior=substitution_prior_bundle.base_frequency_prior
        )
    if target_name == "gamma-alpha":
        if substitution_prior_bundle.gamma_alpha_prior is None:
            return None
        return build_substitution_parameter_prior_bundle(
            gamma_alpha_prior=substitution_prior_bundle.gamma_alpha_prior
        )
    if target_name == "invariant-proportion":
        if substitution_prior_bundle.invariant_proportion_prior is None:
            return None
        return build_substitution_parameter_prior_bundle(
            invariant_proportion_prior=substitution_prior_bundle.invariant_proportion_prior
        )
    raise AssertionError(f"unsupported target {target_name!r}")
