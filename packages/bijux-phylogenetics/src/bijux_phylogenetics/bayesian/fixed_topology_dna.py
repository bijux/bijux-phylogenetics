from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math

from bijux_phylogenetics.bayesian.branch_length_priors import (
    BranchLengthPriorModel,
    evaluate_tree_branch_length_log_prior,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsRunReport,
    propose_base_frequency_simplex_move,
    propose_branch_length_scaling_move,
    propose_gtr_exchangeability_move,
    propose_kappa_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
    build_bayesian_prior_component_state,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    SubstitutionParameterPriorBundle,
    evaluate_substitution_parameter_log_prior,
)
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_k80_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.dna import validate_positive_kappa
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    estimate_empirical_dna_base_frequencies_from_records,
    normalize_dna_likelihood_records,
    validate_dna_observation_policy,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    DNA_EXCHANGEABILITY_LABELS,
    parameterize_dna_base_frequency_simplex,
    parameterize_dna_exchangeability_simplex,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXED_TOPOLOGY_DNA_SUBSTITUTION_MODELS = ("JC69", "K80", "HKY85", "GTR")
_SUPPORTED_FIXED_TOPOLOGY_DNA_MODELS = ("K80", "HKY85", "GTR")
_ACTIVE_PARAMETER_TARGETS_BY_MODEL = {
    "JC69": (),
    "K80": ("kappa",),
    "HKY85": ("base-frequencies", "kappa"),
    "GTR": ("base-frequencies", "exchangeabilities"),
}
_DEFAULT_INITIAL_KAPPA = 2.0
_DEFAULT_INITIAL_EXCHANGEABILITIES = dict.fromkeys(DNA_EXCHANGEABILITY_LABELS, 1.0)
_MINIMUM_SIMPLEX_COMPONENT = 1e-6


@dataclass(frozen=True, slots=True)
class FixedTopologyDnaModelDefinition:
    """One validated fixed-topology DNA posterior model definition."""

    substitution_model_name: str
    branch_length_prior: BranchLengthPriorModel
    substitution_parameter_prior_bundle: SubstitutionParameterPriorBundle
    initial_kappa: float | None = None
    initial_base_frequencies: dict[str, float] | None = None
    initial_exchangeabilities: dict[str, float] | None = None

    @property
    def active_parameter_targets(self) -> tuple[str, ...]:
        return _ACTIVE_PARAMETER_TARGETS_BY_MODEL[self.substitution_model_name]


@dataclass(frozen=True, slots=True)
class FixedTopologyDnaProposalSchedule:
    """One validated proposal schedule for fixed-topology DNA posterior sampling."""

    substitution_model_name: str
    branch_length_move_weight: float
    branch_length_log_scale_standard_deviation: float
    kappa_move_weight: float
    kappa_log_scale_standard_deviation: float | None
    base_frequency_move_weight: float
    base_frequency_coordinate_standard_deviation: float | None
    exchangeability_move_weight: float
    exchangeability_coordinate_standard_deviation: float | None


@dataclass(frozen=True, slots=True)
class FixedTopologyDnaPosteriorRow:
    """One sampled posterior row from a fixed-topology DNA chain."""

    sample_index: int
    iteration_index: int
    topology_id: str
    substitution_model_name: str
    total_log_prior: float
    log_likelihood: float
    posterior_log_score: float
    prior_component_log_priors: dict[str, float]
    branch_lengths: dict[str, float]
    scalar_parameters: dict[str, float]
    vector_parameters: dict[str, dict[str, float]]


@dataclass(frozen=True, slots=True)
class FixedTopologyDnaRunReport:
    """One completed fixed-topology DNA posterior run."""

    model_definition: FixedTopologyDnaModelDefinition
    proposal_schedule: FixedTopologyDnaProposalSchedule
    observation_policy: str
    chain_report: MetropolisHastingsRunReport
    posterior_rows: list[FixedTopologyDnaPosteriorRow]


def build_fixed_topology_dna_model_definition(
    *,
    substitution_model_name: str,
    branch_length_prior: BranchLengthPriorModel,
    substitution_parameter_prior_bundle: SubstitutionParameterPriorBundle,
    initial_kappa: float | None = None,
    initial_base_frequencies: Sequence[float] | dict[str, float] | None = None,
    initial_exchangeabilities: Sequence[float] | dict[str, float] | None = None,
) -> FixedTopologyDnaModelDefinition:
    """Build one validated fixed-topology DNA posterior model definition."""
    validated_substitution_model_name = _validate_substitution_model_name(
        substitution_model_name
    )
    if validated_substitution_model_name not in _SUPPORTED_FIXED_TOPOLOGY_DNA_MODELS:
        raise PhylogeneticsError(
            "fixed-topology DNA posterior model requires at least one sampled substitution parameter and does not support JC69",
            code="fixed_topology_dna_model_without_free_parameters",
            details={"substitution_model_name": validated_substitution_model_name},
        )
    if not isinstance(branch_length_prior, BranchLengthPriorModel):
        raise PhylogeneticsError(
            "fixed-topology DNA posterior model requires one BranchLengthPriorModel",
            code="fixed_topology_dna_branch_length_prior_type_invalid",
        )
    if branch_length_prior.family == "fixed":
        raise PhylogeneticsError(
            "fixed-topology DNA posterior model does not support fixed branch-length priors because Metropolis-Hastings proposals would leave finite prior support",
            code="fixed_topology_dna_branch_length_prior_family_invalid",
        )
    if not isinstance(
        substitution_parameter_prior_bundle,
        SubstitutionParameterPriorBundle,
    ):
        raise PhylogeneticsError(
            "fixed-topology DNA posterior model requires one SubstitutionParameterPriorBundle",
            code="fixed_topology_dna_substitution_prior_bundle_type_invalid",
        )
    _validate_supported_substitution_priors(
        substitution_model_name=validated_substitution_model_name,
        prior_bundle=substitution_parameter_prior_bundle,
    )
    resolved_initial_kappa = (
        float(
            format(
                validate_positive_kappa(
                    initial_kappa,
                    model_name="fixed-topology DNA posterior model",
                ),
                ".15g",
            )
        )
        if initial_kappa is not None
        else None
    )
    resolved_initial_base_frequencies = (
        parameterize_dna_base_frequency_simplex(
            initial_base_frequencies
        ).constrained_mapping()
        if initial_base_frequencies is not None
        else None
    )
    resolved_initial_exchangeabilities = (
        parameterize_dna_exchangeability_simplex(
            initial_exchangeabilities
        ).constrained_mapping()
        if initial_exchangeabilities is not None
        else None
    )
    active_targets = _ACTIVE_PARAMETER_TARGETS_BY_MODEL[
        validated_substitution_model_name
    ]
    if "kappa" not in active_targets and resolved_initial_kappa is not None:
        raise PhylogeneticsError(
            "fixed-topology DNA posterior model received one initial kappa value for a substitution model that does not use kappa",
            code="fixed_topology_dna_initial_kappa_unused",
            details={"substitution_model_name": validated_substitution_model_name},
        )
    if (
        "base-frequencies" not in active_targets
        and resolved_initial_base_frequencies is not None
    ):
        raise PhylogeneticsError(
            "fixed-topology DNA posterior model received one initial base-frequency vector for a substitution model that does not use base frequencies",
            code="fixed_topology_dna_initial_base_frequencies_unused",
            details={"substitution_model_name": validated_substitution_model_name},
        )
    if (
        "exchangeabilities" not in active_targets
        and resolved_initial_exchangeabilities is not None
    ):
        raise PhylogeneticsError(
            "fixed-topology DNA posterior model received one initial exchangeability vector for a substitution model that does not use exchangeabilities",
            code="fixed_topology_dna_initial_exchangeabilities_unused",
            details={"substitution_model_name": validated_substitution_model_name},
        )
    return FixedTopologyDnaModelDefinition(
        substitution_model_name=validated_substitution_model_name,
        branch_length_prior=branch_length_prior,
        substitution_parameter_prior_bundle=substitution_parameter_prior_bundle,
        initial_kappa=resolved_initial_kappa,
        initial_base_frequencies=resolved_initial_base_frequencies,
        initial_exchangeabilities=resolved_initial_exchangeabilities,
    )


def build_fixed_topology_dna_proposal_schedule(
    *,
    model_definition: FixedTopologyDnaModelDefinition,
    branch_length_move_weight: float,
    branch_length_log_scale_standard_deviation: float,
    kappa_move_weight: float = 0.0,
    kappa_log_scale_standard_deviation: float | None = None,
    base_frequency_move_weight: float = 0.0,
    base_frequency_coordinate_standard_deviation: float | None = None,
    exchangeability_move_weight: float = 0.0,
    exchangeability_coordinate_standard_deviation: float | None = None,
) -> FixedTopologyDnaProposalSchedule:
    """Build one validated proposal schedule for fixed-topology DNA sampling."""
    if not isinstance(model_definition, FixedTopologyDnaModelDefinition):
        raise PhylogeneticsError(
            "fixed-topology DNA proposal schedule requires one FixedTopologyDnaModelDefinition",
            code="fixed_topology_dna_proposal_schedule_model_definition_type_invalid",
        )
    validated_branch_length_move_weight = _validate_nonnegative_finite_float(
        value=branch_length_move_weight,
        field_name="branch_length_move_weight",
        owner_name="fixed-topology DNA proposal schedule",
    )
    if validated_branch_length_move_weight <= 0.0:
        raise PhylogeneticsError(
            "fixed-topology DNA proposal schedule requires 'branch_length_move_weight' to be greater than zero",
            code="fixed_topology_dna_branch_length_move_weight_invalid",
        )
    validated_branch_length_log_scale_standard_deviation = (
        _validate_positive_finite_float(
            value=branch_length_log_scale_standard_deviation,
            field_name="branch_length_log_scale_standard_deviation",
            owner_name="fixed-topology DNA proposal schedule",
        )
    )
    validated_kappa_move_weight = _validate_nonnegative_finite_float(
        value=kappa_move_weight,
        field_name="kappa_move_weight",
        owner_name="fixed-topology DNA proposal schedule",
    )
    validated_base_frequency_move_weight = _validate_nonnegative_finite_float(
        value=base_frequency_move_weight,
        field_name="base_frequency_move_weight",
        owner_name="fixed-topology DNA proposal schedule",
    )
    validated_exchangeability_move_weight = _validate_nonnegative_finite_float(
        value=exchangeability_move_weight,
        field_name="exchangeability_move_weight",
        owner_name="fixed-topology DNA proposal schedule",
    )
    active_targets = set(model_definition.active_parameter_targets)
    validated_kappa_log_scale_standard_deviation = (
        _validate_optional_positive_finite_float(
            value=kappa_log_scale_standard_deviation,
            field_name="kappa_log_scale_standard_deviation",
            owner_name="fixed-topology DNA proposal schedule",
        )
    )
    validated_base_frequency_coordinate_standard_deviation = (
        _validate_optional_positive_finite_float(
            value=base_frequency_coordinate_standard_deviation,
            field_name="base_frequency_coordinate_standard_deviation",
            owner_name="fixed-topology DNA proposal schedule",
        )
    )
    validated_exchangeability_coordinate_standard_deviation = (
        _validate_optional_positive_finite_float(
            value=exchangeability_coordinate_standard_deviation,
            field_name="exchangeability_coordinate_standard_deviation",
            owner_name="fixed-topology DNA proposal schedule",
        )
    )
    _validate_parameter_move_activation(
        target_name="kappa",
        active_targets=active_targets,
        move_weight=validated_kappa_move_weight,
        standard_deviation=validated_kappa_log_scale_standard_deviation,
    )
    _validate_parameter_move_activation(
        target_name="base-frequencies",
        active_targets=active_targets,
        move_weight=validated_base_frequency_move_weight,
        standard_deviation=validated_base_frequency_coordinate_standard_deviation,
    )
    _validate_parameter_move_activation(
        target_name="exchangeabilities",
        active_targets=active_targets,
        move_weight=validated_exchangeability_move_weight,
        standard_deviation=validated_exchangeability_coordinate_standard_deviation,
    )
    return FixedTopologyDnaProposalSchedule(
        substitution_model_name=model_definition.substitution_model_name,
        branch_length_move_weight=validated_branch_length_move_weight,
        branch_length_log_scale_standard_deviation=(
            validated_branch_length_log_scale_standard_deviation
        ),
        kappa_move_weight=validated_kappa_move_weight,
        kappa_log_scale_standard_deviation=validated_kappa_log_scale_standard_deviation,
        base_frequency_move_weight=validated_base_frequency_move_weight,
        base_frequency_coordinate_standard_deviation=(
            validated_base_frequency_coordinate_standard_deviation
        ),
        exchangeability_move_weight=validated_exchangeability_move_weight,
        exchangeability_coordinate_standard_deviation=(
            validated_exchangeability_coordinate_standard_deviation
        ),
    )


def run_fixed_topology_dna_metropolis_hastings(
    *,
    tree: PhyloTree,
    records: Sequence[AlignmentRecord],
    model_definition: FixedTopologyDnaModelDefinition,
    proposal_schedule: FixedTopologyDnaProposalSchedule,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
    observation_policy: str = "reject",
) -> FixedTopologyDnaRunReport:
    """Run one fixed-topology DNA posterior sampler over branch lengths and model parameters."""
    if not isinstance(tree, PhyloTree):
        raise PhylogeneticsError(
            "fixed-topology DNA posterior runner requires one PhyloTree",
            code="fixed_topology_dna_tree_type_invalid",
        )
    if not isinstance(model_definition, FixedTopologyDnaModelDefinition):
        raise PhylogeneticsError(
            "fixed-topology DNA posterior runner requires one FixedTopologyDnaModelDefinition",
            code="fixed_topology_dna_model_definition_type_invalid",
        )
    if not isinstance(proposal_schedule, FixedTopologyDnaProposalSchedule):
        raise PhylogeneticsError(
            "fixed-topology DNA posterior runner requires one FixedTopologyDnaProposalSchedule",
            code="fixed_topology_dna_proposal_schedule_type_invalid",
        )
    if (
        proposal_schedule.substitution_model_name
        != model_definition.substitution_model_name
    ):
        raise PhylogeneticsError(
            "fixed-topology DNA posterior runner requires the proposal schedule and model definition to use the same substitution model",
            code="fixed_topology_dna_model_schedule_mismatch",
            details={
                "model_definition": model_definition.substitution_model_name,
                "proposal_schedule": proposal_schedule.substitution_model_name,
            },
        )
    validated_observation_policy = validate_dna_observation_policy(
        observation_policy,
        owner_name="fixed-topology DNA posterior runner",
    )
    normalized_records = normalize_dna_likelihood_records(
        list(records),
        model_name=model_definition.substitution_model_name,
        observation_policy=validated_observation_policy,
    )
    fixed_tree = tree.copy()
    fixed_tree.rooted = tree.rooted
    initial_model_parameters = _build_initial_model_parameters(
        model_definition=model_definition,
        records=normalized_records,
        observation_policy=validated_observation_policy,
    )
    initial_state = score_bayesian_phylogenetic_state(
        tree=fixed_tree,
        model_parameters=initial_model_parameters,
        update_prior_components=lambda state: (
            _build_fixed_topology_dna_prior_components(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=None,
            )
        ),
        update_log_likelihood=lambda state: _evaluate_fixed_topology_dna_log_likelihood(
            state=state,
            records=normalized_records,
            model_definition=model_definition,
            fixed_topology_id=None,
            observation_policy=validated_observation_policy,
        ),
    )
    fixed_topology_id = initial_state.tree.topology_id
    chain_report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: _propose_fixed_topology_dna_state(
            current_state=current_state,
            rng=rng,
            proposal_schedule=proposal_schedule,
        ),
        update_prior_components=lambda state: (
            _build_fixed_topology_dna_prior_components(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=fixed_topology_id,
            )
        ),
        update_log_likelihood=lambda state: _evaluate_fixed_topology_dna_log_likelihood(
            state=state,
            records=normalized_records,
            model_definition=model_definition,
            fixed_topology_id=fixed_topology_id,
            observation_policy=validated_observation_policy,
        ),
        iteration_count=iteration_count,
        sample_every=sample_every,
        seed=seed,
    )
    posterior_rows = _build_fixed_topology_dna_posterior_rows(
        chain_report=chain_report,
        fixed_topology_id=fixed_topology_id,
        substitution_model_name=model_definition.substitution_model_name,
    )
    return FixedTopologyDnaRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        observation_policy=validated_observation_policy,
        chain_report=chain_report,
        posterior_rows=posterior_rows,
    )


def _build_initial_model_parameters(
    *,
    model_definition: FixedTopologyDnaModelDefinition,
    records: list[AlignmentRecord],
    observation_policy: str,
) -> object:
    scalar_parameters: dict[str, float] = {}
    vector_parameters: dict[str, dict[str, float]] = {}
    if "kappa" in model_definition.active_parameter_targets:
        scalar_parameters["kappa"] = (
            model_definition.initial_kappa or _DEFAULT_INITIAL_KAPPA
        )
    if "base-frequencies" in model_definition.active_parameter_targets:
        vector_parameters["base-frequencies"] = (
            model_definition.initial_base_frequencies
            or parameterize_dna_base_frequency_simplex(
                _stabilize_positive_simplex_mapping(
                    estimate_empirical_dna_base_frequencies_from_records(
                        records,
                        model_name=model_definition.substitution_model_name,
                        observation_policy=observation_policy,
                    ),
                )
            ).constrained_mapping()
        )
    if "exchangeabilities" in model_definition.active_parameter_targets:
        vector_parameters["exchangeabilities"] = (
            model_definition.initial_exchangeabilities
            or parameterize_dna_exchangeability_simplex(
                _DEFAULT_INITIAL_EXCHANGEABILITIES
            ).constrained_mapping()
        )
    return build_bayesian_model_parameter_state(
        categorical_parameters={
            "substitution-model": model_definition.substitution_model_name
        },
        scalar_parameters=scalar_parameters,
        vector_parameters=vector_parameters,
    )


def _build_fixed_topology_dna_prior_components(
    *,
    state: BayesianPhylogeneticState,
    model_definition: FixedTopologyDnaModelDefinition,
    fixed_topology_id: str | None,
) -> list[BayesianPriorComponentState]:
    model_name = _require_fixed_topology_dna_state_consistency(
        state=state,
        model_definition=model_definition,
        fixed_topology_id=fixed_topology_id,
    )
    branch_length_prior_report = evaluate_tree_branch_length_log_prior(
        state.tree.to_tree(),
        model_definition.branch_length_prior,
    )
    substitution_parameter_prior_report = evaluate_substitution_parameter_log_prior(
        prior_bundle=model_definition.substitution_parameter_prior_bundle,
        kappa=state.model_parameters.scalar_parameters.get("kappa"),
        base_frequencies=state.model_parameters.vector_parameters.get(
            "base-frequencies"
        ),
        exchangeabilities=state.model_parameters.vector_parameters.get(
            "exchangeabilities"
        ),
    )
    prior_components = [
        build_bayesian_prior_component_state(
            component_name="branch-lengths",
            family=branch_length_prior_report.family,
            log_prior=branch_length_prior_report.total_log_prior,
            parameter_values=branch_length_prior_report.parameter_values,
        )
    ]
    prior_components.extend(
        build_bayesian_prior_component_state(
            component_name=f"substitution:{row.target_name}",
            family=row.family,
            log_prior=row.log_prior_contribution,
            parameter_values=row.hyperparameter_values,
        )
        for row in substitution_parameter_prior_report.rows
    )
    if model_name not in _SUPPORTED_FIXED_TOPOLOGY_DNA_MODELS:
        raise AssertionError(f"unsupported fixed-topology DNA model {model_name}")
    return prior_components


def _evaluate_fixed_topology_dna_log_likelihood(
    *,
    state: BayesianPhylogeneticState,
    records: list[AlignmentRecord],
    model_definition: FixedTopologyDnaModelDefinition,
    fixed_topology_id: str | None,
    observation_policy: str,
) -> float:
    model_name = _require_fixed_topology_dna_state_consistency(
        state=state,
        model_definition=model_definition,
        fixed_topology_id=fixed_topology_id,
    )
    tree = state.tree.to_tree()
    if model_name == "K80":
        return evaluate_k80_tree_likelihood(
            tree,
            records,
            kappa=state.model_parameters.scalar_parameters["kappa"],
            observation_policy=observation_policy,
        ).log_likelihood
    if model_name == "HKY85":
        return evaluate_hky85_tree_likelihood(
            tree,
            records,
            kappa=state.model_parameters.scalar_parameters["kappa"],
            base_frequencies=state.model_parameters.vector_parameters[
                "base-frequencies"
            ],
            observation_policy=observation_policy,
        ).log_likelihood
    if model_name == "GTR":
        return evaluate_gtr_tree_likelihood(
            tree,
            records,
            exchangeabilities=state.model_parameters.vector_parameters[
                "exchangeabilities"
            ],
            base_frequencies=state.model_parameters.vector_parameters[
                "base-frequencies"
            ],
            observation_policy=observation_policy,
        ).log_likelihood
    raise AssertionError(f"unsupported fixed-topology DNA model {model_name}")


def _propose_fixed_topology_dna_state(
    *,
    current_state: BayesianPhylogeneticState,
    rng,
    proposal_schedule: FixedTopologyDnaProposalSchedule,
):
    weighted_moves = [
        (
            proposal_schedule.branch_length_move_weight,
            lambda: propose_branch_length_scaling_move(
                current_state,
                rng,
                log_scale_standard_deviation=(
                    proposal_schedule.branch_length_log_scale_standard_deviation
                ),
            ),
        )
    ]
    if proposal_schedule.kappa_move_weight > 0.0:
        kappa_log_scale_standard_deviation = require_present(
            proposal_schedule.kappa_log_scale_standard_deviation,
            owner_name="fixed-topology DNA proposal schedule",
            field_name="kappa_log_scale_standard_deviation",
        )
        weighted_moves.append(
            (
                proposal_schedule.kappa_move_weight,
                lambda: propose_kappa_move(
                    current_state,
                    rng,
                    log_scale_standard_deviation=kappa_log_scale_standard_deviation,
                ),
            )
        )
    if proposal_schedule.base_frequency_move_weight > 0.0:
        base_frequency_coordinate_standard_deviation = require_present(
            proposal_schedule.base_frequency_coordinate_standard_deviation,
            owner_name="fixed-topology DNA proposal schedule",
            field_name="base_frequency_coordinate_standard_deviation",
        )
        weighted_moves.append(
            (
                proposal_schedule.base_frequency_move_weight,
                lambda: propose_base_frequency_simplex_move(
                    current_state,
                    rng,
                    unconstrained_coordinate_standard_deviation=(
                        base_frequency_coordinate_standard_deviation
                    ),
                ),
            )
        )
    if proposal_schedule.exchangeability_move_weight > 0.0:
        exchangeability_coordinate_standard_deviation = require_present(
            proposal_schedule.exchangeability_coordinate_standard_deviation,
            owner_name="fixed-topology DNA proposal schedule",
            field_name="exchangeability_coordinate_standard_deviation",
        )
        weighted_moves.append(
            (
                proposal_schedule.exchangeability_move_weight,
                lambda: propose_gtr_exchangeability_move(
                    current_state,
                    rng,
                    unconstrained_coordinate_standard_deviation=(
                        exchangeability_coordinate_standard_deviation
                    ),
                ),
            )
        )
    total_weight = math.fsum(weight for weight, _move in weighted_moves)
    move_threshold = rng.random() * total_weight
    cumulative_weight = 0.0
    for weight, move in weighted_moves:
        cumulative_weight += weight
        if move_threshold <= cumulative_weight:
            return move()
    return weighted_moves[-1][1]()


def _build_fixed_topology_dna_posterior_rows(
    *,
    chain_report: MetropolisHastingsRunReport,
    fixed_topology_id: str,
    substitution_model_name: str,
) -> list[FixedTopologyDnaPosteriorRow]:
    posterior_rows: list[FixedTopologyDnaPosteriorRow] = []
    for sample_index, state in enumerate(chain_report.sampled_states):
        if state.tree.topology_id != fixed_topology_id:
            raise PhylogeneticsError(
                "fixed-topology DNA posterior trace detected one topology change in sampled states",
                code="fixed_topology_dna_trace_topology_changed",
                details={
                    "expected_topology_id": fixed_topology_id,
                    "observed_topology_id": state.tree.topology_id,
                    "sample_index": sample_index,
                },
            )
        posterior_rows.append(
            FixedTopologyDnaPosteriorRow(
                sample_index=sample_index,
                iteration_index=sample_index * chain_report.sample_every,
                topology_id=state.tree.topology_id,
                substitution_model_name=substitution_model_name,
                total_log_prior=state.total_log_prior,
                log_likelihood=state.log_likelihood,
                posterior_log_score=state.posterior_log_score,
                prior_component_log_priors={
                    component.component_name: component.log_prior
                    for component in state.prior_components
                },
                branch_lengths={
                    branch_row.branch_id: branch_row.branch_length
                    for branch_row in state.tree.branch_rows
                },
                scalar_parameters=dict(state.model_parameters.scalar_parameters),
                vector_parameters={
                    parameter_name: dict(component_values)
                    for parameter_name, component_values in state.model_parameters.vector_parameters.items()
                },
            )
        )
    return posterior_rows


def _require_fixed_topology_dna_state_consistency(
    *,
    state: BayesianPhylogeneticState,
    model_definition: FixedTopologyDnaModelDefinition,
    fixed_topology_id: str | None,
) -> str:
    model_name = state.model_parameters.categorical_parameters.get("substitution-model")
    if model_name != model_definition.substitution_model_name:
        raise PhylogeneticsError(
            "fixed-topology DNA posterior model requires every sampled state to preserve the configured substitution-model label",
            code="fixed_topology_dna_state_model_label_invalid",
            details={
                "expected_model_name": model_definition.substitution_model_name,
                "observed_model_name": model_name,
            },
        )
    if fixed_topology_id is not None and state.tree.topology_id != fixed_topology_id:
        raise PhylogeneticsError(
            "fixed-topology DNA posterior model requires topology to remain unchanged across sampled states",
            code="fixed_topology_dna_state_topology_changed",
            details={
                "expected_topology_id": fixed_topology_id,
                "observed_topology_id": state.tree.topology_id,
            },
        )
    return model_name


def _validate_supported_substitution_priors(
    *,
    substitution_model_name: str,
    prior_bundle: SubstitutionParameterPriorBundle,
) -> None:
    expected_targets = set(_ACTIVE_PARAMETER_TARGETS_BY_MODEL[substitution_model_name])
    provided_targets = {
        target_name
        for target_name, prior_model in (
            ("kappa", prior_bundle.kappa_prior),
            ("exchangeabilities", prior_bundle.exchangeability_prior),
            ("base-frequencies", prior_bundle.base_frequency_prior),
            ("gamma-alpha", prior_bundle.gamma_alpha_prior),
            ("invariant-proportion", prior_bundle.invariant_proportion_prior),
        )
        if prior_model is not None
    }
    unsupported_targets = provided_targets - expected_targets
    if unsupported_targets:
        raise PhylogeneticsError(
            "fixed-topology DNA posterior model received substitution priors for parameters that the chosen substitution model does not sample",
            code="fixed_topology_dna_prior_targets_invalid",
            details={
                "substitution_model_name": substitution_model_name,
                "unsupported_targets": sorted(unsupported_targets),
            },
        )
    missing_targets = expected_targets - provided_targets
    if missing_targets:
        raise PhylogeneticsError(
            "fixed-topology DNA posterior model requires explicit substitution priors for every sampled substitution parameter",
            code="fixed_topology_dna_prior_targets_missing",
            details={
                "substitution_model_name": substitution_model_name,
                "missing_targets": sorted(missing_targets),
            },
        )
    for target_name, prior_model in (
        ("kappa", prior_bundle.kappa_prior),
        ("exchangeabilities", prior_bundle.exchangeability_prior),
        ("base-frequencies", prior_bundle.base_frequency_prior),
    ):
        if prior_model is not None and prior_model.family == "fixed":
            raise PhylogeneticsError(
                "fixed-topology DNA posterior model does not support fixed substitution priors because Metropolis-Hastings proposals would leave finite prior support",
                code="fixed_topology_dna_prior_family_invalid",
                details={
                    "substitution_model_name": substitution_model_name,
                    "target_name": target_name,
                },
            )


def _validate_parameter_move_activation(
    *,
    target_name: str,
    active_targets: set[str],
    move_weight: float,
    standard_deviation: float | None,
) -> None:
    if target_name in active_targets:
        if move_weight <= 0.0:
            raise PhylogeneticsError(
                "fixed-topology DNA proposal schedule requires a positive move weight for every sampled substitution parameter",
                code="fixed_topology_dna_parameter_move_weight_invalid",
                details={"target_name": target_name},
            )
        if standard_deviation is None:
            raise PhylogeneticsError(
                "fixed-topology DNA proposal schedule requires one finite positive proposal scale for every sampled substitution parameter",
                code="fixed_topology_dna_parameter_move_scale_missing",
                details={"target_name": target_name},
            )
        return
    if move_weight > 0.0 or standard_deviation is not None:
        raise PhylogeneticsError(
            "fixed-topology DNA proposal schedule cannot configure proposals for substitution parameters that the chosen substitution model does not sample",
            code="fixed_topology_dna_parameter_move_unused",
            details={"target_name": target_name},
        )


def _validate_substitution_model_name(value: str) -> str:
    normalized_value = value.strip().upper()
    if normalized_value not in FIXED_TOPOLOGY_DNA_SUBSTITUTION_MODELS:
        raise PhylogeneticsError(
            "fixed-topology DNA posterior model requires one supported substitution model name",
            code="fixed_topology_dna_substitution_model_name_invalid",
            details={
                "substitution_model_name": value,
                "allowed_models": list(FIXED_TOPOLOGY_DNA_SUBSTITUTION_MODELS),
            },
        )
    return normalized_value


def _validate_nonnegative_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    try:
        validated_value = float(value)
    except (TypeError, ValueError) as error:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be a finite float",
            code="fixed_topology_dna_float_type_invalid",
            details={"field_name": field_name},
        ) from error
    if not math.isfinite(validated_value) or validated_value < 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be nonnegative and finite",
            code="fixed_topology_dna_nonnegative_float_invalid",
            details={"field_name": field_name, "value": value},
        )
    return float(format(validated_value, ".15g"))


def _validate_positive_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    validated_value = _validate_nonnegative_finite_float(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )
    if validated_value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be greater than zero",
            code="fixed_topology_dna_positive_float_invalid",
            details={"field_name": field_name, "value": value},
        )
    return validated_value


def _validate_optional_positive_finite_float(
    *,
    value: float | None,
    field_name: str,
    owner_name: str,
) -> float | None:
    if value is None:
        return None
    return _validate_positive_finite_float(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )


def _stabilize_positive_simplex_mapping(
    raw_values: Sequence[float],
) -> dict[str, float]:
    stabilized_values = [
        max(float(component_value), _MINIMUM_SIMPLEX_COMPONENT)
        for component_value in raw_values
    ]
    total = math.fsum(stabilized_values)
    return {
        component_name: float(format(component_value / total, ".15g"))
        for component_name, component_value in zip(
            ("A", "C", "G", "T"),
            stabilized_values,
            strict=True,
        )
    }


__all__ = [
    "FIXED_TOPOLOGY_DNA_SUBSTITUTION_MODELS",
    "FixedTopologyDnaModelDefinition",
    "FixedTopologyDnaPosteriorRow",
    "FixedTopologyDnaProposalSchedule",
    "FixedTopologyDnaRunReport",
    "build_fixed_topology_dna_model_definition",
    "build_fixed_topology_dna_proposal_schedule",
    "run_fixed_topology_dna_metropolis_hastings",
]
