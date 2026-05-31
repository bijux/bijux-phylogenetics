from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math
import random

from bijux_phylogenetics.bayesian.branch_length_priors import (
    BranchLengthPriorModel,
    evaluate_tree_branch_length_log_prior,
)
from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    PositiveSubstitutionParameterPriorModel,
    ProbabilitySubstitutionParameterPriorModel,
    SimplexSubstitutionParameterPriorModel,
    SubstitutionParameterPriorBundle,
    SubstitutionParameterPriorRow,
    evaluate_substitution_parameter_log_prior,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    TreeTopologyPriorModel,
    evaluate_tree_topology_log_prior,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

_DNA_BASE_FREQUENCY_COMPONENT_NAMES = ("A", "C", "G", "T")


@dataclass(frozen=True, slots=True)
class PriorOnlySampledBranchRow:
    """One sampled branch-length row from one prior-only phylogenetic state."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    branch_length: float
    log_prior_contribution: float


@dataclass(frozen=True, slots=True)
class PriorOnlySubstitutionParameterState:
    """One sampled substitution-parameter state drawn directly from priors."""

    kappa: float | None = None
    exchangeabilities: dict[str, float] | None = None
    base_frequencies: dict[str, float] | None = None
    gamma_alpha: float | None = None
    invariant_proportion: float | None = None


@dataclass(frozen=True, slots=True)
class PriorOnlyPhylogeneticSample:
    """One complete prior-only phylogenetic state sample."""

    sample_index: int
    tree_newick: str
    topology_id: str
    tree_topology_prior_family: str
    branch_length_prior_family: str
    topology_log_prior: float
    branch_length_log_prior: float
    substitution_log_prior: float
    total_log_prior: float
    branch_rows: list[PriorOnlySampledBranchRow]
    substitution_prior_rows: list[SubstitutionParameterPriorRow]
    substitution_parameter_state: PriorOnlySubstitutionParameterState


@dataclass(frozen=True, slots=True)
class PriorOnlyPhylogeneticSimulationReport:
    """One seeded prior-only simulation run over phylogenetic states."""

    sample_count: int
    seed: int
    tree_topology_prior_family: str
    branch_length_prior_family: str
    substitution_prior_count: int
    samples: list[PriorOnlyPhylogeneticSample]


def sample_prior_only_phylogenetic_state(
    *,
    tree_topology_prior: TreeTopologyPriorModel,
    branch_length_prior: BranchLengthPriorModel,
    substitution_parameter_prior: SubstitutionParameterPriorBundle | None = None,
    seed: int = 0,
) -> PriorOnlyPhylogeneticSample:
    """Sample one complete phylogenetic state from priors without likelihood data."""
    simulation_report = simulate_prior_only_phylogenetic_states(
        tree_topology_prior=tree_topology_prior,
        branch_length_prior=branch_length_prior,
        substitution_parameter_prior=substitution_parameter_prior,
        sample_count=1,
        seed=seed,
    )
    return simulation_report.samples[0]


def simulate_prior_only_phylogenetic_states(
    *,
    tree_topology_prior: TreeTopologyPriorModel,
    branch_length_prior: BranchLengthPriorModel,
    substitution_parameter_prior: SubstitutionParameterPriorBundle | None = None,
    sample_count: int = 1,
    seed: int = 0,
) -> PriorOnlyPhylogeneticSimulationReport:
    """Sample seeded phylogenetic states directly from configured priors."""
    validated_sample_count = _validate_prior_only_sample_count(sample_count)
    validated_seed = _validate_prior_only_seed(seed)
    rng = random.Random(validated_seed)  # nosec B311
    samples = [
        _sample_prior_only_phylogenetic_state(
            sample_index=sample_index,
            tree_topology_prior=tree_topology_prior,
            branch_length_prior=branch_length_prior,
            substitution_parameter_prior=substitution_parameter_prior,
            rng=rng,
        )
        for sample_index in range(1, validated_sample_count + 1)
    ]
    return PriorOnlyPhylogeneticSimulationReport(
        sample_count=validated_sample_count,
        seed=validated_seed,
        tree_topology_prior_family=tree_topology_prior.family,
        branch_length_prior_family=branch_length_prior.family,
        substitution_prior_count=(
            substitution_parameter_prior.prior_count()
            if substitution_parameter_prior is not None
            else 0
        ),
        samples=samples,
    )


def _sample_prior_only_phylogenetic_state(
    *,
    sample_index: int,
    tree_topology_prior: TreeTopologyPriorModel,
    branch_length_prior: BranchLengthPriorModel,
    substitution_parameter_prior: SubstitutionParameterPriorBundle | None,
    rng: random.Random,
) -> PriorOnlyPhylogeneticSample:
    sampled_tree = _sample_rooted_labeled_bifurcating_tree(
        taxa=tree_topology_prior.taxa,
        rng=rng,
    )
    _assign_sampled_branch_lengths(
        tree=sampled_tree,
        prior_model=branch_length_prior,
        rng=rng,
    )
    topology_report = evaluate_tree_topology_log_prior(
        sampled_tree, tree_topology_prior
    )
    branch_length_report = evaluate_tree_branch_length_log_prior(
        sampled_tree,
        branch_length_prior,
    )
    substitution_parameter_state = PriorOnlySubstitutionParameterState()
    substitution_prior_rows: list[SubstitutionParameterPriorRow] = []
    substitution_log_prior = 0.0
    if substitution_parameter_prior is not None:
        substitution_parameter_state = _sample_substitution_parameter_state(
            prior_bundle=substitution_parameter_prior,
            rng=rng,
        )
        substitution_report = evaluate_substitution_parameter_log_prior(
            prior_bundle=substitution_parameter_prior,
            kappa=substitution_parameter_state.kappa,
            exchangeabilities=substitution_parameter_state.exchangeabilities,
            base_frequencies=substitution_parameter_state.base_frequencies,
            gamma_alpha=substitution_parameter_state.gamma_alpha,
            invariant_proportion=substitution_parameter_state.invariant_proportion,
        )
        substitution_prior_rows = list(substitution_report.rows)
        substitution_log_prior = substitution_report.total_log_prior
    branch_rows = [
        PriorOnlySampledBranchRow(
            branch_id=row.branch_id,
            child_name=row.child_name,
            descendant_taxa=list(row.descendant_taxa),
            branch_length=row.branch_length,
            log_prior_contribution=row.log_density,
        )
        for row in branch_length_report.branch_rows
    ]
    total_log_prior = _round_float(
        math.fsum(
            [
                topology_report.log_prior,
                branch_length_report.total_log_prior,
                substitution_log_prior,
            ]
        )
    )
    return PriorOnlyPhylogeneticSample(
        sample_index=sample_index,
        tree_newick=sampled_tree.to_newick(),
        topology_id=topology_report.topology_id,
        tree_topology_prior_family=topology_report.family,
        branch_length_prior_family=branch_length_report.family,
        topology_log_prior=topology_report.log_prior,
        branch_length_log_prior=branch_length_report.total_log_prior,
        substitution_log_prior=substitution_log_prior,
        total_log_prior=total_log_prior,
        branch_rows=branch_rows,
        substitution_prior_rows=substitution_prior_rows,
        substitution_parameter_state=substitution_parameter_state,
    )


def _sample_rooted_labeled_bifurcating_tree(
    *,
    taxa: Sequence[str],
    rng: random.Random,
) -> PhyloTree:
    if len(taxa) < 2:
        raise PhylogeneticsError(
            "prior-only tree sampling requires at least two taxa",
            code="prior_only_tree_sampling_requires_two_or_more_taxa",
        )
    tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(name=taxa[0]),
                TreeNode(name=taxa[1]),
            ]
        ),
        rooted=True,
    )
    for taxon in taxa[2:]:
        insertion_positions = [None, *list(tree.iter_edges())]
        chosen_position = insertion_positions[rng.randrange(len(insertion_positions))]
        if chosen_position is None:
            tree.root = TreeNode(children=[tree.root, TreeNode(name=taxon)])
            tree.refresh()
            continue
        parent, child = chosen_position
        inserted_internal_node = TreeNode(children=[child, TreeNode(name=taxon)])
        child_index = parent.children.index(child)
        parent.children[child_index] = inserted_internal_node
        inserted_internal_node.parent = parent
        tree.refresh()
    return tree


def _assign_sampled_branch_lengths(
    *,
    tree: PhyloTree,
    prior_model: BranchLengthPriorModel,
    rng: random.Random,
) -> None:
    for node in tree.iter_nodes(order="preorder"):
        if node is tree.root:
            continue
        node.branch_length = _round_float(
            _sample_branch_length_value(prior_model=prior_model, rng=rng)
        )


def _sample_substitution_parameter_state(
    *,
    prior_bundle: SubstitutionParameterPriorBundle,
    rng: random.Random,
) -> PriorOnlySubstitutionParameterState:
    kappa = (
        _round_float(_sample_positive_prior_value(prior_bundle.kappa_prior, rng=rng))
        if prior_bundle.kappa_prior is not None
        else None
    )
    exchangeabilities = (
        _sample_simplex_prior_values(prior_bundle.exchangeability_prior, rng=rng)
        if prior_bundle.exchangeability_prior is not None
        else None
    )
    base_frequencies = (
        _sample_simplex_prior_values(prior_bundle.base_frequency_prior, rng=rng)
        if prior_bundle.base_frequency_prior is not None
        else None
    )
    gamma_alpha = (
        _round_float(
            _sample_positive_prior_value(prior_bundle.gamma_alpha_prior, rng=rng)
        )
        if prior_bundle.gamma_alpha_prior is not None
        else None
    )
    invariant_proportion = (
        _round_float(
            _sample_probability_prior_value(
                prior_bundle.invariant_proportion_prior,
                rng=rng,
            )
        )
        if prior_bundle.invariant_proportion_prior is not None
        else None
    )
    return PriorOnlySubstitutionParameterState(
        kappa=kappa,
        exchangeabilities=exchangeabilities,
        base_frequencies=base_frequencies,
        gamma_alpha=gamma_alpha,
        invariant_proportion=invariant_proportion,
    )


def _sample_branch_length_value(
    *,
    prior_model: BranchLengthPriorModel,
    rng: random.Random,
) -> float:
    if prior_model.family == "exponential":
        rate = require_present(
            prior_model.rate,
            owner_name="prior-only branch-length sampling",
            field_name="rate",
        )
        return rng.expovariate(rate)
    if prior_model.family == "gamma":
        shape = require_present(
            prior_model.shape,
            owner_name="prior-only branch-length sampling",
            field_name="shape",
        )
        scale = require_present(
            prior_model.scale,
            owner_name="prior-only branch-length sampling",
            field_name="scale",
        )
        return rng.gammavariate(shape, scale)
    if prior_model.family == "lognormal":
        log_mean = require_present(
            prior_model.log_mean,
            owner_name="prior-only branch-length sampling",
            field_name="log_mean",
        )
        log_standard_deviation = require_present(
            prior_model.log_standard_deviation,
            owner_name="prior-only branch-length sampling",
            field_name="log_standard_deviation",
        )
        return rng.lognormvariate(
            log_mean,
            log_standard_deviation,
        )
    if prior_model.family == "fixed":
        return require_present(
            prior_model.fixed_value,
            owner_name="prior-only branch-length sampling",
            field_name="fixed_value",
        )
    raise PhylogeneticsError(
        "prior-only simulation does not support the configured branch-length prior family",
        code="prior_only_branch_length_prior_family_invalid",
        details={"family": prior_model.family},
    )


def _sample_positive_prior_value(
    prior_model: PositiveSubstitutionParameterPriorModel,
    *,
    rng: random.Random,
) -> float:
    if prior_model.family == "exponential":
        rate = require_present(
            prior_model.rate,
            owner_name="prior-only positive substitution-parameter sampling",
            field_name="rate",
        )
        return rng.expovariate(rate)
    if prior_model.family == "gamma":
        shape = require_present(
            prior_model.shape,
            owner_name="prior-only positive substitution-parameter sampling",
            field_name="shape",
        )
        scale = require_present(
            prior_model.scale,
            owner_name="prior-only positive substitution-parameter sampling",
            field_name="scale",
        )
        return rng.gammavariate(shape, scale)
    if prior_model.family == "lognormal":
        log_mean = require_present(
            prior_model.log_mean,
            owner_name="prior-only positive substitution-parameter sampling",
            field_name="log_mean",
        )
        log_standard_deviation = require_present(
            prior_model.log_standard_deviation,
            owner_name="prior-only positive substitution-parameter sampling",
            field_name="log_standard_deviation",
        )
        return rng.lognormvariate(
            log_mean,
            log_standard_deviation,
        )
    if prior_model.family == "fixed":
        return require_present(
            prior_model.fixed_value,
            owner_name="prior-only positive substitution-parameter sampling",
            field_name="fixed_value",
        )
    raise PhylogeneticsError(
        "prior-only simulation does not support the configured positive substitution-parameter prior family",
        code="prior_only_positive_substitution_parameter_prior_family_invalid",
        details={"family": prior_model.family},
    )


def _sample_probability_prior_value(
    prior_model: ProbabilitySubstitutionParameterPriorModel,
    *,
    rng: random.Random,
) -> float:
    if prior_model.family == "beta":
        alpha = require_present(
            prior_model.alpha,
            owner_name="prior-only probability substitution-parameter sampling",
            field_name="alpha",
        )
        beta = require_present(
            prior_model.beta,
            owner_name="prior-only probability substitution-parameter sampling",
            field_name="beta",
        )
        numerator = rng.gammavariate(alpha, 1.0)
        denominator_remainder = rng.gammavariate(beta, 1.0)
        denominator = numerator + denominator_remainder
        if denominator <= 0.0:
            raise PhylogeneticsError(
                "beta prior sampling produced one non-positive normalization constant",
                code="prior_only_beta_prior_normalization_invalid",
            )
        return numerator / denominator
    if prior_model.family == "fixed":
        return require_present(
            prior_model.fixed_value,
            owner_name="prior-only probability substitution-parameter sampling",
            field_name="fixed_value",
        )
    raise PhylogeneticsError(
        "prior-only simulation does not support the configured probability substitution-parameter prior family",
        code="prior_only_probability_substitution_parameter_prior_family_invalid",
        details={"family": prior_model.family},
    )


def _sample_simplex_prior_values(
    prior_model: SimplexSubstitutionParameterPriorModel,
    *,
    rng: random.Random,
) -> dict[str, float]:
    if prior_model.family == "fixed":
        fixed_values = require_present(
            prior_model.fixed_values,
            owner_name="prior-only simplex substitution-parameter sampling",
            field_name="fixed_values",
        )
        return {
            component_name: _round_float(component_value)
            for component_name, component_value in zip(
                prior_model.component_names,
                fixed_values,
                strict=True,
            )
        }
    if prior_model.family == "dirichlet":
        concentration_parameters = require_present(
            prior_model.concentration_parameters,
            owner_name="prior-only simplex substitution-parameter sampling",
            field_name="concentration_parameters",
        )
        gamma_draws = [
            rng.gammavariate(concentration_parameter, 1.0)
            for concentration_parameter in concentration_parameters
        ]
        total = math.fsum(gamma_draws)
        if total <= 0.0:
            raise PhylogeneticsError(
                "dirichlet prior sampling produced one non-positive normalization constant",
                code="prior_only_dirichlet_prior_normalization_invalid",
            )
        return {
            component_name: _round_float(component_value / total)
            for component_name, component_value in zip(
                prior_model.component_names,
                gamma_draws,
                strict=True,
            )
        }
    raise PhylogeneticsError(
        "prior-only simulation does not support the configured simplex substitution-parameter prior family",
        code="prior_only_simplex_substitution_parameter_prior_family_invalid",
        details={"family": prior_model.family},
    )


def _validate_prior_only_sample_count(sample_count: int) -> int:
    if isinstance(sample_count, bool) or not isinstance(sample_count, int):
        raise PhylogeneticsError(
            "prior-only simulation requires sample_count to be one integer",
            code="prior_only_simulation_sample_count_type_invalid",
        )
    if sample_count <= 0:
        raise PhylogeneticsError(
            "prior-only simulation requires sample_count to be positive",
            code="prior_only_simulation_sample_count_invalid",
        )
    return sample_count


def _validate_prior_only_seed(seed: int) -> int:
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise PhylogeneticsError(
            "prior-only simulation requires seed to be one integer",
            code="prior_only_simulation_seed_type_invalid",
        )
    return seed


def _round_float(value: float) -> float:
    return float(format(value, ".15g"))
