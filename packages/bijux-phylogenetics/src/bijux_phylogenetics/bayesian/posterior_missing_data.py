from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math
from statistics import mean, median

import numpy

from bijux_phylogenetics.ancestral.discrete.likelihood.likelihood_math import (
    branch_length as discrete_trait_branch_length,
)
from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.comparative._math import invert_matrix, stable_covariance
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    iterate_partition_sites,
    slice_partition_sequence,
)
from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_STATE_ORDER,
    UNIFORM_DNA_ROOT_PRIOR,
)
from bijux_phylogenetics.phylo.likelihood.f81 import f81_transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.gtr import gtr_transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.hky85 import (
    hky85_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.jc69 import jc69_transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.k80 import k80_transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.posteriors import (
    compute_marginal_state_posteriors,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .brownian_continuous_trait import BrownianContinuousTraitRunReport
from .discrete_trait_mk import (
    DiscreteTraitMkModelDefinition,
    DiscreteTraitMkRunReport,
    _build_initial_transition_rate_rows,
    _build_rate_matrix_from_transition_rows,
    _resolve_root_prior_for_state,
    _select_most_likely_state,
    _validate_nonblank_state_name,
)
from .discrete_trait_rate_parameters import resolve_discrete_trait_rate_rows
from .fixed_topology_dna import FixedTopologyDnaRunReport
from .fixed_topology_partitioned_dna import FixedTopologyPartitionedDnaRunReport
from .joint_topology_dna import JointTopologyDnaRunReport
from .ornstein_uhlenbeck_continuous_trait import (
    OrnsteinUhlenbeckContinuousTraitRunReport,
)
from .partition_model_priors import PartitionSubstitutionModelDefinition
from .partition_model_state import (
    resolve_partition_parameter_linkage_plan_from_model_parameters,
    resolve_partition_parameter_states_from_model_parameters,
)
from .posterior_ancestral_sequences import (
    _validate_alignment_records,
    _validate_low_confidence_state_symbol,
    _validate_optional_partition_surface,
    _validate_probability_threshold,
    _validate_sampled_states,
)
from .posterior_ancestral_traits import (
    _build_continuous_mixture_draws,
    _build_tree_depth_index,
    _evaluate_brownian_covariance,
    _evaluate_ou_covariance,
)
from .posterior_sets.diagnostics import highest_posterior_density_interval
from .state import BayesianPhylogeneticState

_DNA_STATES = ("A", "C", "G", "T")


@dataclass(frozen=True, slots=True)
class PosteriorMissingNucleotideDefinition:
    """Validated posterior missing-state summary configuration for DNA alignments."""

    records: tuple[AlignmentRecord, ...]
    missing_state_symbols: tuple[str, ...]
    posterior_probability_threshold: float
    low_confidence_state_symbol: str
    locus_partitions: tuple[LocusPartition, ...] | None = None
    partition_models: tuple[PartitionSubstitutionModelDefinition, ...] | None = None


@dataclass(frozen=True, slots=True)
class PosteriorMissingNucleotideStateProbabilityRow:
    """One taxon-site-state posterior probability for one masked nucleotide."""

    taxon: str
    site_position: int
    state: str
    posterior_probability: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorMissingNucleotideSiteSummaryRow:
    """One taxon-site posterior summary across all DNA states."""

    taxon: str
    site_position: int
    observed_symbol: str
    consensus_state: str
    exported_state: str
    max_posterior_probability: float
    low_confidence: bool
    posterior_probability_a: float
    posterior_probability_c: float
    posterior_probability_g: float
    posterior_probability_t: float


@dataclass(frozen=True, slots=True)
class PosteriorMissingNucleotideSequenceRecord:
    """One imputed nucleotide sequence for one taxon with masked tip sites."""

    identifier: str
    masked_site_count: int
    sequence: str


@dataclass(frozen=True, slots=True)
class PosteriorMissingNucleotideReport:
    """Posterior missing-state summary for masked nucleotide observations."""

    sample_count: int
    site_count: int
    taxon_count: int
    masked_site_count: int
    distinct_topology_count: int
    sampled_substitution_models: list[str]
    state_probability_rows: list[PosteriorMissingNucleotideStateProbabilityRow]
    site_summary_rows: list[PosteriorMissingNucleotideSiteSummaryRow]
    sequence_records: list[PosteriorMissingNucleotideSequenceRecord]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class PosteriorMissingDiscreteTraitStateProbabilityRow:
    """One taxon-state posterior probability for one masked discrete trait."""

    taxon: str
    state: str
    posterior_probability: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorMissingDiscreteTraitTaxonSummaryRow:
    """One taxon-level posterior summary for one masked discrete trait."""

    taxon: str
    observed_symbol: str
    most_likely_state: str
    max_posterior_probability: float
    posterior_entropy: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorMissingDiscreteTraitReport:
    """Posterior missing-state summary for masked discrete traits."""

    sample_count: int
    distinct_topology_count: int
    sampled_transition_models: list[str]
    state_order: list[str]
    state_probability_rows: list[PosteriorMissingDiscreteTraitStateProbabilityRow]
    taxon_summary_rows: list[PosteriorMissingDiscreteTraitTaxonSummaryRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class PosteriorMissingContinuousTraitTaxonSummaryRow:
    """One taxon-level posterior summary for one masked continuous trait."""

    taxon: str
    posterior_mean: float
    posterior_median: float
    posterior_hpd_95_lower: float
    posterior_hpd_95_upper: float
    mean_conditional_standard_deviation: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorMissingContinuousTraitReport:
    """Posterior missing-value summary for masked continuous traits."""

    sample_count: int
    distinct_topology_count: int
    sampled_trait_models: list[str]
    taxon_summary_rows: list[PosteriorMissingContinuousTraitTaxonSummaryRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class PosteriorMissingDiscreteTraitDefinition:
    """Validated posterior missing-state configuration for discrete traits."""

    observed_tip_states: dict[str, str]
    missing_taxa: tuple[str, ...]
    masked_symbols_by_taxon: dict[str, str]
    missing_state_symbols: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PosteriorMissingContinuousTraitDefinition:
    """Validated posterior missing-value configuration for continuous traits."""

    observed_tip_values: dict[str, float]
    missing_taxa: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _PosteriorMissingContinuousTraitDistribution:
    taxon: str
    conditional_mean: float
    conditional_standard_deviation: float


def build_posterior_missing_nucleotide_definition(
    *,
    records: list[AlignmentRecord] | tuple[AlignmentRecord, ...],
    missing_state_symbols: tuple[str, ...] | list[str] = ("N", "?"),
    posterior_probability_threshold: float = 0.5,
    low_confidence_state_symbol: str = "N",
    locus_partitions: list[LocusPartition] | tuple[LocusPartition, ...] | None = None,
    partition_models: (
        list[PartitionSubstitutionModelDefinition]
        | tuple[PartitionSubstitutionModelDefinition, ...]
        | None
    ) = None,
) -> PosteriorMissingNucleotideDefinition:
    """Build one validated posterior missing-nucleotide summary definition."""
    validated_records = _validate_alignment_records(records)
    validated_missing_symbols = _validate_missing_state_symbols(missing_state_symbols)
    validated_locus_partitions, validated_partition_models = (
        _validate_optional_partition_surface(
            locus_partitions=locus_partitions,
            partition_models=partition_models,
            records=validated_records,
        )
    )
    return PosteriorMissingNucleotideDefinition(
        records=validated_records,
        missing_state_symbols=validated_missing_symbols,
        posterior_probability_threshold=_validate_probability_threshold(
            value=posterior_probability_threshold,
            field_name="posterior_probability_threshold",
            owner_name="posterior missing-nucleotide definition",
        ),
        low_confidence_state_symbol=_validate_low_confidence_state_symbol(
            low_confidence_state_symbol
        ),
        locus_partitions=validated_locus_partitions,
        partition_models=validated_partition_models,
    )


def build_posterior_missing_discrete_trait_definition(
    *,
    tip_states: Mapping[str, str],
    missing_state_symbols: tuple[str, ...] | list[str] = ("?", "NA"),
) -> PosteriorMissingDiscreteTraitDefinition:
    """Build one validated posterior missing-state definition for discrete traits."""
    if not isinstance(tip_states, Mapping):
        raise PhylogeneticsError(
            "posterior missing discrete-trait definition requires one mapping of tip states",
            code="posterior_missing_discrete_trait_tip_states_type_invalid",
        )
    validated_missing_symbols = _validate_missing_discrete_trait_state_symbols(
        missing_state_symbols
    )
    observed_tip_states: dict[str, str] = {}
    missing_taxa: list[str] = []
    masked_symbols_by_taxon: dict[str, str] = {}
    for raw_taxon, raw_state in tip_states.items():
        taxon = _validate_nonblank_state_name(
            raw_taxon,
            field_name=f"tip_states[{raw_taxon!r}]",
            owner_name="posterior missing discrete-trait definition",
        )
        state = _validate_nonblank_state_name(
            raw_state,
            field_name=f"tip_states[{raw_taxon!r}]",
            owner_name="posterior missing discrete-trait definition",
        )
        if state in validated_missing_symbols:
            missing_taxa.append(taxon)
            masked_symbols_by_taxon[taxon] = state
            continue
        observed_tip_states[taxon] = state
    if set(observed_tip_states.values()) & set(validated_missing_symbols):
        raise PhylogeneticsError(
            "posterior missing discrete-trait definition requires missing-state symbols distinct from observed discrete states",
            code="posterior_missing_discrete_trait_symbol_conflicts_with_state",
            details={
                "observed_states": sorted(set(observed_tip_states.values())),
                "missing_state_symbols": list(validated_missing_symbols),
            },
        )
    if not missing_taxa:
        raise PhylogeneticsError(
            "posterior missing discrete-trait definition requires at least one masked taxon",
            code="posterior_missing_discrete_trait_masked_taxa_missing",
        )
    return PosteriorMissingDiscreteTraitDefinition(
        observed_tip_states=observed_tip_states,
        missing_taxa=tuple(sorted(missing_taxa)),
        masked_symbols_by_taxon=masked_symbols_by_taxon,
        missing_state_symbols=validated_missing_symbols,
    )


def build_posterior_missing_continuous_trait_definition(
    *,
    tip_values: Mapping[str, float | None],
) -> PosteriorMissingContinuousTraitDefinition:
    """Build one validated posterior missing-value definition for continuous traits."""
    if not isinstance(tip_values, Mapping):
        raise PhylogeneticsError(
            "posterior missing continuous-trait definition requires one mapping of tip values",
            code="posterior_missing_continuous_trait_tip_values_type_invalid",
        )
    observed_tip_values: dict[str, float] = {}
    missing_taxa: list[str] = []
    for raw_taxon, raw_value in tip_values.items():
        taxon = _validate_nonblank_state_name(
            raw_taxon,
            field_name=f"tip_values[{raw_taxon!r}]",
            owner_name="posterior missing continuous-trait definition",
        )
        if raw_value is None:
            missing_taxa.append(taxon)
            continue
        value = float(raw_value)
        if math.isnan(value):
            missing_taxa.append(taxon)
            continue
        if not math.isfinite(value):
            raise PhylogeneticsError(
                "posterior missing continuous-trait definition requires finite observed tip values",
                code="posterior_missing_continuous_trait_tip_value_nonfinite",
                details={"taxon": taxon, "tip_value": raw_value},
            )
        observed_tip_values[taxon] = value
    if not missing_taxa:
        raise PhylogeneticsError(
            "posterior missing continuous-trait definition requires at least one masked taxon",
            code="posterior_missing_continuous_trait_masked_taxa_missing",
        )
    return PosteriorMissingContinuousTraitDefinition(
        observed_tip_values=observed_tip_values,
        missing_taxa=tuple(sorted(missing_taxa)),
    )


def summarize_fixed_topology_dna_posterior_missing_states(
    *,
    run_report: FixedTopologyDnaRunReport,
    definition: PosteriorMissingNucleotideDefinition,
) -> PosteriorMissingNucleotideReport:
    """Summarize masked nucleotide tip states across a fixed-topology DNA chain."""
    if not isinstance(run_report, FixedTopologyDnaRunReport):
        raise PhylogeneticsError(
            "posterior missing-nucleotide summary requires one FixedTopologyDnaRunReport",
            code="posterior_missing_nucleotide_fixed_topology_run_report_type_invalid",
        )
    return summarize_nucleotide_posterior_missing_states(
        sampled_states=run_report.chain_report.sampled_states,
        definition=definition,
    )


def summarize_joint_topology_dna_posterior_missing_states(
    *,
    run_report: JointTopologyDnaRunReport,
    definition: PosteriorMissingNucleotideDefinition,
) -> PosteriorMissingNucleotideReport:
    """Summarize masked nucleotide tip states across a joint-topology DNA chain."""
    if not isinstance(run_report, JointTopologyDnaRunReport):
        raise PhylogeneticsError(
            "posterior missing-nucleotide summary requires one JointTopologyDnaRunReport",
            code="posterior_missing_nucleotide_joint_topology_run_report_type_invalid",
        )
    return summarize_nucleotide_posterior_missing_states(
        sampled_states=run_report.chain_report.sampled_states,
        definition=definition,
    )


def summarize_fixed_topology_partitioned_dna_posterior_missing_states(
    *,
    run_report: FixedTopologyPartitionedDnaRunReport,
    definition: PosteriorMissingNucleotideDefinition,
) -> PosteriorMissingNucleotideReport:
    """Summarize masked nucleotide tip states across a partitioned DNA chain."""
    if not isinstance(run_report, FixedTopologyPartitionedDnaRunReport):
        raise PhylogeneticsError(
            "posterior missing-nucleotide summary requires one FixedTopologyPartitionedDnaRunReport",
            code="posterior_missing_nucleotide_partitioned_run_report_type_invalid",
        )
    return summarize_nucleotide_posterior_missing_states(
        sampled_states=run_report.chain_report.sampled_states,
        definition=definition,
    )


def summarize_nucleotide_posterior_missing_states(
    *,
    sampled_states: list[BayesianPhylogeneticState]
    | tuple[BayesianPhylogeneticState, ...],
    definition: PosteriorMissingNucleotideDefinition,
) -> PosteriorMissingNucleotideReport:
    """Summarize masked nucleotide observations across posterior sampled states."""
    if not isinstance(definition, PosteriorMissingNucleotideDefinition):
        raise PhylogeneticsError(
            "posterior missing-nucleotide summary requires one PosteriorMissingNucleotideDefinition",
            code="posterior_missing_nucleotide_definition_type_invalid",
        )
    validated_sampled_states = _validate_sampled_states(sampled_states)
    masked_symbols_by_taxon_site = _masked_symbols_by_taxon_site(
        records=definition.records,
        missing_state_symbols=definition.missing_state_symbols,
    )
    if not masked_symbols_by_taxon_site:
        raise PhylogeneticsError(
            "posterior missing-nucleotide summary requires at least one masked alignment site",
            code="posterior_missing_nucleotide_masked_site_missing",
        )
    posterior_sum_by_key: dict[tuple[str, int, str], float] = {}
    sampled_substitution_models: set[str] = set()
    topology_ids = {
        sampled_state.tree.topology_id for sampled_state in validated_sampled_states
    }
    for sampled_state in validated_sampled_states:
        sample_probabilities, model_label = (
            _evaluate_posterior_sample_missing_nucleotides(
                sampled_state=sampled_state,
                definition=definition,
                masked_symbols_by_taxon_site=masked_symbols_by_taxon_site,
            )
        )
        sampled_substitution_models.add(model_label)
        for taxon, site_position, state, posterior_probability in sample_probabilities:
            accumulation_key = (taxon, site_position, state)
            posterior_sum_by_key[accumulation_key] = float(
                format(
                    posterior_sum_by_key.get(accumulation_key, 0.0)
                    + posterior_probability,
                    ".15g",
                )
            )
    sample_count = len(validated_sampled_states)
    state_probability_rows = _build_missing_nucleotide_state_probability_rows(
        posterior_sum_by_key=posterior_sum_by_key,
        masked_symbols_by_taxon_site=masked_symbols_by_taxon_site,
        sample_count=sample_count,
    )
    site_summary_rows = _build_missing_nucleotide_site_summary_rows(
        state_probability_rows=state_probability_rows,
        masked_symbols_by_taxon_site=masked_symbols_by_taxon_site,
        posterior_probability_threshold=definition.posterior_probability_threshold,
        low_confidence_state_symbol=definition.low_confidence_state_symbol,
    )
    return PosteriorMissingNucleotideReport(
        sample_count=sample_count,
        site_count=len(definition.records[0].sequence),
        taxon_count=len(definition.records),
        masked_site_count=len(masked_symbols_by_taxon_site),
        distinct_topology_count=len(topology_ids),
        sampled_substitution_models=sorted(sampled_substitution_models),
        state_probability_rows=state_probability_rows,
        site_summary_rows=site_summary_rows,
        sequence_records=_build_missing_nucleotide_sequence_records(
            records=definition.records,
            site_summary_rows=site_summary_rows,
        ),
        warnings=_build_missing_nucleotide_warnings(
            distinct_topology_count=len(topology_ids),
        ),
    )


def summarize_discrete_trait_mk_posterior_missing_states(
    *,
    run_report: DiscreteTraitMkRunReport,
    definition: PosteriorMissingDiscreteTraitDefinition,
) -> PosteriorMissingDiscreteTraitReport:
    """Summarize masked discrete-trait tip states across one Mk posterior chain."""
    if not isinstance(run_report, DiscreteTraitMkRunReport):
        raise PhylogeneticsError(
            "posterior missing discrete-trait summary requires one DiscreteTraitMkRunReport",
            code="posterior_missing_discrete_trait_run_report_type_invalid",
        )
    if not isinstance(definition, PosteriorMissingDiscreteTraitDefinition):
        raise PhylogeneticsError(
            "posterior missing discrete-trait summary requires one PosteriorMissingDiscreteTraitDefinition",
            code="posterior_missing_discrete_trait_definition_type_invalid",
        )
    validated_sampled_states = _validate_sampled_states(
        run_report.chain_report.sampled_states
    )
    if set(run_report.taxa) != (
        set(definition.observed_tip_states) | set(definition.missing_taxa)
    ):
        raise PhylogeneticsError(
            "posterior missing discrete-trait summary requires tip states to match the sampled tree taxon set exactly",
            code="posterior_missing_discrete_trait_taxa_mismatch",
            details={
                "run_report_taxa": sorted(run_report.taxa),
                "definition_taxa": sorted(
                    set(definition.observed_tip_states) | set(definition.missing_taxa)
                ),
            },
        )
    state_counts = {
        state_name: sum(
            1
            for candidate_state in definition.observed_tip_states.values()
            if candidate_state == state_name
        )
        for state_name in run_report.state_order
    }
    topology_ids = {
        sampled_state.tree.topology_id for sampled_state in validated_sampled_states
    }
    probability_sum_by_key: dict[tuple[str, str], float] = {}
    for sampled_state in validated_sampled_states:
        rate_matrix = _build_discrete_trait_rate_matrix_for_missing_state(
            sampled_state=sampled_state,
            model_definition=run_report.model_definition,
            state_order=run_report.state_order,
        )
        root_prior = _resolve_root_prior_for_state(
            state_order=run_report.state_order,
            state_counts=state_counts,
            model_definition=run_report.model_definition,
        )
        posterior_by_node = (
            _estimate_marginal_discrete_state_probabilities_with_missing_tips(
                tree=sampled_state.tree.to_tree(),
                observed_tip_states=definition.observed_tip_states,
                missing_taxa=definition.missing_taxa,
                state_order=run_report.state_order,
                rate_matrix=rate_matrix,
                root_prior=root_prior,
            )
        )
        for taxon in definition.missing_taxa:
            for state in run_report.state_order:
                accumulation_key = (taxon, state)
                probability_sum_by_key[accumulation_key] = float(
                    format(
                        probability_sum_by_key.get(accumulation_key, 0.0)
                        + posterior_by_node[taxon][state],
                        ".15g",
                    )
                )
    state_probability_rows = _build_missing_discrete_trait_state_probability_rows(
        probability_sum_by_key=probability_sum_by_key,
        missing_taxa=definition.missing_taxa,
        state_order=run_report.state_order,
        sample_count=len(validated_sampled_states),
    )
    return PosteriorMissingDiscreteTraitReport(
        sample_count=len(validated_sampled_states),
        distinct_topology_count=len(topology_ids),
        sampled_transition_models=[run_report.model_definition.transition_model_name],
        state_order=list(run_report.state_order),
        state_probability_rows=state_probability_rows,
        taxon_summary_rows=_build_missing_discrete_trait_taxon_summary_rows(
            definition=definition,
            state_probability_rows=state_probability_rows,
            state_order=run_report.state_order,
        ),
        warnings=_build_posterior_missing_warnings(
            distinct_topology_count=len(topology_ids),
            summarized_surface="masked discrete-trait tip-state probabilities",
        ),
    )


def summarize_brownian_continuous_trait_posterior_missing_values(
    *,
    run_report: BrownianContinuousTraitRunReport,
    definition: PosteriorMissingContinuousTraitDefinition,
) -> PosteriorMissingContinuousTraitReport:
    """Summarize masked Brownian trait values across posterior sampled states."""
    if not isinstance(run_report, BrownianContinuousTraitRunReport):
        raise PhylogeneticsError(
            "posterior missing Brownian trait summary requires one BrownianContinuousTraitRunReport",
            code="posterior_missing_brownian_trait_run_report_type_invalid",
        )
    return _summarize_continuous_trait_posterior_missing_values(
        sampled_states=run_report.chain_report.sampled_states,
        taxa=run_report.taxa,
        definition=definition,
        sampled_trait_models=["brownian"],
        distribution_builder=_evaluate_brownian_missing_tip_distributions,
    )


def summarize_ornstein_uhlenbeck_continuous_trait_posterior_missing_values(
    *,
    run_report: OrnsteinUhlenbeckContinuousTraitRunReport,
    definition: PosteriorMissingContinuousTraitDefinition,
) -> PosteriorMissingContinuousTraitReport:
    """Summarize masked OU trait values across posterior sampled states."""
    if not isinstance(run_report, OrnsteinUhlenbeckContinuousTraitRunReport):
        raise PhylogeneticsError(
            "posterior missing OU trait summary requires one OrnsteinUhlenbeckContinuousTraitRunReport",
            code="posterior_missing_ou_trait_run_report_type_invalid",
        )
    return _summarize_continuous_trait_posterior_missing_values(
        sampled_states=run_report.chain_report.sampled_states,
        taxa=run_report.taxa,
        definition=definition,
        sampled_trait_models=["ornstein-uhlenbeck"],
        distribution_builder=_evaluate_ou_missing_tip_distributions,
    )


def summarize_continuous_trait_posterior_missing_values(
    *,
    run_report: BrownianContinuousTraitRunReport
    | OrnsteinUhlenbeckContinuousTraitRunReport,
    definition: PosteriorMissingContinuousTraitDefinition,
) -> PosteriorMissingContinuousTraitReport:
    """Summarize masked continuous trait values from one supported trait chain."""
    if isinstance(run_report, BrownianContinuousTraitRunReport):
        return summarize_brownian_continuous_trait_posterior_missing_values(
            run_report=run_report,
            definition=definition,
        )
    if isinstance(run_report, OrnsteinUhlenbeckContinuousTraitRunReport):
        return summarize_ornstein_uhlenbeck_continuous_trait_posterior_missing_values(
            run_report=run_report,
            definition=definition,
        )
    raise PhylogeneticsError(
        "posterior missing continuous-trait summary requires one BrownianContinuousTraitRunReport or OrnsteinUhlenbeckContinuousTraitRunReport",
        code="posterior_missing_continuous_trait_run_report_type_invalid",
    )


def _validate_missing_state_symbols(
    missing_state_symbols: tuple[str, ...] | list[str],
) -> tuple[str, ...]:
    validated_symbols = tuple(
        symbol.strip().upper() for symbol in missing_state_symbols
    )
    if not validated_symbols:
        raise PhylogeneticsError(
            "posterior missing-nucleotide definition requires at least one missing-state symbol",
            code="posterior_missing_nucleotide_symbols_empty",
        )
    if any(len(symbol) != 1 for symbol in validated_symbols):
        raise PhylogeneticsError(
            "posterior missing-nucleotide definition requires single-character missing-state symbols",
            code="posterior_missing_nucleotide_symbol_length_invalid",
            details={"missing_state_symbols": list(validated_symbols)},
        )
    if set(validated_symbols) & set(_DNA_STATES):
        raise PhylogeneticsError(
            "posterior missing-nucleotide definition requires missing-state symbols distinct from A, C, G, and T",
            code="posterior_missing_nucleotide_symbol_conflicts_with_dna_state",
            details={"missing_state_symbols": list(validated_symbols)},
        )
    if len(set(validated_symbols)) != len(validated_symbols):
        raise PhylogeneticsError(
            "posterior missing-nucleotide definition requires unique missing-state symbols",
            code="posterior_missing_nucleotide_symbols_duplicated",
            details={"missing_state_symbols": list(validated_symbols)},
        )
    return validated_symbols


def _masked_symbols_by_taxon_site(
    *,
    records: tuple[AlignmentRecord, ...],
    missing_state_symbols: tuple[str, ...],
) -> dict[tuple[str, int], str]:
    masked_symbols_by_taxon_site: dict[tuple[str, int], str] = {}
    valid_observed_symbols = set(DNA_STATE_ORDER) | set(missing_state_symbols)
    for record in records:
        normalized_sequence = record.sequence.upper()
        invalid_symbols = sorted(set(normalized_sequence) - valid_observed_symbols)
        if invalid_symbols:
            raise PhylogeneticsError(
                "posterior missing-nucleotide summary requires every alignment character to be one DNA state or one configured missing-state symbol",
                code="posterior_missing_nucleotide_alignment_symbol_invalid",
                details={
                    "taxon": record.identifier,
                    "invalid_symbols": invalid_symbols,
                    "allowed_dna_states": list(DNA_STATE_ORDER),
                    "missing_state_symbols": list(missing_state_symbols),
                },
            )
        for site_index, symbol in enumerate(normalized_sequence, start=1):
            if symbol in missing_state_symbols:
                masked_symbols_by_taxon_site[(record.identifier, site_index)] = symbol
    return masked_symbols_by_taxon_site


def _evaluate_posterior_sample_missing_nucleotides(
    *,
    sampled_state: BayesianPhylogeneticState,
    definition: PosteriorMissingNucleotideDefinition,
    masked_symbols_by_taxon_site: dict[tuple[str, int], str],
) -> tuple[list[tuple[str, int, str, float]], str]:
    if (
        definition.partition_models is not None
        and definition.locus_partitions is not None
    ):
        return _evaluate_partitioned_posterior_sample_missing_nucleotides(
            sampled_state=sampled_state,
            definition=definition,
            masked_symbols_by_taxon_site=masked_symbols_by_taxon_site,
        )
    return _evaluate_unpartitioned_posterior_sample_missing_nucleotides(
        sampled_state=sampled_state,
        definition=definition,
        masked_symbols_by_taxon_site=masked_symbols_by_taxon_site,
    )


def _evaluate_unpartitioned_posterior_sample_missing_nucleotides(
    *,
    sampled_state: BayesianPhylogeneticState,
    definition: PosteriorMissingNucleotideDefinition,
    masked_symbols_by_taxon_site: dict[tuple[str, int], str],
) -> tuple[list[tuple[str, int, str, float]], str]:
    model_name = sampled_state.model_parameters.categorical_parameters.get(
        "substitution-model"
    )
    if model_name in {None, "partitioned-dna"}:
        raise PhylogeneticsError(
            "posterior missing-nucleotide summary requires a named non-partitioned substitution model when no partition surface is configured",
            code="posterior_missing_nucleotide_model_name_missing",
            details={"observed_model_name": model_name},
        )
    probabilities = _evaluate_missing_nucleotide_probabilities_for_records(
        tree=sampled_state.tree.to_tree(),
        records=definition.records,
        missing_state_symbols=definition.missing_state_symbols,
        model_name=model_name,
        scalar_parameters=sampled_state.model_parameters.scalar_parameters,
        vector_parameters=sampled_state.model_parameters.vector_parameters,
        masked_symbols_by_taxon_site=masked_symbols_by_taxon_site,
    )
    return probabilities, model_name


def _evaluate_partitioned_posterior_sample_missing_nucleotides(
    *,
    sampled_state: BayesianPhylogeneticState,
    definition: PosteriorMissingNucleotideDefinition,
    masked_symbols_by_taxon_site: dict[tuple[str, int], str],
) -> tuple[list[tuple[str, int, str, float]], str]:
    model_name = sampled_state.model_parameters.categorical_parameters.get(
        "substitution-model"
    )
    if model_name != "partitioned-dna":
        raise PhylogeneticsError(
            "posterior missing-nucleotide summary with partition surfaces requires the partitioned DNA substitution-model label in every sampled state",
            code="posterior_missing_nucleotide_partitioned_model_name_invalid",
            details={"observed_model_name": model_name},
        )
    partition_models = require_present(
        definition.partition_models,
        owner_name="posterior missing-nucleotide partitioned evaluation",
        field_name="partition_models",
    )
    require_present(
        definition.locus_partitions,
        owner_name="posterior missing-nucleotide partitioned evaluation",
        field_name="locus_partitions",
    )
    partition_names = tuple(
        partition_model.partition_name for partition_model in partition_models
    )
    linkage_plan = resolve_partition_parameter_linkage_plan_from_model_parameters(
        model_parameters=sampled_state.model_parameters,
        partition_names=partition_names,
    )
    partition_states = resolve_partition_parameter_states_from_model_parameters(
        model_parameters=sampled_state.model_parameters,
        partition_models=partition_models,
        linkage_plan=linkage_plan,
    )
    partition_state_by_name = {
        partition_state.partition_name: partition_state
        for partition_state in partition_states
    }
    probabilities: list[tuple[str, int, str, float]] = []
    model_fragments: list[str] = []
    for locus_partition, partition_model in zip(
        definition.locus_partitions,
        definition.partition_models,
        strict=True,
    ):
        partition_state = partition_state_by_name.get(partition_model.partition_name)
        if partition_state is None:
            raise PhylogeneticsError(
                "posterior missing-nucleotide summary requires one realized parameter state for every configured partition",
                code="posterior_missing_nucleotide_partition_state_missing",
                details={"partition_name": partition_model.partition_name},
            )
        global_site_positions = _global_partition_site_positions(locus_partition)
        partition_masked_symbols_by_taxon_site = {
            (taxon, local_site_index): symbol
            for (taxon, site_position), symbol in masked_symbols_by_taxon_site.items()
            for local_site_index, global_site_position in enumerate(
                global_site_positions,
                start=1,
            )
            if global_site_position == site_position
        }
        if not partition_masked_symbols_by_taxon_site:
            model_fragments.append(
                f"{partition_model.partition_name}={partition_model.model_name}"
            )
            continue
        partition_records = tuple(
            AlignmentRecord(
                identifier=record.identifier,
                sequence=slice_partition_sequence(record.sequence, locus_partition),
            )
            for record in definition.records
        )
        partition_probabilities = (
            _evaluate_missing_nucleotide_probabilities_for_records(
                tree=sampled_state.tree.to_tree(),
                records=partition_records,
                missing_state_symbols=definition.missing_state_symbols,
                model_name=partition_model.base_model_name,
                scalar_parameters={"kappa": partition_state.kappa}
                if partition_state.kappa is not None
                else {},
                vector_parameters={
                    key: value
                    for key, value in (
                        ("base-frequencies", partition_state.base_frequencies),
                        ("exchangeabilities", partition_state.exchangeabilities),
                    )
                    if value is not None
                },
                masked_symbols_by_taxon_site=partition_masked_symbols_by_taxon_site,
            )
        )
        for (
            taxon,
            local_site_position,
            state,
            posterior_probability,
        ) in partition_probabilities:
            probabilities.append(
                (
                    taxon,
                    global_site_positions[local_site_position - 1],
                    state,
                    posterior_probability,
                )
            )
        model_fragments.append(
            f"{partition_model.partition_name}={partition_model.model_name}"
        )
    return probabilities, f"partitioned-dna[{','.join(model_fragments)}]"


def _evaluate_missing_nucleotide_probabilities_for_records(
    *,
    tree,
    records: tuple[AlignmentRecord, ...],
    missing_state_symbols: tuple[str, ...],
    model_name: str,
    scalar_parameters: dict[str, float],
    vector_parameters: dict[str, dict[str, float]],
    masked_symbols_by_taxon_site: dict[tuple[str, int], str],
) -> list[tuple[str, int, str, float]]:
    record_ids = [record.identifier for record in records]
    if set(record_ids) != set(tree.tip_names):
        raise PhylogeneticsError(
            "posterior missing-nucleotide summary requires the alignment taxon set to match every sampled tree",
            code="posterior_missing_nucleotide_taxa_mismatch",
            details={
                "record_ids": sorted(record_ids),
                "tree_taxa": sorted(tree.tip_names),
            },
        )
    root_prior, transition_by_child_id = _build_dna_posterior_surface(
        tree=tree,
        model_name=model_name,
        scalar_parameters=scalar_parameters,
        vector_parameters=vector_parameters,
    )
    tip_lookup = {
        node.name: node for node in tree.iter_leaves() if node.name is not None
    }
    symbols_by_taxon = {
        record.identifier: record.sequence.upper() for record in records
    }
    masked_taxa_by_site: dict[int, list[str]] = {}
    for taxon, site_position in masked_symbols_by_taxon_site:
        masked_taxa_by_site.setdefault(site_position, []).append(taxon)
    probabilities: list[tuple[str, int, str, float]] = []
    site_count = len(records[0].sequence)
    for site_position in range(1, site_count + 1):
        masked_taxa = masked_taxa_by_site.get(site_position)
        if not masked_taxa:
            continue
        posterior_pass = compute_marginal_state_posteriors(
            tree,
            state_count=len(DNA_STATE_ORDER),
            leaf_likelihood=lambda node, site_position=site_position: (
                _dna_leaf_likelihood_vector(
                    taxon=node.name,
                    symbols_by_taxon=symbols_by_taxon,
                    site_position=site_position,
                    missing_state_symbols=missing_state_symbols,
                )
            ),
            transition_matrix_for_child=lambda child: transition_by_child_id[
                child.node_id or ""
            ],
            root_prior=root_prior,
        )
        for taxon in sorted(masked_taxa):
            posterior = posterior_pass.posterior_for_node(tip_lookup[taxon])
            for state, posterior_probability in zip(
                DNA_STATE_ORDER,
                posterior,
                strict=True,
            ):
                probabilities.append(
                    (
                        taxon,
                        site_position,
                        state,
                        float(format(posterior_probability, ".15g")),
                    )
                )
    return probabilities


def _build_dna_posterior_surface(
    *,
    tree,
    model_name: str,
    scalar_parameters: dict[str, float],
    vector_parameters: dict[str, dict[str, float]],
) -> tuple[numpy.ndarray, dict[str, numpy.ndarray]]:
    if model_name == "JC69":
        root_prior = numpy.array(UNIFORM_DNA_ROOT_PRIOR, dtype=float)
        return root_prior, {
            child.node_id or "": jc69_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0)
            )
            for _parent, child in tree.iter_edges()
        }
    if model_name == "K80":
        root_prior = numpy.array(UNIFORM_DNA_ROOT_PRIOR, dtype=float)
        kappa = scalar_parameters["kappa"]
        return root_prior, {
            child.node_id or "": k80_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                kappa=kappa,
            )
            for _parent, child in tree.iter_edges()
        }
    if model_name == "F81":
        base_frequencies = vector_parameters["base-frequencies"]
        root_prior = numpy.array(
            [float(base_frequencies[state]) for state in DNA_STATE_ORDER],
            dtype=float,
        )
        return root_prior, {
            child.node_id or "": f81_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                base_frequencies=base_frequencies,
            )
            for _parent, child in tree.iter_edges()
        }
    if model_name == "HKY85":
        base_frequencies = vector_parameters["base-frequencies"]
        kappa = scalar_parameters["kappa"]
        root_prior = numpy.array(
            [float(base_frequencies[state]) for state in DNA_STATE_ORDER],
            dtype=float,
        )
        return root_prior, {
            child.node_id or "": hky85_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                base_frequencies=base_frequencies,
                kappa=kappa,
            )
            for _parent, child in tree.iter_edges()
        }
    if model_name == "GTR":
        base_frequencies = vector_parameters["base-frequencies"]
        exchangeabilities = vector_parameters["exchangeabilities"]
        root_prior = numpy.array(
            [float(base_frequencies[state]) for state in DNA_STATE_ORDER],
            dtype=float,
        )
        return root_prior, {
            child.node_id or "": gtr_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                exchangeabilities=exchangeabilities,
                base_frequencies=base_frequencies,
            )
            for _parent, child in tree.iter_edges()
        }
    raise PhylogeneticsError(
        "posterior missing-nucleotide summary encountered one unsupported substitution model",
        code="posterior_missing_nucleotide_model_unsupported",
        details={"model_name": model_name},
    )


def _dna_leaf_likelihood_vector(
    *,
    taxon: str | None,
    symbols_by_taxon: dict[str, str],
    site_position: int,
    missing_state_symbols: tuple[str, ...],
) -> numpy.ndarray:
    if taxon is None:
        raise PhylogeneticsError(
            "posterior missing-nucleotide summary requires every sampled tip to have one name",
            code="posterior_missing_nucleotide_tip_name_missing",
        )
    state_symbol = symbols_by_taxon[taxon][site_position - 1]
    if state_symbol in missing_state_symbols:
        return numpy.ones(len(DNA_STATE_ORDER), dtype=float)
    likelihood = numpy.zeros(len(DNA_STATE_ORDER), dtype=float)
    likelihood[DNA_STATE_ORDER.index(state_symbol)] = 1.0
    return likelihood


def _build_missing_nucleotide_state_probability_rows(
    *,
    posterior_sum_by_key: dict[tuple[str, int, str], float],
    masked_symbols_by_taxon_site: dict[tuple[str, int], str],
    sample_count: int,
) -> list[PosteriorMissingNucleotideStateProbabilityRow]:
    rows: list[PosteriorMissingNucleotideStateProbabilityRow] = []
    for taxon, site_position in sorted(masked_symbols_by_taxon_site):
        for state in DNA_STATE_ORDER:
            rows.append(
                PosteriorMissingNucleotideStateProbabilityRow(
                    taxon=taxon,
                    site_position=site_position,
                    state=state,
                    posterior_probability=float(
                        format(
                            posterior_sum_by_key[(taxon, site_position, state)]
                            / sample_count,
                            ".15g",
                        )
                    ),
                    supporting_sample_count=sample_count,
                    total_sample_count=sample_count,
                )
            )
    return rows


def _build_missing_nucleotide_site_summary_rows(
    *,
    state_probability_rows: list[PosteriorMissingNucleotideStateProbabilityRow],
    masked_symbols_by_taxon_site: dict[tuple[str, int], str],
    posterior_probability_threshold: float,
    low_confidence_state_symbol: str,
) -> list[PosteriorMissingNucleotideSiteSummaryRow]:
    probabilities_by_taxon_site = {
        (taxon, site_position): {
            row.state: row.posterior_probability
            for row in state_probability_rows
            if row.taxon == taxon and row.site_position == site_position
        }
        for taxon, site_position in masked_symbols_by_taxon_site
    }
    rows: list[PosteriorMissingNucleotideSiteSummaryRow] = []
    for taxon, site_position in sorted(masked_symbols_by_taxon_site):
        probabilities = probabilities_by_taxon_site[(taxon, site_position)]
        consensus_state = max(
            DNA_STATE_ORDER,
            key=lambda state: (probabilities[state], -DNA_STATE_ORDER.index(state)),
        )
        max_posterior_probability = probabilities[consensus_state]
        low_confidence = max_posterior_probability < posterior_probability_threshold
        rows.append(
            PosteriorMissingNucleotideSiteSummaryRow(
                taxon=taxon,
                site_position=site_position,
                observed_symbol=masked_symbols_by_taxon_site[(taxon, site_position)],
                consensus_state=consensus_state,
                exported_state=(
                    low_confidence_state_symbol if low_confidence else consensus_state
                ),
                max_posterior_probability=max_posterior_probability,
                low_confidence=low_confidence,
                posterior_probability_a=probabilities["A"],
                posterior_probability_c=probabilities["C"],
                posterior_probability_g=probabilities["G"],
                posterior_probability_t=probabilities["T"],
            )
        )
    return rows


def _build_missing_nucleotide_sequence_records(
    *,
    records: tuple[AlignmentRecord, ...],
    site_summary_rows: list[PosteriorMissingNucleotideSiteSummaryRow],
) -> list[PosteriorMissingNucleotideSequenceRecord]:
    summary_by_taxon_site = {
        (row.taxon, row.site_position): row for row in site_summary_rows
    }
    sequence_records: list[PosteriorMissingNucleotideSequenceRecord] = []
    for record in records:
        masked_site_positions = sorted(
            site_position
            for taxon, site_position in summary_by_taxon_site
            if taxon == record.identifier
        )
        if not masked_site_positions:
            continue
        symbols = list(record.sequence.upper())
        for site_position in masked_site_positions:
            symbols[site_position - 1] = summary_by_taxon_site[
                (record.identifier, site_position)
            ].exported_state
        sequence_records.append(
            PosteriorMissingNucleotideSequenceRecord(
                identifier=record.identifier,
                masked_site_count=len(masked_site_positions),
                sequence="".join(symbols),
            )
        )
    return sequence_records


def _build_missing_nucleotide_warnings(
    *,
    distinct_topology_count: int,
) -> list[str]:
    return _build_posterior_missing_warnings(
        distinct_topology_count=distinct_topology_count,
        summarized_surface="masked tip-state probabilities",
    )


def _build_posterior_missing_warnings(
    *,
    distinct_topology_count: int,
    summarized_surface: str,
) -> list[str]:
    if distinct_topology_count <= 1:
        return []
    return [
        f"{summarized_surface} were aggregated across posterior samples spanning multiple tree topologies"
    ]


def _global_partition_site_positions(
    locus_partition: LocusPartition,
) -> list[int]:
    return [
        site_position
        for segment in locus_partition.segments
        for site_position in iterate_partition_sites(segment)
    ]


def _validate_missing_discrete_trait_state_symbols(
    missing_state_symbols: tuple[str, ...] | list[str],
) -> tuple[str, ...]:
    validated_symbols = tuple(
        _validate_nonblank_state_name(
            symbol,
            field_name="missing_state_symbols",
            owner_name="posterior missing discrete-trait definition",
        )
        for symbol in missing_state_symbols
    )
    if not validated_symbols:
        raise PhylogeneticsError(
            "posterior missing discrete-trait definition requires at least one missing-state symbol",
            code="posterior_missing_discrete_trait_symbols_empty",
        )
    if len(set(validated_symbols)) != len(validated_symbols):
        raise PhylogeneticsError(
            "posterior missing discrete-trait definition requires unique missing-state symbols",
            code="posterior_missing_discrete_trait_symbols_duplicated",
            details={"missing_state_symbols": list(validated_symbols)},
        )
    return validated_symbols


def _build_discrete_trait_rate_matrix_for_missing_state(
    *,
    sampled_state: BayesianPhylogeneticState,
    model_definition: DiscreteTraitMkModelDefinition,
    state_order: list[str],
) -> numpy.ndarray:
    allowed_transition_pairs = {
        (left_index, right_index)
        for left_index in range(len(state_order))
        for right_index in range(len(state_order))
        if left_index != right_index
    }
    transition_rate_rows = _build_initial_transition_rate_rows(
        state_order=state_order,
        model_definition=model_definition,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    try:
        resolved_transition_rate_rows = (
            sampled_state.model_parameters.vector_parameters["discrete-trait-rates"]
        )
    except KeyError as error:
        raise PhylogeneticsError(
            "posterior missing discrete-trait summary requires one 'discrete-trait-rates' vector parameter in every sampled state",
            code="posterior_missing_discrete_trait_rate_parameters_missing",
        ) from error

    return _build_rate_matrix_from_transition_rows(
        state_order=state_order,
        transition_rate_rows=resolve_discrete_trait_rate_rows(
            model=model_definition.transition_model_name,
            transition_rate_rows=transition_rate_rows,
            parameter_values=resolved_transition_rate_rows,
        ),
    )


def _estimate_marginal_discrete_state_probabilities_with_missing_tips(
    *,
    tree,
    observed_tip_states: dict[str, str],
    missing_taxa: tuple[str, ...],
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
) -> dict[str, dict[str, float]]:
    state_index = {state: index for index, state in enumerate(state_order)}
    tip_taxa = set(tree.tip_names)
    expected_taxa = set(observed_tip_states) | set(missing_taxa)
    if tip_taxa != expected_taxa:
        raise PhylogeneticsError(
            "posterior missing discrete-trait summary requires the tip-state definition to match every sampled tree taxon set",
            code="posterior_missing_discrete_trait_sampled_tree_taxa_mismatch",
            details={
                "tree_taxa": sorted(tip_taxa),
                "definition_taxa": sorted(expected_taxa),
            },
        )
    transition_evaluator = build_transition_matrix_evaluator(rate_matrix)
    posterior_pass = compute_marginal_state_posteriors(
        tree,
        state_count=len(state_order),
        leaf_likelihood=lambda node: _discrete_trait_leaf_likelihood_vector(
            taxon=node.name,
            observed_tip_states=observed_tip_states,
            missing_taxa=missing_taxa,
            state_index=state_index,
            state_count=len(state_order),
        ),
        transition_matrix_for_child=lambda child: (
            transition_evaluator.transition_probability_matrix(
                discrete_trait_branch_length(child)
            )
        ),
        root_prior=root_prior,
    )
    tip_lookup = {
        node.name: node for node in tree.iter_leaves() if node.name is not None
    }
    return {
        taxon: {
            state: float(format(probability, ".15g"))
            for state, probability in zip(
                state_order,
                posterior_pass.posterior_for_node(tip_lookup[taxon]),
                strict=True,
            )
        }
        for taxon in missing_taxa
    }


def _discrete_trait_leaf_likelihood_vector(
    *,
    taxon: str | None,
    observed_tip_states: dict[str, str],
    missing_taxa: tuple[str, ...],
    state_index: dict[str, int],
    state_count: int,
) -> numpy.ndarray:
    if taxon is None:
        raise PhylogeneticsError(
            "posterior missing discrete-trait summary requires every sampled tip to have one name",
            code="posterior_missing_discrete_trait_tip_name_missing",
        )
    if taxon in missing_taxa:
        return numpy.ones(state_count, dtype=float)
    likelihood = numpy.zeros(state_count, dtype=float)
    likelihood[state_index[observed_tip_states[taxon]]] = 1.0
    return likelihood


def _build_missing_discrete_trait_state_probability_rows(
    *,
    probability_sum_by_key: dict[tuple[str, str], float],
    missing_taxa: tuple[str, ...],
    state_order: list[str],
    sample_count: int,
) -> list[PosteriorMissingDiscreteTraitStateProbabilityRow]:
    rows: list[PosteriorMissingDiscreteTraitStateProbabilityRow] = []
    for taxon in missing_taxa:
        for state in state_order:
            rows.append(
                PosteriorMissingDiscreteTraitStateProbabilityRow(
                    taxon=taxon,
                    state=state,
                    posterior_probability=float(
                        format(
                            probability_sum_by_key[(taxon, state)] / sample_count,
                            ".15g",
                        )
                    ),
                    supporting_sample_count=sample_count,
                    total_sample_count=sample_count,
                )
            )
    return rows


def _build_missing_discrete_trait_taxon_summary_rows(
    *,
    definition: PosteriorMissingDiscreteTraitDefinition,
    state_probability_rows: list[PosteriorMissingDiscreteTraitStateProbabilityRow],
    state_order: list[str],
) -> list[PosteriorMissingDiscreteTraitTaxonSummaryRow]:
    rows_by_taxon: dict[
        str, dict[str, PosteriorMissingDiscreteTraitStateProbabilityRow]
    ] = {}
    for row in state_probability_rows:
        rows_by_taxon.setdefault(row.taxon, {})[row.state] = row
    summary_rows: list[PosteriorMissingDiscreteTraitTaxonSummaryRow] = []
    for taxon in definition.missing_taxa:
        ordered_rows = [rows_by_taxon[taxon][state] for state in state_order]
        most_likely_state = _select_most_likely_state(
            {row.state: row.posterior_probability for row in ordered_rows},
            state_order=state_order,
        )
        posterior_entropy = -math.fsum(
            row.posterior_probability * math.log(row.posterior_probability)
            for row in ordered_rows
            if row.posterior_probability > 0.0
        )
        summary_rows.append(
            PosteriorMissingDiscreteTraitTaxonSummaryRow(
                taxon=taxon,
                observed_symbol=definition.masked_symbols_by_taxon[taxon],
                most_likely_state=most_likely_state,
                max_posterior_probability=max(
                    row.posterior_probability for row in ordered_rows
                ),
                posterior_entropy=float(format(posterior_entropy, ".15g")),
                supporting_sample_count=ordered_rows[0].supporting_sample_count,
                total_sample_count=ordered_rows[0].total_sample_count,
            )
        )
    return summary_rows


def _summarize_continuous_trait_posterior_missing_values(
    *,
    sampled_states: list[BayesianPhylogeneticState]
    | tuple[BayesianPhylogeneticState, ...],
    taxa: list[str],
    definition: PosteriorMissingContinuousTraitDefinition,
    sampled_trait_models: list[str],
    distribution_builder,
) -> PosteriorMissingContinuousTraitReport:
    if not isinstance(definition, PosteriorMissingContinuousTraitDefinition):
        raise PhylogeneticsError(
            "posterior missing continuous-trait summary requires one PosteriorMissingContinuousTraitDefinition",
            code="posterior_missing_continuous_trait_definition_type_invalid",
        )
    validated_sampled_states = _validate_sampled_states(sampled_states)
    if set(taxa) != set(definition.observed_tip_values) | set(definition.missing_taxa):
        raise PhylogeneticsError(
            "posterior missing continuous-trait summary requires tip values to match the sampled tree taxon set exactly",
            code="posterior_missing_continuous_trait_taxa_mismatch",
            details={
                "run_report_taxa": sorted(taxa),
                "definition_taxa": sorted(
                    set(definition.observed_tip_values) | set(definition.missing_taxa)
                ),
            },
        )
    topology_ids = {
        sampled_state.tree.topology_id for sampled_state in validated_sampled_states
    }
    distributions_by_taxon: dict[
        str, list[_PosteriorMissingContinuousTraitDistribution]
    ] = {taxon: [] for taxon in definition.missing_taxa}
    for sampled_state in validated_sampled_states:
        for distribution in distribution_builder(
            sampled_state=sampled_state,
            observed_tip_values=definition.observed_tip_values,
            missing_taxa=definition.missing_taxa,
        ):
            distributions_by_taxon[distribution.taxon].append(distribution)
    return PosteriorMissingContinuousTraitReport(
        sample_count=len(validated_sampled_states),
        distinct_topology_count=len(topology_ids),
        sampled_trait_models=sorted(set(sampled_trait_models)),
        taxon_summary_rows=_build_missing_continuous_trait_taxon_summary_rows(
            distributions_by_taxon=distributions_by_taxon,
            sample_count=len(validated_sampled_states),
        ),
        warnings=_build_posterior_missing_warnings(
            distinct_topology_count=len(topology_ids),
            summarized_surface="masked continuous-trait tip values",
        ),
    )


def _evaluate_brownian_missing_tip_distributions(
    *,
    sampled_state: BayesianPhylogeneticState,
    observed_tip_values: dict[str, float],
    missing_taxa: tuple[str, ...],
) -> list[_PosteriorMissingContinuousTraitDistribution]:
    return _evaluate_continuous_missing_tip_distributions(
        sampled_state=sampled_state,
        observed_tip_values=observed_tip_values,
        missing_taxa=missing_taxa,
        location=float(sampled_state.model_parameters.scalar_parameters["root-state"]),
        scale=float(sampled_state.model_parameters.scalar_parameters["sigma-squared"]),
        covariance_evaluator=_evaluate_brownian_covariance,
    )


def _evaluate_ou_missing_tip_distributions(
    *,
    sampled_state: BayesianPhylogeneticState,
    observed_tip_values: dict[str, float],
    missing_taxa: tuple[str, ...],
) -> list[_PosteriorMissingContinuousTraitDistribution]:
    alpha = float(sampled_state.model_parameters.scalar_parameters["alpha"])
    return _evaluate_continuous_missing_tip_distributions(
        sampled_state=sampled_state,
        observed_tip_values=observed_tip_values,
        missing_taxa=missing_taxa,
        location=float(sampled_state.model_parameters.scalar_parameters["optimum"]),
        scale=float(sampled_state.model_parameters.scalar_parameters["sigma-squared"]),
        covariance_evaluator=lambda left_node, right_node, depth_index: (
            _evaluate_ou_covariance(
                left_node=left_node,
                right_node=right_node,
                depth_index=depth_index,
                alpha=alpha,
            )
        ),
    )


def _evaluate_continuous_missing_tip_distributions(
    *,
    sampled_state: BayesianPhylogeneticState,
    observed_tip_values: dict[str, float],
    missing_taxa: tuple[str, ...],
    location: float,
    scale: float,
    covariance_evaluator,
) -> list[_PosteriorMissingContinuousTraitDistribution]:
    sampled_tree = sampled_state.tree.to_tree()
    if sampled_tree.rooted is not True:
        raise PhylogeneticsError(
            "posterior missing continuous-trait summary requires rooted sampled trees",
            code="posterior_missing_continuous_trait_tree_rooting_invalid",
        )
    expected_taxa = set(observed_tip_values) | set(missing_taxa)
    if set(sampled_tree.tip_names) != expected_taxa:
        raise PhylogeneticsError(
            "posterior missing continuous-trait summary requires each sampled tree to match the observed tip-value taxon set",
            code="posterior_missing_continuous_trait_sampled_tree_taxa_mismatch",
            details={
                "tree_taxa": sorted(sampled_tree.tip_names),
                "definition_taxa": sorted(expected_taxa),
            },
        )
    tip_lookup = {
        node.name: node for node in sampled_tree.iter_leaves() if node.name is not None
    }
    ordered_missing_taxa = list(missing_taxa)
    ordered_missing_nodes = [tip_lookup[taxon] for taxon in ordered_missing_taxa]
    ordered_observed_taxa = sorted(observed_tip_values)
    ordered_observed_nodes = [tip_lookup[taxon] for taxon in ordered_observed_taxa]
    depth_index = _build_tree_depth_index(sampled_tree)
    missing_covariance = stable_covariance(
        [
            [
                covariance_evaluator(left_node, right_node, depth_index)
                for right_node in ordered_missing_nodes
            ]
            for left_node in ordered_missing_nodes
        ]
    )
    if not ordered_observed_nodes:
        conditional_means = [location] * len(ordered_missing_nodes)
        conditional_covariance = missing_covariance
    else:
        observed_covariance = stable_covariance(
            [
                [
                    covariance_evaluator(left_node, right_node, depth_index)
                    for right_node in ordered_observed_nodes
                ]
                for left_node in ordered_observed_nodes
            ]
        )
        cross_covariance = [
            [
                covariance_evaluator(left_node, right_node, depth_index)
                for right_node in ordered_observed_nodes
            ]
            for left_node in ordered_missing_nodes
        ]
        inverse_observed_covariance = numpy.array(
            invert_matrix(observed_covariance),
            dtype=float,
        )
        observed_residuals = numpy.array(
            [observed_tip_values[taxon] - location for taxon in ordered_observed_taxa],
            dtype=float,
        )
        conditional_mean_offsets = numpy.array(cross_covariance, dtype=float) @ (
            inverse_observed_covariance @ observed_residuals
        )
        conditional_means = [
            float(format(location + offset, ".15g"))
            for offset in conditional_mean_offsets
        ]
        conditional_covariance = (
            numpy.array(missing_covariance, dtype=float)
            - numpy.array(cross_covariance, dtype=float)
            @ inverse_observed_covariance
            @ numpy.array(cross_covariance, dtype=float).T
        )
    return [
        _PosteriorMissingContinuousTraitDistribution(
            taxon=taxon,
            conditional_mean=float(format(conditional_means[index], ".15g")),
            conditional_standard_deviation=float(
                format(
                    math.sqrt(
                        max(float(conditional_covariance[index][index]), 0.0) * scale
                    ),
                    ".15g",
                )
            ),
        )
        for index, taxon in enumerate(ordered_missing_taxa)
    ]


def _build_missing_continuous_trait_taxon_summary_rows(
    *,
    distributions_by_taxon: dict[
        str, list[_PosteriorMissingContinuousTraitDistribution]
    ],
    sample_count: int,
) -> list[PosteriorMissingContinuousTraitTaxonSummaryRow]:
    summary_rows: list[PosteriorMissingContinuousTraitTaxonSummaryRow] = []
    for taxon in sorted(distributions_by_taxon):
        distributions = distributions_by_taxon[taxon]
        mixture_draws = _build_continuous_mixture_draws(distributions)
        hpd_95_lower, hpd_95_upper = highest_posterior_density_interval(mixture_draws)
        summary_rows.append(
            PosteriorMissingContinuousTraitTaxonSummaryRow(
                taxon=taxon,
                posterior_mean=float(
                    format(
                        mean(
                            distribution.conditional_mean
                            for distribution in distributions
                        ),
                        ".15g",
                    )
                ),
                posterior_median=float(format(median(mixture_draws), ".15g")),
                posterior_hpd_95_lower=float(format(hpd_95_lower, ".15g")),
                posterior_hpd_95_upper=float(format(hpd_95_upper, ".15g")),
                mean_conditional_standard_deviation=float(
                    format(
                        mean(
                            distribution.conditional_standard_deviation
                            for distribution in distributions
                        ),
                        ".15g",
                    )
                ),
                supporting_sample_count=sample_count,
                total_sample_count=sample_count,
            )
        )
    return summary_rows


__all__ = [
    "PosteriorMissingContinuousTraitDefinition",
    "PosteriorMissingContinuousTraitReport",
    "PosteriorMissingContinuousTraitTaxonSummaryRow",
    "PosteriorMissingDiscreteTraitDefinition",
    "PosteriorMissingDiscreteTraitReport",
    "PosteriorMissingDiscreteTraitStateProbabilityRow",
    "PosteriorMissingDiscreteTraitTaxonSummaryRow",
    "PosteriorMissingNucleotideDefinition",
    "PosteriorMissingNucleotideReport",
    "PosteriorMissingNucleotideSequenceRecord",
    "PosteriorMissingNucleotideSiteSummaryRow",
    "PosteriorMissingNucleotideStateProbabilityRow",
    "build_posterior_missing_continuous_trait_definition",
    "build_posterior_missing_discrete_trait_definition",
    "build_posterior_missing_nucleotide_definition",
    "summarize_brownian_continuous_trait_posterior_missing_values",
    "summarize_continuous_trait_posterior_missing_values",
    "summarize_discrete_trait_mk_posterior_missing_states",
    "summarize_fixed_topology_dna_posterior_missing_states",
    "summarize_fixed_topology_partitioned_dna_posterior_missing_states",
    "summarize_joint_topology_dna_posterior_missing_states",
    "summarize_nucleotide_posterior_missing_states",
    "summarize_ornstein_uhlenbeck_continuous_trait_posterior_missing_values",
]
