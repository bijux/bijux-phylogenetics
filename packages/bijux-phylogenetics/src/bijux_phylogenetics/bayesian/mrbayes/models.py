from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.bayesian.posterior_sets.burnin import (
    BurninSensitivityCladeShift,
    BurninSensitivityParameterShift,
)
from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet


@dataclass(slots=True)
class MrBayesPreparationReport:
    alignment_path: Path
    nexus_path: Path
    partition_path: Path | None
    taxon_count: int
    character_count: int
    inferred_alphabet: AlignmentAlphabet
    partition_count: int
    partition_names: list[str]
    partition_data_types: list[str]
    partition_warnings: list[str]
    model: str
    rates: str
    ngen: int
    nchains: int
    samplefreq: int
    printfreq: int
    burnin_fraction: float


@dataclass(slots=True)
class MrBayesTraceRow:
    generation: int
    values: dict[str, float]


@dataclass(slots=True)
class MrBayesTraceReport:
    path: Path
    row_count: int
    columns: list[str]
    rows: list[MrBayesTraceRow]


@dataclass(slots=True)
class EffectiveSampleSize:
    parameter: str
    sample_count: int
    effective_sample_size: float


@dataclass(slots=True)
class MrBayesESSReport:
    path: Path
    sample_count: int
    effective_sample_sizes: list[EffectiveSampleSize]


@dataclass(slots=True)
class MrBayesParameterSummary:
    parameter: str
    sample_count: int
    effective_sample_size: float
    mean: float
    median: float
    standard_deviation: float
    minimum: float
    maximum: float
    hpd_95_lower: float
    hpd_95_upper: float
    first_half_mean: float
    second_half_mean: float
    standardized_mean_shift: float


@dataclass(slots=True)
class MrBayesParameterDiagnosticsReport:
    path: Path
    burnin_fraction: float
    burnin_row_count: int
    kept_row_count: int
    first_kept_generation: int
    last_kept_generation: int
    parameter_summaries: list[MrBayesParameterSummary]


@dataclass(slots=True)
class MrBayesPosteriorDecompositionRow:
    generation: int
    log_posterior: float
    log_likelihood: float
    log_prior: float
    decomposition_delta: float
    decomposition_valid: bool


@dataclass(slots=True)
class MrBayesPosteriorDecompositionReport:
    path: Path
    burnin_fraction: float
    burnin_row_count: int
    kept_row_count: int
    first_kept_generation: int
    last_kept_generation: int
    posterior_term_source: str
    likelihood_term_source: str
    prior_term_source: str
    identity_tolerance: float
    verified: bool
    maximum_absolute_delta: float
    rows: list[MrBayesPosteriorDecompositionRow]


@dataclass(slots=True)
class MrBayesConvergenceReport:
    path: Path
    sample_count: int
    converged: bool
    ess_threshold: float
    mean_shift_threshold: float
    warnings: list[dict[str, object]]
    parameter_summaries: list[dict[str, object]]


@dataclass(slots=True)
class MrBayesPosteriorSummaryReport:
    source_path: Path
    filtered_tree_set_path: Path
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    rooted_topology_count: int
    shared_taxa: list[str]
    consensus_newick: str
    clade_frequency_count: int


@dataclass(slots=True)
class MrBayesBurninSensitivitySlice:
    burnin_fraction: float
    burnin_tree_count: int
    kept_tree_count: int
    rooted_topology_count: int
    clade_frequency_count: int
    consensus_newick: str
    kept_row_count: int | None
    first_kept_generation: int | None
    last_kept_generation: int | None
    lnl_mean: float | None
    tree_length_mean: float | None


@dataclass(slots=True)
class MrBayesBurninSensitivityReport:
    posterior_tree_path: Path
    trace_path: Path | None
    slices: list[MrBayesBurninSensitivitySlice]
    changed_consensus_count: int
    parameter_shifts: list[BurninSensitivityParameterShift]
    clade_shifts: list[BurninSensitivityCladeShift]
    unstable_parameter_count: int
    unstable_clade_count: int
    warnings: list[str]


@dataclass(slots=True)
class MrBayesPosteriorTreeSample:
    tree_name: str
    generation: int | None
    rooted: bool | None
    tip_names: list[str]
    newick: str


@dataclass(slots=True)
class MrBayesPosteriorTreeSetReport:
    path: Path
    tree_count: int
    rooted_tree_count: int
    sampled_generations: list[int]
    tip_names: list[str]
    trees: list[MrBayesPosteriorTreeSample]


@dataclass(slots=True)
class MrBayesMcmcRow:
    generation: int
    values: dict[str, float | None]


@dataclass(slots=True)
class MrBayesMcmcReport:
    path: Path
    row_count: int
    columns: list[str]
    comment_lines: list[str]
    rows: list[MrBayesMcmcRow]


@dataclass(slots=True)
class MrBayesConsensusTreeReport:
    path: Path
    tree_name: str
    rooted: bool | None
    tip_names: list[str]
    consensus_newick: str
    annotated_node_count: int
    minimum_posterior_probability: float | None
    maximum_posterior_probability: float | None
    minimum_posterior_probability_percent: float | None
    maximum_posterior_probability_percent: float | None
