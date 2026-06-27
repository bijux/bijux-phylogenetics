from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    iterate_partition_sites,
    slice_partition_sequence,
    validate_locus_partitions,
)
from bijux_phylogenetics.phylo.likelihood.models import SiteLogLikelihoodRow
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.sites import (
    expanded_site_log_likelihood_rows_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

FiniteStateLeafVector = Callable[[str, str], numpy.ndarray]
FiniteStateTransitionMatrix = Callable[[float], numpy.ndarray]
FiniteStatePatternCompiler = Callable[
    [list[AlignmentRecord], Path | None],
    CompressedAlignmentSitePatterns,
]
FiniteStateAlignmentLengthResolver = Callable[[list[AlignmentRecord]], int]


@dataclass(frozen=True, slots=True)
class FiniteStateLikelihoodSurface:
    """One reusable finite-state CTMC model surface over rooted-tree tips."""

    model_name: str
    state_labels: tuple[str, ...]
    root_prior: numpy.ndarray
    root_prior_source: str
    observation_policy: str
    transition_matrix_for_branch_length: FiniteStateTransitionMatrix
    leaf_likelihood_vector_for_state: FiniteStateLeafVector
    normalization_policy: str | None = None
    parameter_values: Mapping[str, float] = field(default_factory=dict)

    @property
    def state_count(self) -> int:
        return len(self.state_labels)


@dataclass(slots=True)
class CompiledFiniteStateLikelihoodContext:
    """One tree-bound finite-state likelihood context with cached branches."""

    tree: PhyloTree
    surface: FiniteStateLikelihoodSurface
    transition_matrix_by_node_id: dict[str, numpy.ndarray]


@dataclass(frozen=True, slots=True)
class CompiledFiniteStatePatternTerm:
    """One compiled pattern plus its stable source positions and effective weight."""

    pattern_id: str
    states: tuple[str, ...]
    pattern_site_count: int
    source_site_positions: tuple[int, ...]
    effective_weight: float


@dataclass(frozen=True, slots=True)
class CompiledFiniteStatePartitionContext:
    """One compiled partition block inside a reusable finite-state analysis."""

    partition_name: str
    partition_data_type: str | None
    alignment_length: int
    pattern_count: int
    taxon_order: list[str]
    compiled_pattern_terms: tuple[CompiledFiniteStatePatternTerm, ...]


@dataclass(frozen=True, slots=True)
class CompiledFiniteStateLikelihoodAnalysis:
    """One reusable tree plus alignment compilation for repeated likelihood work."""

    source_path: Path | None
    alignment_length: int
    pattern_count: int
    partition_count: int
    site_weight_total: float
    missing_data_policy: str
    context: CompiledFiniteStateLikelihoodContext
    partition_contexts: tuple[CompiledFiniteStatePartitionContext, ...]


def compile_finite_state_likelihood_context(
    tree: PhyloTree,
    surface: FiniteStateLikelihoodSurface,
) -> CompiledFiniteStateLikelihoodContext:
    """Compile one rooted tree plus one finite-state model into cached branches."""
    transition_matrix_by_node_id: dict[str, numpy.ndarray] = {}
    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        transition_matrix_by_node_id[child.node_id] = (
            surface.transition_matrix_for_branch_length(
                max(float(child.branch_length or 0.0), 0.0)
            )
        )
    return CompiledFiniteStateLikelihoodContext(
        tree=tree,
        surface=surface,
        transition_matrix_by_node_id=transition_matrix_by_node_id,
    )


def compile_finite_state_likelihood_analysis(
    tree: PhyloTree,
    surface: FiniteStateLikelihoodSurface,
    *,
    records: list[AlignmentRecord],
    pattern_compiler: FiniteStatePatternCompiler,
    missing_data_policy: str,
    alignment_length_resolver: FiniteStateAlignmentLengthResolver | None = None,
    partition_map: tuple[LocusPartition, ...] | None = None,
    site_weight_by_position: Mapping[int, float] | None = None,
    source_path: Path | None = None,
) -> CompiledFiniteStateLikelihoodAnalysis:
    """Compile one rooted tree, model surface, and alignment into reusable state."""
    validate_explicit_branch_lengths(tree, model_name=surface.model_name)
    resolved_alignment_length = (
        _validated_alignment_length(records)
        if alignment_length_resolver is None
        else alignment_length_resolver(records)
    )
    alignment_length = _validated_compiled_alignment_length(
        resolved_alignment_length,
        records=records,
    )
    validated_site_weight_by_position = _validated_site_weight_by_position(
        site_weight_by_position,
        alignment_length=alignment_length,
    )
    if partition_map is not None:
        validate_locus_partitions(
            partition_map,
            alignment_length=alignment_length,
        )
    partition_contexts = _compile_partition_contexts(
        records,
        alignment_length=alignment_length,
        pattern_compiler=pattern_compiler,
        partition_map=partition_map,
        site_weight_by_position=validated_site_weight_by_position,
        source_path=source_path,
    )
    context = compile_finite_state_likelihood_context(tree, surface)
    _validate_partition_taxa_against_tree(tree, partition_contexts, surface=surface)
    return CompiledFiniteStateLikelihoodAnalysis(
        source_path=source_path,
        alignment_length=alignment_length,
        pattern_count=sum(
            partition_context.pattern_count
            for partition_context in partition_contexts
        ),
        partition_count=len(partition_contexts),
        site_weight_total=sum(
            term.effective_weight
            for partition_context in partition_contexts
            for term in partition_context.compiled_pattern_terms
        ),
        missing_data_policy=missing_data_policy,
        context=context,
        partition_contexts=partition_contexts,
    )


def evaluate_finite_state_site_log_likelihood(
    context: CompiledFiniteStateLikelihoodContext,
    states: tuple[str, ...],
    *,
    taxon_order: list[str],
) -> float:
    """Evaluate one aligned site against one compiled finite-state context."""
    states_by_taxon = dict(zip(taxon_order, states, strict=True))
    pruning_pass = postorder_conditional_likelihoods(
        context.tree,
        state_count=context.surface.state_count,
        leaf_likelihood=lambda node: _leaf_likelihood_vector(
            context.surface,
            states_by_taxon=states_by_taxon,
            node_name=node.name,
        ),
        transition_matrix_for_child=lambda child: _transition_matrix_for_child(
            context,
            child_node_id=child.node_id,
        ),
    )
    return log_likelihood_from_root_prior(
        context.tree,
        pruning_pass,
        root_prior=context.surface.root_prior,
    )


def evaluate_finite_state_tree_log_likelihood_from_patterns(
    context: CompiledFiniteStateLikelihoodContext,
    compressed_patterns: CompressedAlignmentSitePatterns,
) -> float:
    """Evaluate one compressed alignment against one compiled finite-state model."""
    total = 0.0
    for pattern in compressed_patterns.patterns:
        total += pattern.weight * evaluate_finite_state_site_log_likelihood(
            context,
            pattern.states,
            taxon_order=compressed_patterns.taxon_order,
        )
    return total


def evaluate_finite_state_site_log_likelihood_rows_from_patterns(
    context: CompiledFiniteStateLikelihoodContext,
    compressed_patterns: CompressedAlignmentSitePatterns,
) -> tuple[list[SiteLogLikelihoodRow], float]:
    """Expand one finite-state pattern set back to per-site rows and total likelihood."""
    return expanded_site_log_likelihood_rows_from_patterns(
        compressed_patterns,
        site_log_likelihood=lambda states: evaluate_finite_state_site_log_likelihood(
            context,
            states,
            taxon_order=compressed_patterns.taxon_order,
        ),
    )


def evaluate_finite_state_tree_log_likelihood(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    surface: FiniteStateLikelihoodSurface,
) -> float:
    """Evaluate one compressed alignment on one shared finite-state surface."""
    validate_explicit_branch_lengths(tree, model_name=surface.model_name)
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name=surface.model_name,
    )
    context = compile_finite_state_likelihood_context(tree, surface)
    return evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )


def evaluate_compiled_finite_state_likelihood_analysis(
    analysis: CompiledFiniteStateLikelihoodAnalysis,
) -> float:
    """Evaluate one reusable compiled finite-state analysis without reparsing inputs."""
    total = 0.0
    for partition_context in analysis.partition_contexts:
        for term in partition_context.compiled_pattern_terms:
            if term.effective_weight == 0.0:
                continue
            total += term.effective_weight * evaluate_finite_state_site_log_likelihood(
                analysis.context,
                term.states,
                taxon_order=partition_context.taxon_order,
            )
    return total


def _leaf_likelihood_vector(
    surface: FiniteStateLikelihoodSurface,
    *,
    states_by_taxon: dict[str, str],
    node_name: str | None,
) -> numpy.ndarray:
    if node_name is None:
        raise ValueError(
            f"{surface.model_name} likelihood requires named tree tips for alignment lookup"
        )
    return surface.leaf_likelihood_vector_for_state(
        node_name,
        states_by_taxon[node_name],
    )


def _transition_matrix_for_child(
    context: CompiledFiniteStateLikelihoodContext,
    *,
    child_node_id: str | None,
) -> numpy.ndarray:
    if child_node_id is None:
        raise ValueError("tree node is missing a stable node_id")
    return context.transition_matrix_by_node_id[child_node_id]


def _compile_partition_contexts(
    records: list[AlignmentRecord],
    *,
    alignment_length: int,
    pattern_compiler: FiniteStatePatternCompiler,
    partition_map: tuple[LocusPartition, ...] | None,
    site_weight_by_position: Mapping[int, float],
    source_path: Path | None,
) -> tuple[CompiledFiniteStatePartitionContext, ...]:
    if partition_map is None:
        return (
            _compile_partition_context(
                partition_name="full-alignment",
                partition_data_type=None,
                records=records,
                source_site_positions=tuple(range(1, alignment_length + 1)),
                pattern_compiler=pattern_compiler,
                site_weight_by_position=site_weight_by_position,
                source_path=source_path,
            ),
        )
    return tuple(
        _compile_partition_context(
            partition_name=partition.name,
            partition_data_type=partition.data_type,
            records=[
                AlignmentRecord(
                    identifier=record.identifier,
                    sequence=slice_partition_sequence(record.sequence, partition),
                )
                for record in records
            ],
            source_site_positions=tuple(
                site
                for segment in partition.segments
                for site in iterate_partition_sites(segment)
            ),
            pattern_compiler=pattern_compiler,
            site_weight_by_position=site_weight_by_position,
            source_path=None,
        )
        for partition in partition_map
    )


def _compile_partition_context(
    *,
    partition_name: str,
    partition_data_type: str | None,
    records: list[AlignmentRecord],
    source_site_positions: tuple[int, ...],
    pattern_compiler: FiniteStatePatternCompiler,
    site_weight_by_position: Mapping[int, float],
    source_path: Path | None,
) -> CompiledFiniteStatePartitionContext:
    compressed_patterns = pattern_compiler(records, source_path)
    compiled_terms = tuple(
        _compile_pattern_term(
            pattern_id=pattern.pattern_id,
            states=pattern.states,
            pattern_site_count=pattern.weight,
            local_site_positions=pattern.site_positions,
            source_site_positions=source_site_positions,
            site_weight_by_position=site_weight_by_position,
        )
        for pattern in compressed_patterns.patterns
    )
    return CompiledFiniteStatePartitionContext(
        partition_name=partition_name,
        partition_data_type=partition_data_type,
        alignment_length=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        taxon_order=compressed_patterns.taxon_order,
        compiled_pattern_terms=compiled_terms,
    )


def _compile_pattern_term(
    *,
    pattern_id: str,
    states: tuple[str, ...],
    pattern_site_count: int,
    local_site_positions: list[int],
    source_site_positions: tuple[int, ...],
    site_weight_by_position: Mapping[int, float],
) -> CompiledFiniteStatePatternTerm:
    resolved_source_positions = tuple(
        source_site_positions[local_site_position - 1]
        for local_site_position in local_site_positions
    )
    return CompiledFiniteStatePatternTerm(
        pattern_id=pattern_id,
        states=states,
        pattern_site_count=pattern_site_count,
        source_site_positions=resolved_source_positions,
        effective_weight=sum(
            site_weight_by_position[site_position]
            for site_position in resolved_source_positions
        ),
    )


def _validated_alignment_length(records: list[AlignmentRecord]) -> int:
    if not records:
        raise ValueError(
            "compiled finite-state likelihood analysis requires records"
        )
    alignment_lengths = {len(record.sequence) for record in records}
    if len(alignment_lengths) != 1:
        raise ValueError(
            "compiled finite-state likelihood analysis requires equal-length aligned records"
        )
    return next(iter(alignment_lengths))


def _validated_compiled_alignment_length(
    alignment_length: int,
    *,
    records: list[AlignmentRecord],
) -> int:
    if not records:
        raise ValueError(
            "compiled finite-state likelihood analysis requires records"
        )
    if alignment_length < 1:
        raise ValueError(
            "compiled finite-state likelihood analysis requires a positive alignment length"
        )
    return alignment_length


def _validated_site_weight_by_position(
    site_weight_by_position: Mapping[int, float] | None,
    *,
    alignment_length: int,
) -> dict[int, float]:
    validated = {position: 1.0 for position in range(1, alignment_length + 1)}
    if site_weight_by_position is None:
        return validated
    for position, raw_weight in site_weight_by_position.items():
        if position < 1 or position > alignment_length:
            raise ValueError(
                "compiled finite-state likelihood analysis received a site weight outside the alignment"
            )
        weight = float(raw_weight)
        if not numpy.isfinite(weight) or weight < 0.0:
            raise ValueError(
                "compiled finite-state likelihood analysis requires finite nonnegative site weights"
            )
        validated[position] = weight
    return validated


def _validate_partition_taxa_against_tree(
    tree: PhyloTree,
    partition_contexts: tuple[CompiledFiniteStatePartitionContext, ...],
    *,
    surface: FiniteStateLikelihoodSurface,
) -> None:
    for partition_context in partition_contexts:
        _validate_tree_taxa_against_order(
            tree,
            taxon_order=partition_context.taxon_order,
            model_name=surface.model_name,
        )


def _validate_tree_taxa_against_order(
    tree: PhyloTree,
    *,
    taxon_order: list[str],
    model_name: str,
) -> None:
    validate_tree_taxa_against_patterns(
        tree,
        CompressedAlignmentSitePatterns(
            source_path=None,
            taxon_order=taxon_order,
            alignment_length=1,
            pattern_count=1,
            patterns=[],
        ),
        model_name=model_name,
    )
