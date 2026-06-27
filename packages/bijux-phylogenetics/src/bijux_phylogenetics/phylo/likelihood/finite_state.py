from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
import hashlib
import math
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
from bijux_phylogenetics.phylo.likelihood.models import (
    CompiledFiniteStateSiteLogLikelihoodRow,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    FiniteStatePruningPass,
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
    refresh_postorder_conditional_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.sites import (
    expanded_site_log_likelihood_rows_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

FiniteStateLeafVector = Callable[[str, str], numpy.ndarray]
FiniteStateTransitionMatrix = Callable[[float], numpy.ndarray]
FiniteStatePatternCompiler = Callable[
    [list[AlignmentRecord], Path | None],
    CompressedAlignmentSitePatterns,
]
FiniteStateAlignmentLengthResolver = Callable[[list[AlignmentRecord]], int]
FiniteStatePartitionRecordBuilder = Callable[
    [list[AlignmentRecord], LocusPartition],
    list[AlignmentRecord],
]
FiniteStatePartitionSourceSitePositionsResolver = Callable[
    [LocusPartition],
    tuple[int, ...],
]


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
    transition_rate_matrix: numpy.ndarray | None = None
    normalization_policy: str | None = None
    parameter_values: Mapping[str, float] = field(default_factory=dict)

    @property
    def state_count(self) -> int:
        return len(self.state_labels)


@dataclass(frozen=True, slots=True)
class FiniteStateTransitionModelState:
    """One stable finite-state transition identity for shared matrix caching."""

    model_name: str
    state_labels: tuple[str, ...]
    observation_policy: str
    normalization_policy: str | None
    root_prior_source: str
    root_prior_digest: str
    parameter_items: tuple[tuple[str, float], ...]
    transition_rate_matrix_digest: str


@dataclass(frozen=True, slots=True)
class FiniteStateTransitionCacheKey:
    """One shared finite-state transition cache key for a model state and branch."""

    model_state: FiniteStateTransitionModelState
    branch_length: float


@dataclass(frozen=True, slots=True)
class FiniteStateSiteEvaluationKey:
    """One stable cached site-evaluation identity inside a compiled context."""

    taxon_order: tuple[str, ...]
    states: tuple[str, ...]


@dataclass(slots=True)
class CachedFiniteStateSiteEvaluation:
    """One cached pruning pass and score for one compiled site pattern."""

    pruning_pass: FiniteStatePruningPass
    log_likelihood: float


_SHARED_FINITE_STATE_TRANSITION_CACHE: dict[
    FiniteStateTransitionCacheKey,
    numpy.ndarray,
] = {}


@dataclass(slots=True)
class CompiledFiniteStateLikelihoodContext:
    """One tree-bound finite-state likelihood context with cached branches."""

    tree: PhyloTree
    surface: FiniteStateLikelihoodSurface
    transition_matrix_by_node_id: dict[str, numpy.ndarray]
    transition_model_state: FiniteStateTransitionModelState | None = None
    transition_cache_key_by_node_id: dict[str, FiniteStateTransitionCacheKey] = field(
        default_factory=dict
    )
    node_by_id: dict[str, TreeNode] = field(default_factory=dict)
    ancestor_path_by_child_node_id: dict[str, tuple[str, ...]] = field(
        default_factory=dict
    )
    site_evaluation_by_key: dict[
        FiniteStateSiteEvaluationKey,
        CachedFiniteStateSiteEvaluation,
    ] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class CompiledFiniteStatePatternTerm:
    """One compiled pattern plus its stable source positions and effective weight."""

    pattern_id: str
    states: tuple[str, ...]
    pattern_site_count: int
    source_site_positions: tuple[int, ...]
    source_site_weights: tuple[float, ...]
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


def clear_shared_finite_state_transition_cache() -> None:
    """Drop every shared finite-state transition matrix cached across contexts."""
    _SHARED_FINITE_STATE_TRANSITION_CACHE.clear()


def shared_finite_state_transition_cache_entry_count() -> int:
    """Report how many shared finite-state transition matrices are cached."""
    return len(_SHARED_FINITE_STATE_TRANSITION_CACHE)


def compile_finite_state_likelihood_context(
    tree: PhyloTree,
    surface: FiniteStateLikelihoodSurface,
) -> CompiledFiniteStateLikelihoodContext:
    """Compile one rooted tree plus one finite-state model into cached branches."""
    transition_matrix_by_node_id: dict[str, numpy.ndarray] = {}
    transition_cache_key_by_node_id: dict[str, FiniteStateTransitionCacheKey] = {}
    node_by_id = _node_index(tree)
    ancestor_path_by_child_node_id: dict[str, tuple[str, ...]] = {}
    transition_model_state = _transition_model_state(surface)
    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        branch_length = max(float(child.branch_length or 0.0), 0.0)
        transition_cache_key = (
            None
            if transition_model_state is None
            else FiniteStateTransitionCacheKey(
                model_state=transition_model_state,
                branch_length=branch_length,
            )
        )
        transition_matrix_by_node_id[child.node_id] = _compiled_transition_matrix(
            surface,
            branch_length=branch_length,
            transition_cache_key=transition_cache_key,
        )
        if transition_cache_key is not None:
            transition_cache_key_by_node_id[child.node_id] = transition_cache_key
        ancestor_path_by_child_node_id[child.node_id] = _ancestor_path_from_child(child)
    return CompiledFiniteStateLikelihoodContext(
        tree=tree,
        surface=surface,
        transition_matrix_by_node_id=transition_matrix_by_node_id,
        transition_model_state=transition_model_state,
        transition_cache_key_by_node_id=transition_cache_key_by_node_id,
        node_by_id=node_by_id,
        ancestor_path_by_child_node_id=ancestor_path_by_child_node_id,
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
    partition_alignment_length: int | None = None,
    partition_record_builder: FiniteStatePartitionRecordBuilder | None = None,
    partition_source_site_positions_resolver: (
        FiniteStatePartitionSourceSitePositionsResolver | None
    ) = None,
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
            alignment_length=(
                alignment_length
                if partition_alignment_length is None
                else partition_alignment_length
            ),
        )
    partition_contexts = _compile_partition_contexts(
        records,
        alignment_length=alignment_length,
        pattern_compiler=pattern_compiler,
        partition_map=partition_map,
        partition_record_builder=partition_record_builder,
        partition_source_site_positions_resolver=(
            partition_source_site_positions_resolver
        ),
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
    site_evaluation_key = FiniteStateSiteEvaluationKey(
        taxon_order=tuple(taxon_order),
        states=states,
    )
    cached_site_evaluation = context.site_evaluation_by_key.get(site_evaluation_key)
    if cached_site_evaluation is not None:
        return cached_site_evaluation.log_likelihood
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
    log_likelihood = log_likelihood_from_root_prior(
        context.tree,
        pruning_pass,
        root_prior=context.surface.root_prior,
    )
    context.site_evaluation_by_key[site_evaluation_key] = CachedFiniteStateSiteEvaluation(
        pruning_pass=pruning_pass,
        log_likelihood=log_likelihood,
    )
    return log_likelihood


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


def evaluate_compiled_finite_state_site_log_likelihood_rows(
    analysis: CompiledFiniteStateLikelihoodAnalysis,
) -> tuple[list[CompiledFiniteStateSiteLogLikelihoodRow], float]:
    """Expand one compiled finite-state analysis back to weighted source-site rows."""
    rows: list[CompiledFiniteStateSiteLogLikelihoodRow] = []
    weighted_total = 0.0
    for partition_context in analysis.partition_contexts:
        for term in partition_context.compiled_pattern_terms:
            log_likelihood = evaluate_finite_state_site_log_likelihood(
                analysis.context,
                term.states,
                taxon_order=partition_context.taxon_order,
            )
            weighted_total += term.effective_weight * log_likelihood
            rows.extend(
                CompiledFiniteStateSiteLogLikelihoodRow(
                    partition_name=partition_context.partition_name,
                    partition_data_type=partition_context.partition_data_type,
                    pattern_id=term.pattern_id,
                    pattern_site_count=term.pattern_site_count,
                    source_site_position=source_site_position,
                    site_states=term.states,
                    site_weight=site_weight,
                    log_likelihood=log_likelihood,
                )
                for source_site_position, site_weight in zip(
                    term.source_site_positions,
                    term.source_site_weights,
                    strict=True,
                )
            )
    rows.sort(key=lambda row: (row.source_site_position, row.partition_name, row.pattern_id))
    return rows, weighted_total


def validate_compiled_finite_state_site_log_likelihood_reconstruction(
    rows: list[CompiledFiniteStateSiteLogLikelihoodRow],
    *,
    expected_total_log_likelihood: float,
    expected_site_count: int,
    expected_pattern_count: int,
    owner_name: str,
) -> None:
    """Reject compiled source-site rows that do not reconstruct the declared total."""
    if len(rows) != expected_site_count:
        raise ValueError(
            f"{owner_name} emitted {len(rows)} source-site rows for {expected_site_count} declared sites"
        )
    unique_positions = {row.source_site_position for row in rows}
    if len(unique_positions) != expected_site_count:
        raise ValueError(
            f"{owner_name} emitted duplicate or missing source-site positions"
        )
    unique_pattern_count = len(
        {(row.partition_name, row.pattern_id) for row in rows}
    )
    if unique_pattern_count != expected_pattern_count:
        raise ValueError(
            f"{owner_name} emitted {unique_pattern_count} partition-pattern rows for {expected_pattern_count} declared patterns"
        )
    partition_pattern_terms: dict[
        tuple[str, str],
        tuple[int, tuple[str, ...], float],
    ] = {}
    for row in rows:
        prior = partition_pattern_terms.get((row.partition_name, row.pattern_id))
        current = (row.pattern_site_count, row.site_states, row.log_likelihood)
        if prior is None:
            partition_pattern_terms[(row.partition_name, row.pattern_id)] = current
            continue
        if prior != current:
            raise ValueError(
                f"{owner_name} emitted inconsistent rows for partition '{row.partition_name}' pattern '{row.pattern_id}'"
            )
    reconstructed_total = math.fsum(
        row.site_weight * row.log_likelihood for row in rows
    )
    if not math.isclose(
        reconstructed_total,
        expected_total_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            f"{owner_name} source-site rows did not reconstruct the declared total likelihood"
        )


def update_compiled_finite_state_context_branch_length(
    context: CompiledFiniteStateLikelihoodContext,
    *,
    child_node_id: str,
    branch_length: float,
) -> None:
    """Update one compiled branch length and refresh only the affected cached path."""
    node = context.node_by_id.get(child_node_id)
    if node is None:
        raise KeyError(f"compiled context does not contain node_id '{child_node_id}'")
    if node.parent is None:
        raise ValueError("compiled context cannot update the root branch length")
    normalized_branch_length = _normalized_branch_length(branch_length)
    existing_branch_length = _normalized_branch_length(float(node.branch_length or 0.0))
    if existing_branch_length == normalized_branch_length:
        node.branch_length = normalized_branch_length
        return
    node.branch_length = normalized_branch_length
    transition_cache_key = (
        None
        if context.transition_model_state is None
        else FiniteStateTransitionCacheKey(
            model_state=context.transition_model_state,
            branch_length=normalized_branch_length,
        )
    )
    context.transition_matrix_by_node_id[child_node_id] = _compiled_transition_matrix(
        context.surface,
        branch_length=normalized_branch_length,
        transition_cache_key=transition_cache_key,
    )
    if transition_cache_key is None:
        context.transition_cache_key_by_node_id.pop(child_node_id, None)
    else:
        context.transition_cache_key_by_node_id[child_node_id] = transition_cache_key
    affected_node_ids = context.ancestor_path_by_child_node_id[child_node_id]
    for cached_site_evaluation in context.site_evaluation_by_key.values():
        refresh_postorder_conditional_likelihoods(
            context.tree,
            cached_site_evaluation.pruning_pass,
            affected_node_ids=affected_node_ids,
            transition_matrix_for_child=lambda child: _transition_matrix_for_child(
                context,
                child_node_id=child.node_id,
            ),
        )
        cached_site_evaluation.log_likelihood = log_likelihood_from_root_prior(
            context.tree,
            cached_site_evaluation.pruning_pass,
            root_prior=context.surface.root_prior,
        )


def update_compiled_finite_state_analysis_branch_length(
    analysis: CompiledFiniteStateLikelihoodAnalysis,
    *,
    child_node_id: str,
    branch_length: float,
) -> None:
    """Update one compiled-analysis branch length and refresh affected site caches."""
    update_compiled_finite_state_context_branch_length(
        analysis.context,
        child_node_id=child_node_id,
        branch_length=branch_length,
    )


def update_compiled_finite_state_context_topology(
    context: CompiledFiniteStateLikelihoodContext,
    *,
    tree: PhyloTree,
) -> None:
    """Update one compiled context to a new topology while reusing unchanged subtrees."""
    moved_tree = tree.copy().refresh()
    validate_explicit_branch_lengths(moved_tree, model_name=context.surface.model_name)
    _validate_topology_update_taxa(context.tree, moved_tree, surface=context.surface)
    old_tree = context.tree
    old_transition_by_branch_signature = _transition_by_branch_signature(context)
    transition_matrix_by_node_id: dict[str, numpy.ndarray] = {}
    transition_cache_key_by_node_id: dict[str, FiniteStateTransitionCacheKey] = {}
    for _parent, child in moved_tree.iter_edges():
        if child.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        branch_length = _normalized_branch_length(float(child.branch_length or 0.0))
        branch_signature = frozenset(child.descendant_taxa)
        reused_transition = old_transition_by_branch_signature.get(
            (branch_signature, branch_length)
        )
        if reused_transition is not None:
            transition_matrix_by_node_id[child.node_id] = reused_transition
        else:
            transition_cache_key = (
                None
                if context.transition_model_state is None
                else FiniteStateTransitionCacheKey(
                    model_state=context.transition_model_state,
                    branch_length=branch_length,
                )
            )
            transition_matrix_by_node_id[child.node_id] = _compiled_transition_matrix(
                context.surface,
                branch_length=branch_length,
                transition_cache_key=transition_cache_key,
            )
        if context.transition_model_state is not None:
            transition_cache_key_by_node_id[child.node_id] = FiniteStateTransitionCacheKey(
                model_state=context.transition_model_state,
                branch_length=branch_length,
            )
    context.tree = moved_tree
    context.node_by_id = _node_index(moved_tree)
    context.ancestor_path_by_child_node_id = {
        child.node_id: _ancestor_path_from_child(child)
        for _parent, child in moved_tree.iter_edges()
        if child.node_id is not None
    }
    context.transition_matrix_by_node_id = transition_matrix_by_node_id
    context.transition_cache_key_by_node_id = transition_cache_key_by_node_id
    for site_evaluation_key, cached_site_evaluation in context.site_evaluation_by_key.items():
        updated_pruning_pass = _updated_topology_pruning_pass(
            old_tree=old_tree,
            old_pruning_pass=cached_site_evaluation.pruning_pass,
            moved_tree=moved_tree,
            context=context,
            site_evaluation_key=site_evaluation_key,
        )
        cached_site_evaluation.pruning_pass = updated_pruning_pass
        cached_site_evaluation.log_likelihood = log_likelihood_from_root_prior(
            moved_tree,
            updated_pruning_pass,
            root_prior=context.surface.root_prior,
        )


def update_compiled_finite_state_analysis_topology(
    analysis: CompiledFiniteStateLikelihoodAnalysis,
    *,
    tree: PhyloTree,
) -> None:
    """Update one compiled analysis to a new topology while reusing unchanged subtrees."""
    update_compiled_finite_state_context_topology(
        analysis.context,
        tree=tree,
    )


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


def _compiled_transition_matrix(
    surface: FiniteStateLikelihoodSurface,
    *,
    branch_length: float,
    transition_cache_key: FiniteStateTransitionCacheKey | None,
) -> numpy.ndarray:
    if transition_cache_key is None:
        return surface.transition_matrix_for_branch_length(branch_length)
    cached_matrix = _SHARED_FINITE_STATE_TRANSITION_CACHE.get(transition_cache_key)
    if cached_matrix is not None:
        return cached_matrix
    resolved_matrix = numpy.asarray(
        surface.transition_matrix_for_branch_length(branch_length),
        dtype=float,
    )
    resolved_matrix.setflags(write=False)
    _SHARED_FINITE_STATE_TRANSITION_CACHE[transition_cache_key] = resolved_matrix
    return resolved_matrix


def _transition_model_state(
    surface: FiniteStateLikelihoodSurface,
) -> FiniteStateTransitionModelState | None:
    if surface.transition_rate_matrix is None:
        return None
    return FiniteStateTransitionModelState(
        model_name=surface.model_name,
        state_labels=surface.state_labels,
        observation_policy=surface.observation_policy,
        normalization_policy=surface.normalization_policy,
        root_prior_source=surface.root_prior_source,
        root_prior_digest=_float_array_digest(surface.root_prior),
        parameter_items=tuple(
            (key, float(value))
            for key, value in sorted(surface.parameter_values.items())
        ),
        transition_rate_matrix_digest=_float_array_digest(surface.transition_rate_matrix),
    )


def _float_array_digest(values: numpy.ndarray) -> str:
    normalized = numpy.asarray(values, dtype=float)
    return hashlib.sha256(normalized.tobytes(order="C")).hexdigest()


def _node_index(tree: PhyloTree) -> dict[str, TreeNode]:
    node_by_id: dict[str, TreeNode] = {}
    for node in tree.iter_nodes(order="preorder"):
        if node.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        node_by_id[node.node_id] = node
    return node_by_id


def _ancestor_path_from_child(child: TreeNode) -> tuple[str, ...]:
    path: list[str] = []
    current = child.parent
    while current is not None:
        if current.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        path.append(current.node_id)
        current = current.parent
    return tuple(path)


def _normalized_branch_length(branch_length: float) -> float:
    if not math.isfinite(branch_length):
        raise ValueError("branch length must be finite")
    return max(float(branch_length), 0.0)


def _validate_topology_update_taxa(
    original_tree: PhyloTree,
    moved_tree: PhyloTree,
    *,
    surface: FiniteStateLikelihoodSurface,
) -> None:
    original_taxa = sorted(original_tree.tip_names)
    moved_taxa = sorted(moved_tree.tip_names)
    if original_taxa != moved_taxa:
        raise ValueError(
            f"{surface.model_name} topology update requires identical tree taxa"
        )


def _transition_by_branch_signature(
    context: CompiledFiniteStateLikelihoodContext,
) -> dict[tuple[frozenset[str], float], numpy.ndarray]:
    transition_by_branch_signature: dict[
        tuple[frozenset[str], float],
        numpy.ndarray,
    ] = {}
    for _parent, child in context.tree.iter_edges():
        if child.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        transition_by_branch_signature[
            (
                frozenset(child.descendant_taxa),
                _normalized_branch_length(float(child.branch_length or 0.0)),
            )
        ] = context.transition_matrix_by_node_id[child.node_id]
    return transition_by_branch_signature


def _updated_topology_pruning_pass(
    *,
    old_tree: PhyloTree,
    old_pruning_pass: FiniteStatePruningPass,
    moved_tree: PhyloTree,
    context: CompiledFiniteStateLikelihoodContext,
    site_evaluation_key: FiniteStateSiteEvaluationKey,
) -> FiniteStatePruningPass:
    old_node_signature_map = _node_signature_by_node_id(old_tree)
    old_conditional_by_signature = {
        signature: old_pruning_pass.conditional_by_node_id[node_id]
        for node_id, signature in old_node_signature_map.items()
    }
    old_subtree_scale_by_signature = {
        signature: old_pruning_pass.subtree_log_scaling_by_node_id[node_id]
        for node_id, signature in old_node_signature_map.items()
    }
    leaf_conditional_by_name = {
        node.name: old_pruning_pass.conditional_by_node_id[node.node_id]
        for node in old_tree.iter_leaves()
        if node.name is not None and node.node_id is not None
    }
    conditional_by_node_id: dict[str, numpy.ndarray] = {}
    subtree_log_scaling_by_node_id: dict[str, float] = {}
    updated_pruning_pass = FiniteStatePruningPass(
        conditional_by_node_id=conditional_by_node_id,
        subtree_log_scaling_by_node_id=subtree_log_scaling_by_node_id,
        state_count=old_pruning_pass.state_count,
        normalized_internal_conditionals=old_pruning_pass.normalized_internal_conditionals,
    )
    states_by_taxon = dict(
        zip(site_evaluation_key.taxon_order, site_evaluation_key.states, strict=True)
    )
    for node in moved_tree.iter_nodes(order="postorder"):
        if node.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        if node.is_leaf():
            if node.name is None:
                raise ValueError(
                    f"{context.surface.model_name} topology update requires named tree tips"
                )
            conditional_by_node_id[node.node_id] = leaf_conditional_by_name.get(
                node.name,
                context.surface.leaf_likelihood_vector_for_state(
                    node.name,
                    states_by_taxon[node.name],
                ),
            )
            subtree_log_scaling_by_node_id[node.node_id] = 0.0
            continue
        node_signature = _node_topology_signature(node)
        reusable_conditional = old_conditional_by_signature.get(node_signature)
        reusable_scale = old_subtree_scale_by_signature.get(node_signature)
        if reusable_conditional is not None and reusable_scale is not None:
            conditional_by_node_id[node.node_id] = reusable_conditional
            subtree_log_scaling_by_node_id[node.node_id] = reusable_scale
            continue
        refresh_postorder_conditional_likelihoods(
            moved_tree,
            updated_pruning_pass,
            affected_node_ids=(node.node_id,),
            transition_matrix_for_child=lambda child: _transition_matrix_for_child(
                context,
                child_node_id=child.node_id,
            ),
        )
    return updated_pruning_pass


def _node_signature_by_node_id(tree: PhyloTree) -> dict[str, tuple[object, ...]]:
    signature_by_node_id: dict[str, tuple[object, ...]] = {}
    for node in tree.iter_nodes(order="postorder"):
        if node.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        signature_by_node_id[node.node_id] = _node_topology_signature(
            node,
            child_signatures=signature_by_node_id,
        )
    return signature_by_node_id


def _node_topology_signature(
    node: TreeNode,
    *,
    child_signatures: dict[str, tuple[object, ...]] | None = None,
) -> tuple[object, ...]:
    if node.is_leaf():
        return ("leaf", node.name)
    child_components: list[tuple[float, tuple[object, ...]]] = []
    for child in node.children:
        if child.node_id is None:
            raise ValueError("tree child node is missing a stable node_id")
        if child_signatures is None:
            child_signature = _node_topology_signature(child)
        else:
            child_signature = child_signatures[child.node_id]
        child_components.append(
            (
                _normalized_branch_length(float(child.branch_length or 0.0)),
                child_signature,
            )
        )
    return (
        "internal",
        tuple(sorted(child_components, key=repr)),
    )


def _compile_partition_contexts(
    records: list[AlignmentRecord],
    *,
    alignment_length: int,
    pattern_compiler: FiniteStatePatternCompiler,
    partition_map: tuple[LocusPartition, ...] | None,
    partition_record_builder: FiniteStatePartitionRecordBuilder | None,
    partition_source_site_positions_resolver: (
        FiniteStatePartitionSourceSitePositionsResolver | None
    ),
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
    resolved_partition_record_builder = (
        _default_partition_record_builder
        if partition_record_builder is None
        else partition_record_builder
    )
    resolved_partition_source_site_positions_resolver = (
        _default_partition_source_site_positions
        if partition_source_site_positions_resolver is None
        else partition_source_site_positions_resolver
    )
    return tuple(
        _compile_partition_context(
            partition_name=partition.name,
            partition_data_type=partition.data_type,
            records=resolved_partition_record_builder(
                records,
                partition,
            ),
            source_site_positions=resolved_partition_source_site_positions_resolver(
                partition
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
    resolved_source_site_weights = tuple(
        site_weight_by_position[site_position]
        for site_position in resolved_source_positions
    )
    return CompiledFiniteStatePatternTerm(
        pattern_id=pattern_id,
        states=states,
        pattern_site_count=pattern_site_count,
        source_site_positions=resolved_source_positions,
        source_site_weights=resolved_source_site_weights,
        effective_weight=sum(resolved_source_site_weights),
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


def _default_partition_record_builder(
    records: list[AlignmentRecord],
    partition: LocusPartition,
) -> list[AlignmentRecord]:
    return [
        AlignmentRecord(
            identifier=record.identifier,
            sequence=slice_partition_sequence(record.sequence, partition),
        )
        for record in records
    ]


def _default_partition_source_site_positions(
    partition: LocusPartition,
) -> tuple[int, ...]:
    return tuple(
        site
        for segment in partition.segments
        for site in iterate_partition_sites(segment)
    )


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
