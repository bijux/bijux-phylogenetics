from __future__ import annotations

from collections.abc import Callable

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.sites import (
    sum_compressed_site_pattern_log_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    AlignmentTaxonMismatchError,
    InvalidAlignmentError,
)

PROTEIN_STATE_ORDER = (
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
)
PROTEIN_STATE_INDEX = {state: index for index, state in enumerate(PROTEIN_STATE_ORDER)}
UNIFORM_PROTEIN_ROOT_PRIOR = numpy.full(
    len(PROTEIN_STATE_ORDER),
    1.0 / len(PROTEIN_STATE_ORDER),
    dtype=float,
)
PROTEIN_GAP_CHARACTER = "-"
PROTEIN_MISSING_CHARACTER = "?"
_PROTEIN_ALLOWED_POLICY_VALUES = frozenset({"treat-as-missing", "reject"})


def normalize_unambiguous_protein_records(
    records: list[AlignmentRecord],
    *,
    model_name: str,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> list[AlignmentRecord]:
    """Uppercase aligned protein records and apply explicit gap and missing policy."""
    validated_gap_policy = validate_protein_gap_policy(
        gap_policy,
        model_name=model_name,
    )
    validated_missing_policy = validate_protein_missing_policy(
        missing_policy,
        model_name=model_name,
    )
    normalized_records: list[AlignmentRecord] = []
    for record in records:
        normalized_sequence = record.sequence.upper()
        invalid_states: set[str] = set()
        for state in normalized_sequence:
            if state in PROTEIN_STATE_INDEX:
                continue
            if state == PROTEIN_GAP_CHARACTER:
                if validated_gap_policy == "reject":
                    raise InvalidAlignmentError(
                        f"{model_name} likelihood gap policy rejects '-' in record '{record.identifier}'"
                    )
                continue
            if state == PROTEIN_MISSING_CHARACTER:
                if validated_missing_policy == "reject":
                    raise InvalidAlignmentError(
                        f"{model_name} likelihood missing-state policy rejects '?' in record '{record.identifier}'"
                    )
                continue
            invalid_states.add(state)
        if invalid_states:
            joined_states = ", ".join(sorted(invalid_states))
            raise InvalidAlignmentError(
                f"{model_name} likelihood currently requires unambiguous amino-acid states plus explicit gap/missing policy; "
                f"record '{record.identifier}' contains {joined_states}"
            )
        normalized_records.append(
            AlignmentRecord(
                identifier=record.identifier,
                sequence=normalized_sequence,
            )
        )
    return normalized_records


def validate_protein_gap_policy(
    gap_policy: str,
    *,
    model_name: str,
) -> str:
    if gap_policy not in _PROTEIN_ALLOWED_POLICY_VALUES:
        raise ValueError(
            f"{model_name} protein gap policy must be one of {sorted(_PROTEIN_ALLOWED_POLICY_VALUES)}"
        )
    return gap_policy


def validate_protein_missing_policy(
    missing_policy: str,
    *,
    model_name: str,
) -> str:
    if missing_policy not in _PROTEIN_ALLOWED_POLICY_VALUES:
        raise ValueError(
            f"{model_name} protein missing-state policy must be one of {sorted(_PROTEIN_ALLOWED_POLICY_VALUES)}"
        )
    return missing_policy


def validate_protein_root_prior(
    root_prior: numpy.ndarray | list[float] | tuple[float, ...],
    *,
    model_name: str,
) -> numpy.ndarray:
    vector = numpy.asarray(root_prior, dtype=float)
    if vector.shape != (len(PROTEIN_STATE_ORDER),):
        raise InvalidAlignmentError(
            f"{model_name} likelihood requires exactly twenty root-prior entries in protein state order"
        )
    if not numpy.all(numpy.isfinite(vector)):
        raise InvalidAlignmentError(
            f"{model_name} likelihood root prior must contain only finite values"
        )
    if numpy.any(vector < 0.0):
        raise InvalidAlignmentError(
            f"{model_name} likelihood root prior must be nonnegative"
        )
    total = float(vector.sum())
    if total <= 0.0:
        raise InvalidAlignmentError(
            f"{model_name} likelihood root prior must sum to a positive value"
        )
    return vector / total


def validate_empirical_protein_rate_matrix(
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    *,
    model_name: str,
    row_sum_tolerance: float = 1e-9,
) -> numpy.ndarray:
    candidate = numpy.asarray(rate_matrix, dtype=float)
    state_count = len(PROTEIN_STATE_ORDER)
    if candidate.shape != (state_count, state_count):
        raise InvalidAlignmentError(
            f"{model_name} likelihood requires one 20x20 empirical protein rate matrix"
        )
    if not numpy.all(numpy.isfinite(candidate)):
        raise InvalidAlignmentError(
            f"{model_name} likelihood empirical protein rate matrix must contain only finite values"
        )
    off_diagonal = candidate.copy()
    numpy.fill_diagonal(off_diagonal, 0.0)
    if numpy.any(off_diagonal < 0.0):
        raise InvalidAlignmentError(
            f"{model_name} likelihood empirical protein rate matrix must not contain negative off-diagonal rates"
        )
    diagonal = numpy.diag(candidate)
    if numpy.any(diagonal > row_sum_tolerance):
        raise InvalidAlignmentError(
            f"{model_name} likelihood empirical protein rate matrix diagonal must be zero or negative"
        )
    row_sums = candidate.sum(axis=1)
    if not numpy.allclose(
        row_sums,
        numpy.zeros(state_count, dtype=float),
        rtol=0.0,
        atol=row_sum_tolerance,
    ):
        raise InvalidAlignmentError(
            f"{model_name} likelihood empirical protein rate matrix rows must sum to zero"
        )
    if not numpy.any(off_diagonal > 0.0):
        raise InvalidAlignmentError(
            f"{model_name} likelihood empirical protein rate matrix must contain at least one positive off-diagonal rate"
        )
    return candidate.copy()


def protein_leaf_likelihood_vector(
    states_by_taxon: dict[str, str],
    *,
    model_name: str,
    node_name: str | None,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> numpy.ndarray:
    """Return a one-hot or all-states-allowed vector for one protein tip state."""
    if node_name is None:
        raise AlignmentTaxonMismatchError(
            f"{model_name} likelihood requires named tree tips for alignment lookup"
        )
    state = states_by_taxon[node_name]
    if state in PROTEIN_STATE_INDEX:
        vector = numpy.zeros(len(PROTEIN_STATE_ORDER), dtype=float)
        vector[PROTEIN_STATE_INDEX[state]] = 1.0
        return vector
    if state == PROTEIN_GAP_CHARACTER:
        if validate_protein_gap_policy(gap_policy, model_name=model_name) == "reject":
            raise InvalidAlignmentError(
                f"{model_name} likelihood gap policy rejects '-' at taxon '{node_name}'"
            )
        return numpy.ones(len(PROTEIN_STATE_ORDER), dtype=float)
    if state == PROTEIN_MISSING_CHARACTER:
        if (
            validate_protein_missing_policy(
                missing_policy,
                model_name=model_name,
            )
            == "reject"
        ):
            raise InvalidAlignmentError(
                f"{model_name} likelihood missing-state policy rejects '?' at taxon '{node_name}'"
            )
        return numpy.ones(len(PROTEIN_STATE_ORDER), dtype=float)
    raise InvalidAlignmentError(
        f"{model_name} likelihood encountered unsupported amino-acid state '{state}' at taxon '{node_name}'"
    )


def evaluate_fixed_topology_protein_likelihood_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    model_name: str,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...],
    transition_matrix_for_child: Callable[[object], numpy.ndarray],
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> float:
    """Evaluate one fixed-topology protein CTMC likelihood on compressed patterns."""
    validate_explicit_branch_lengths(tree, model_name=model_name)
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name=model_name,
    )
    validated_root_prior = validate_protein_root_prior(
        root_prior,
        model_name=model_name,
    )

    return sum_compressed_site_pattern_log_likelihoods(
        compressed_patterns,
        site_log_likelihood=lambda states: (
            evaluate_fixed_topology_protein_site_log_likelihood(
                tree,
                states,
                taxon_order=compressed_patterns.taxon_order,
                model_name=model_name,
                root_prior=validated_root_prior,
                transition_matrix_for_child=transition_matrix_for_child,
                gap_policy=gap_policy,
                missing_policy=missing_policy,
            )
        ),
    )


def evaluate_fixed_topology_protein_site_log_likelihood(
    tree: PhyloTree,
    states: tuple[str, ...],
    *,
    taxon_order: list[str],
    model_name: str,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...],
    transition_matrix_for_child: Callable[[object], numpy.ndarray],
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> float:
    """Evaluate one protein site log likelihood on one fixed topology."""
    states_by_taxon = dict(zip(taxon_order, states, strict=True))
    validated_root_prior = validate_protein_root_prior(
        root_prior,
        model_name=model_name,
    )
    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=len(PROTEIN_STATE_ORDER),
        leaf_likelihood=lambda node: protein_leaf_likelihood_vector(
            states_by_taxon,
            model_name=model_name,
            node_name=node.name,
            gap_policy=gap_policy,
            missing_policy=missing_policy,
        ),
        transition_matrix_for_child=transition_matrix_for_child,
    )
    return log_likelihood_from_root_prior(
        tree,
        pruning_pass,
        root_prior=validated_root_prior,
    )
