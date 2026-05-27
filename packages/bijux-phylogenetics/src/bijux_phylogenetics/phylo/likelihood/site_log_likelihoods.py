from __future__ import annotations

from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    UNIFORM_DNA_ROOT_PRIOR,
    estimate_empirical_dna_base_frequencies,
    evaluate_fixed_topology_dna_site_log_likelihood,
    normalize_dna_exchangeabilities_by_anchor,
    normalize_unambiguous_dna_records,
    validate_dna_base_frequencies,
    validate_positive_kappa,
)
from bijux_phylogenetics.phylo.likelihood.f81 import (
    f81_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.gtr import (
    gtr_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.hky85 import (
    hky85_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.jc69 import (
    jc69_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.k80 import (
    k80_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    FixedTopologySiteLogLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.sites import (
    expanded_site_log_likelihood_rows_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_NUCLEOTIDE_SITE_LOG_LIKELIHOOD_MODELS = frozenset(
    {"jc69", "k80", "f81", "hky85", "gtr"}
)


def validate_nucleotide_site_log_likelihood_model(model_name: str) -> str:
    normalized_model_name = model_name.strip().lower()
    if normalized_model_name not in _NUCLEOTIDE_SITE_LOG_LIKELIHOOD_MODELS:
        raise ValueError(
            "nucleotide site-log-likelihood model must be one of "
            + ", ".join(sorted(_NUCLEOTIDE_SITE_LOG_LIKELIHOOD_MODELS))
        )
    return normalized_model_name


def evaluate_nucleotide_site_log_likelihoods(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    kappa: float | None = None,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
) -> FixedTopologySiteLogLikelihoodReport:
    """Evaluate one fixed-topology nucleotide likelihood with expanded site rows."""
    normalized_model_name = validate_nucleotide_site_log_likelihood_model(model_name)
    normalized_records = normalize_unambiguous_dna_records(
        records,
        model_name=normalized_model_name.upper(),
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(normalized_records)
    if normalized_model_name == "jc69":
        _reject_irrelevant_parameter("JC69 site log likelihood export", "kappa", kappa)
        _reject_irrelevant_parameter(
            "JC69 site log likelihood export",
            "base_frequencies",
            base_frequencies,
        )
        _reject_irrelevant_parameter(
            "JC69 site log likelihood export",
            "exchangeabilities",
            exchangeabilities,
        )
        return _evaluate_selected_dna_site_log_likelihoods_from_patterns(
            tree,
            compressed_patterns,
            model_name="JC69",
            root_prior=UNIFORM_DNA_ROOT_PRIOR,
            parameter_values={},
            transition_matrix_for_child=lambda child: jc69_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0)
            ),
        )
    if normalized_model_name == "k80":
        _reject_irrelevant_parameter(
            "K80 site log likelihood export",
            "base_frequencies",
            base_frequencies,
        )
        _reject_irrelevant_parameter(
            "K80 site log likelihood export",
            "exchangeabilities",
            exchangeabilities,
        )
        if kappa is None:
            raise ValueError("K80 site log likelihood export requires 'kappa'")
        validated_kappa = validate_positive_kappa(kappa, model_name="K80")
        return _evaluate_selected_dna_site_log_likelihoods_from_patterns(
            tree,
            compressed_patterns,
            model_name="K80",
            root_prior=UNIFORM_DNA_ROOT_PRIOR,
            parameter_values={"kappa": validated_kappa},
            transition_matrix_for_child=lambda child: k80_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                kappa=validated_kappa,
            ),
        )
    if normalized_model_name == "f81":
        _reject_irrelevant_parameter("F81 site log likelihood export", "kappa", kappa)
        _reject_irrelevant_parameter(
            "F81 site log likelihood export",
            "exchangeabilities",
            exchangeabilities,
        )
        stationary = _resolve_dna_base_frequencies(
            normalized_records,
            base_frequencies=base_frequencies,
            model_name="F81",
        )
        return _evaluate_selected_dna_site_log_likelihoods_from_patterns(
            tree,
            compressed_patterns,
            model_name="F81",
            root_prior=stationary,
            parameter_values=_base_frequency_parameter_values(stationary),
            transition_matrix_for_child=lambda child: f81_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                base_frequencies=stationary,
            ),
        )
    if normalized_model_name == "hky85":
        _reject_irrelevant_parameter(
            "HKY85 site log likelihood export",
            "exchangeabilities",
            exchangeabilities,
        )
        if kappa is None:
            raise ValueError("HKY85 site log likelihood export requires 'kappa'")
        stationary = _resolve_dna_base_frequencies(
            normalized_records,
            base_frequencies=base_frequencies,
            model_name="HKY85",
        )
        validated_kappa = validate_positive_kappa(kappa, model_name="HKY85")
        return _evaluate_selected_dna_site_log_likelihoods_from_patterns(
            tree,
            compressed_patterns,
            model_name="HKY85",
            root_prior=stationary,
            parameter_values={
                **_base_frequency_parameter_values(stationary),
                "kappa": validated_kappa,
            },
            transition_matrix_for_child=lambda child: hky85_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                base_frequencies=stationary,
                kappa=validated_kappa,
            ),
        )
    _reject_irrelevant_parameter("GTR site log likelihood export", "kappa", kappa)
    if exchangeabilities is None:
        raise ValueError("GTR site log likelihood export requires 'exchangeabilities'")
    stationary = _resolve_dna_base_frequencies(
        normalized_records,
        base_frequencies=base_frequencies,
        model_name="GTR",
    )
    normalized_exchangeabilities = normalize_dna_exchangeabilities_by_anchor(
        exchangeabilities,
        model_name="GTR",
    )
    return _evaluate_selected_dna_site_log_likelihoods_from_patterns(
        tree,
        compressed_patterns,
        model_name="GTR",
        root_prior=stationary,
        parameter_values={
            **_base_frequency_parameter_values(stationary),
            "exchangeability_ac": float(normalized_exchangeabilities[0]),
            "exchangeability_ag": float(normalized_exchangeabilities[1]),
            "exchangeability_at": float(normalized_exchangeabilities[2]),
            "exchangeability_cg": float(normalized_exchangeabilities[3]),
            "exchangeability_ct": float(normalized_exchangeabilities[4]),
            "exchangeability_gt": float(normalized_exchangeabilities[5]),
        },
        transition_matrix_for_child=lambda child: gtr_transition_probability_matrix(
            max(float(child.branch_length or 0.0), 0.0),
            exchangeabilities=normalized_exchangeabilities,
            base_frequencies=stationary,
        ),
    )


def evaluate_nucleotide_site_log_likelihoods_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    kappa: float | None = None,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
) -> FixedTopologySiteLogLikelihoodReport:
    """Evaluate one fixed-topology nucleotide likelihood from paths with site rows."""
    return evaluate_nucleotide_site_log_likelihoods(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        model_name=model_name,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
    )


def _evaluate_selected_dna_site_log_likelihoods_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    model_name: str,
    root_prior: numpy.ndarray,
    parameter_values: dict[str, float],
    transition_matrix_for_child,
) -> FixedTopologySiteLogLikelihoodReport:
    validate_explicit_branch_lengths(tree, model_name=model_name)
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name=model_name,
    )
    site_log_likelihoods, total_log_likelihood = (
        expanded_site_log_likelihood_rows_from_patterns(
            compressed_patterns,
            site_log_likelihood=lambda states: evaluate_fixed_topology_dna_site_log_likelihood(
                tree,
                states,
                taxon_order=compressed_patterns.taxon_order,
                model_name=model_name,
                root_prior=root_prior,
                transition_matrix_for_child=transition_matrix_for_child,
            ),
        )
    )
    return FixedTopologySiteLogLikelihoodReport(
        model_name=model_name,
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        expansion_policy="expanded-site-rows",
        tree_newick=dumps_newick(tree),
        parameter_values=parameter_values,
        log_likelihood=total_log_likelihood,
        site_log_likelihoods=site_log_likelihoods,
    )


def _resolve_dna_base_frequencies(
    records: list[AlignmentRecord],
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None,
    model_name: str,
) -> numpy.ndarray:
    if base_frequencies is None:
        return estimate_empirical_dna_base_frequencies(records)
    return validate_dna_base_frequencies(base_frequencies, model_name=model_name)


def _base_frequency_parameter_values(
    stationary: numpy.ndarray,
) -> dict[str, float]:
    return {
        "base_frequency_a": float(stationary[0]),
        "base_frequency_c": float(stationary[1]),
        "base_frequency_g": float(stationary[2]),
        "base_frequency_t": float(stationary[3]),
    }


def _reject_irrelevant_parameter(
    owner_name: str,
    parameter_name: str,
    value: object,
) -> None:
    if value is not None:
        raise ValueError(
            f"{owner_name} does not accept '{parameter_name}' because that model does not use it"
        )
