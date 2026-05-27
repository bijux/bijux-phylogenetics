from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import load_newick_tree_set
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import (
    CandidateTreeSiteLikelihoodMatrixReport,
    CandidateTreeSiteLikelihoodRow,
    CandidateTreeSiteLikelihoodSummary,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    resolve_selected_nucleotide_likelihood_specification,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    validate_selected_nucleotide_likelihood_model,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.site_log_likelihoods import (
    evaluate_selected_dna_site_log_likelihoods_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.dna import normalize_unambiguous_dna_records
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import TreeParseError


@dataclass(frozen=True, slots=True)
class _CandidateTreeRecord:
    """One resolved candidate tree with a stable identifier inside matrix evaluation."""

    candidate_tree_id: str
    candidate_tree_label: str
    tree: PhyloTree


def evaluate_nucleotide_candidate_tree_site_likelihood_matrix(
    trees: list[PhyloTree] | Path,
    records: list[AlignmentRecord] | Path,
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
) -> CandidateTreeSiteLikelihoodMatrixReport:
    """Evaluate one shared nucleotide likelihood model across multiple candidate trees."""
    resolved_tree_records, resolved_tree_set_path = resolve_candidate_tree_records(trees)
    resolved_records, resolved_alignment_path = resolve_candidate_tree_alignment_records(
        records
    )
    normalized_model_name = validate_selected_nucleotide_likelihood_model(model_name)
    normalized_records = normalize_unambiguous_dna_records(
        resolved_records,
        model_name=normalized_model_name.upper(),
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(normalized_records)
    specification = resolve_selected_nucleotide_likelihood_specification(
        normalized_records,
        model_name=normalized_model_name,
        owner_name=f"{normalized_model_name.upper()} candidate tree site likelihood matrix",
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
    )

    candidate_trees: list[CandidateTreeSiteLikelihoodSummary] = []
    matrix_rows: list[CandidateTreeSiteLikelihoodRow] = []
    for candidate in resolved_tree_records:
        try:
            tree_report = evaluate_selected_dna_site_log_likelihoods_from_patterns(
                candidate.tree,
                compressed_patterns,
                model_name=specification.model_name,
                root_prior=specification.root_prior,
                parameter_values=specification.parameter_values,
                transition_matrix_for_child=lambda child: (
                    specification.transition_matrix_for_branch_length(
                        max(float(child.branch_length or 0.0), 0.0)
                    )
                ),
            )
        except ValueError as error:
            raise ValueError(
                f"{candidate.candidate_tree_label} is incompatible with the shared alignment/model surface: {error}"
            ) from error
        candidate_trees.append(
            CandidateTreeSiteLikelihoodSummary(
                candidate_tree_id=candidate.candidate_tree_id,
                candidate_tree_label=candidate.candidate_tree_label,
                tree_newick=tree_report.tree_newick,
                log_likelihood=tree_report.log_likelihood,
            )
        )
        matrix_rows.extend(
            CandidateTreeSiteLikelihoodRow(
                candidate_tree_id=candidate.candidate_tree_id,
                candidate_tree_label=candidate.candidate_tree_label,
                tree_newick=tree_report.tree_newick,
                pattern_id=row.pattern_id,
                pattern_weight=row.pattern_weight,
                site_position=row.site_position,
                site_states=row.site_states,
                log_likelihood=row.log_likelihood,
            )
            for row in tree_report.site_log_likelihoods
        )
    return CandidateTreeSiteLikelihoodMatrixReport(
        model_name=specification.model_name,
        tree_set_path=None if resolved_tree_set_path is None else str(resolved_tree_set_path),
        alignment_path=None
        if resolved_alignment_path is None
        else str(resolved_alignment_path),
        taxa=compressed_patterns.taxon_order,
        tree_count=len(candidate_trees),
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        expansion_policy="candidate-tree-expanded-site-rows",
        parameter_values=specification.parameter_values,
        candidate_trees=candidate_trees,
        matrix_rows=matrix_rows,
    )


def evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment(
    tree_set_path: Path,
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
) -> CandidateTreeSiteLikelihoodMatrixReport:
    """Evaluate a candidate-tree site-likelihood matrix from tree-set and alignment paths."""
    return evaluate_nucleotide_candidate_tree_site_likelihood_matrix(
        tree_set_path,
        alignment_path,
        model_name=model_name,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
    )


def resolve_candidate_tree_records(
    trees: list[PhyloTree] | Path,
) -> tuple[list[_CandidateTreeRecord], Path | None]:
    """Resolve one candidate-tree set and assign stable candidate identifiers."""
    resolved_tree_set_path = trees if isinstance(trees, Path) else None
    resolved_trees = load_newick_tree_set(trees) if isinstance(trees, Path) else list(trees)
    if len(resolved_trees) < 2:
        raise TreeParseError(
            "candidate tree site likelihood matrix requires at least two candidate trees"
        )
    return [
        _CandidateTreeRecord(
            candidate_tree_id=f"candidate-tree-{index}",
            candidate_tree_label=f"candidate-tree-{index}",
            tree=tree.copy().refresh(),
        )
        for index, tree in enumerate(resolved_trees, start=1)
    ], resolved_tree_set_path


def resolve_candidate_tree_alignment_records(
    records: list[AlignmentRecord] | Path,
) -> tuple[list[AlignmentRecord], Path | None]:
    """Resolve one candidate-tree alignment input and preserve its path when present."""
    resolved_alignment_path = records if isinstance(records, Path) else None
    resolved_records = (
        load_fasta_alignment(records) if isinstance(records, Path) else list(records)
    )
    return resolved_records, resolved_alignment_path
