from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    evaluate_fixed_topology_dna_site_log_likelihood,
    normalize_unambiguous_dna_records,
)
from bijux_phylogenetics.phylo.likelihood.empirical import (
    _evaluate_empirical_protein_tree_likelihood_from_patterns,
    _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_patterns,
    _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_patterns,
    _evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_patterns,
    validate_empirical_branch_optimization_model,
)
from bijux_phylogenetics.phylo.likelihood.gamma import (
    build_discrete_gamma_rate_categories,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    BranchLikelihoodDiagnosticRow,
    FixedTreeBranchLikelihoodDiagnosticsReport,
    LocalClockLikelihoodReport,
    StrictClockLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    resolve_selected_nucleotide_likelihood_specification,
    validate_selected_nucleotide_likelihood_model,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.poisson import (
    _evaluate_protein_poisson_tree_likelihood_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.protein import (
    UNIFORM_PROTEIN_ROOT_PRIOR,
    normalize_unambiguous_protein_records,
    validate_empirical_protein_rate_matrix,
    validate_protein_root_prior,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_EMPIRICAL_BRANCH_DIAGNOSTIC_MODELS = frozenset(
    {"fixed-rate", "discrete-gamma", "invariant", "discrete-gamma-invariant"}
)
_NEAR_ZERO_BRANCH_LENGTH = 1e-8
_NEGLIGIBLE_DELTA_TOLERANCE = 1e-9

BranchLikelihoodDiagnosticsWritableReport = (
    FixedTreeBranchLikelihoodDiagnosticsReport
    | LocalClockLikelihoodReport
    | StrictClockLikelihoodReport
)


def evaluate_nucleotide_branch_likelihood_diagnostics(
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
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> FixedTreeBranchLikelihoodDiagnosticsReport:
    """Diagnose one fixed-topology nucleotide likelihood by branch-collapse replay."""
    normalized_model_name = validate_selected_nucleotide_likelihood_model(model_name)
    normalized_records = normalize_unambiguous_dna_records(
        records,
        model_name=normalized_model_name.upper(),
    )
    specification = resolve_selected_nucleotide_likelihood_specification(
        normalized_records,
        model_name=normalized_model_name,
        owner_name=f"{normalized_model_name.upper()} branch likelihood diagnostics",
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )

    def evaluate_tree_log_likelihood(candidate_tree: PhyloTree) -> float:
        validate_explicit_branch_lengths(
            candidate_tree,
            model_name=specification.model_name,
        )
        validate_tree_taxa_against_patterns(
            candidate_tree,
            compressed_patterns,
            model_name=specification.model_name,
        )
        total_log_likelihood = 0.0
        for pattern in compressed_patterns.patterns:
            total_log_likelihood += pattern.weight * (
                evaluate_fixed_topology_dna_site_log_likelihood(
                    candidate_tree,
                    pattern.states,
                    taxon_order=compressed_patterns.taxon_order,
                    model_name=specification.model_name,
                    observation_policy=specification.observation_policy,
                    root_prior=specification.root_prior,
                    transition_matrix_for_child=lambda child: (
                        specification.transition_matrix_for_branch_length(
                            max(float(child.branch_length or 0.0), 0.0)
                        )
                    ),
                )
            )
        return total_log_likelihood

    baseline_log_likelihood = evaluate_tree_log_likelihood(tree)
    return summarize_fixed_tree_branch_likelihood_diagnostics(
        tree,
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        model_name=specification.model_name,
        baseline_log_likelihood=baseline_log_likelihood,
        evaluate_tree_log_likelihood=evaluate_tree_log_likelihood,
    )


def evaluate_nucleotide_branch_likelihood_diagnostics_from_alignment(
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
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> FixedTreeBranchLikelihoodDiagnosticsReport:
    normalized_records = load_fasta_alignment(alignment_path)
    return evaluate_nucleotide_branch_likelihood_diagnostics(
        load_tree(tree_path),
        normalized_records,
        model_name=model_name,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )


def evaluate_protein_poisson_branch_likelihood_diagnostics(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> FixedTreeBranchLikelihoodDiagnosticsReport:
    """Diagnose one fixed-topology protein Poisson likelihood by branch collapse."""
    normalized_records = normalize_unambiguous_protein_records(
        records,
        model_name="protein Poisson",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )

    def evaluate_tree_log_likelihood(candidate_tree: PhyloTree) -> float:
        return _evaluate_protein_poisson_tree_likelihood_from_patterns(
            candidate_tree,
            compressed_patterns,
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            ambiguity_policy="reject",
        ).log_likelihood

    baseline_log_likelihood = evaluate_tree_log_likelihood(tree)
    return summarize_fixed_tree_branch_likelihood_diagnostics(
        tree,
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        model_name="protein Poisson",
        baseline_log_likelihood=baseline_log_likelihood,
        evaluate_tree_log_likelihood=evaluate_tree_log_likelihood,
    )


def evaluate_protein_poisson_branch_likelihood_diagnostics_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> FixedTreeBranchLikelihoodDiagnosticsReport:
    return evaluate_protein_poisson_branch_likelihood_diagnostics(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )


def evaluate_empirical_protein_branch_likelihood_diagnostics(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    likelihood_model: str,
    alpha: float | None = None,
    invariant_proportion: float | None = None,
    category_count: int = 4,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> FixedTreeBranchLikelihoodDiagnosticsReport:
    """Diagnose one fixed-topology empirical protein likelihood by branch collapse."""
    validated_likelihood_model = validate_empirical_branch_optimization_model(
        likelihood_model
    )
    if validated_likelihood_model not in _EMPIRICAL_BRANCH_DIAGNOSTIC_MODELS:
        raise ValueError(
            "empirical protein branch diagnostics likelihood_model must be one of "
            f"{sorted(_EMPIRICAL_BRANCH_DIAGNOSTIC_MODELS)}"
        )
    normalized_records = normalize_unambiguous_protein_records(
        records,
        model_name=f"empirical protein matrix {validated_likelihood_model}",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if root_prior is None:
        validated_root_prior = UNIFORM_PROTEIN_ROOT_PRIOR
        root_prior_source = "uniform"
    else:
        validated_root_prior = validate_protein_root_prior(
            root_prior,
            model_name=f"empirical protein matrix {validated_likelihood_model}",
        )
        root_prior_source = "provided"
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name=f"empirical protein matrix {validated_likelihood_model}",
    )
    model_label = _empirical_branch_diagnostic_model_label(validated_likelihood_model)

    def evaluate_tree_log_likelihood(candidate_tree: PhyloTree) -> float:
        if validated_likelihood_model == "fixed-rate":
            return _evaluate_empirical_protein_tree_likelihood_from_patterns(
                candidate_tree,
                compressed_patterns,
                rate_matrix=validated_rate_matrix,
                root_prior=validated_root_prior,
                root_prior_source=root_prior_source,
                matrix_label=matrix_label,
                gap_policy=gap_policy,
                missing_policy=missing_policy,
                ambiguity_policy="reject",
            ).log_likelihood
        if validated_likelihood_model == "discrete-gamma":
            if alpha is None:
                raise ValueError(
                    "empirical protein branch diagnostics with discrete-gamma require alpha"
                )
            return _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_patterns(
                candidate_tree,
                compressed_patterns,
                rate_matrix=validated_rate_matrix,
                alpha=alpha,
                category_count=category_count,
                root_prior=validated_root_prior,
                root_prior_source=root_prior_source,
                matrix_label=matrix_label,
                gap_policy=gap_policy,
                missing_policy=missing_policy,
            ).log_likelihood
        if validated_likelihood_model == "invariant":
            if invariant_proportion is None:
                raise ValueError(
                    "empirical protein branch diagnostics with invariant mixture require invariant_proportion"
                )
            transition_by_node_id = _empirical_transition_by_node_id(
                candidate_tree,
                validated_rate_matrix,
                rate_scale=1.0,
            )
            return _evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_patterns(
                candidate_tree,
                compressed_patterns,
                root_prior=validated_root_prior,
                root_prior_source=root_prior_source,
                matrix_label=matrix_label,
                gap_policy=gap_policy,
                missing_policy=missing_policy,
                invariant_proportion=invariant_proportion,
                initial_invariant_proportion=invariant_proportion,
                lower_invariant_proportion_bound=invariant_proportion,
                upper_invariant_proportion_bound=invariant_proportion,
                transition_by_node_id=transition_by_node_id,
                function_evaluation_count=1,
                converged=True,
            ).log_likelihood
        if invariant_proportion is None:
            raise ValueError(
                "empirical protein branch diagnostics with discrete-gamma-invariant require invariant_proportion"
            )
        if alpha is None:
            raise ValueError(
                "empirical protein branch diagnostics with discrete-gamma-invariant require alpha"
            )
        categories = build_discrete_gamma_rate_categories(
            alpha=alpha,
            category_count=category_count,
        )
        return _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_patterns(
            candidate_tree,
            compressed_patterns,
            validated_rate_matrix=validated_rate_matrix,
            alpha=alpha,
            categories=categories,
            root_prior=validated_root_prior,
            root_prior_source=root_prior_source,
            matrix_label=matrix_label,
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            invariant_proportion=invariant_proportion,
            initial_invariant_proportion=invariant_proportion,
            lower_invariant_proportion_bound=invariant_proportion,
            upper_invariant_proportion_bound=invariant_proportion,
            function_evaluation_count=1,
            converged=True,
            emit_boundary_warnings=False,
        ).log_likelihood

    baseline_log_likelihood = evaluate_tree_log_likelihood(tree)
    return summarize_fixed_tree_branch_likelihood_diagnostics(
        tree,
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        model_name=model_label,
        baseline_log_likelihood=baseline_log_likelihood,
        evaluate_tree_log_likelihood=evaluate_tree_log_likelihood,
    )


def evaluate_empirical_protein_branch_likelihood_diagnostics_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    likelihood_model: str,
    alpha: float | None = None,
    invariant_proportion: float | None = None,
    category_count: int = 4,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> FixedTreeBranchLikelihoodDiagnosticsReport:
    return evaluate_empirical_protein_branch_likelihood_diagnostics(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        rate_matrix=rate_matrix,
        likelihood_model=likelihood_model,
        alpha=alpha,
        invariant_proportion=invariant_proportion,
        category_count=category_count,
        root_prior=root_prior,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )


def summarize_fixed_tree_branch_likelihood_diagnostics(
    tree: PhyloTree,
    *,
    taxa: list[str],
    site_count: int,
    pattern_count: int,
    model_name: str,
    baseline_log_likelihood: float,
    evaluate_tree_log_likelihood,
) -> FixedTreeBranchLikelihoodDiagnosticsReport:
    """Summarize one fixed-tree likelihood surface by replaying branch collapse."""
    branch_diagnostics: list[BranchLikelihoodDiagnosticRow] = []
    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise ValueError("branch likelihood diagnostics require stable node ids")
        collapsed_tree = tree.copy()
        collapsed_branch = collapsed_tree.node_by_id(child.node_id)
        collapsed_branch.branch_length = 0.0
        collapsed_branch_log_likelihood = float(
            evaluate_tree_log_likelihood(collapsed_tree)
        )
        contribution_proxy = _contribution_proxy(
            baseline_log_likelihood=baseline_log_likelihood,
            collapsed_branch_log_likelihood=collapsed_branch_log_likelihood,
        )
        branch_diagnostics.append(
            BranchLikelihoodDiagnosticRow(
                branch_id=child.node_id,
                child_name=child.name,
                descendant_taxa=child.descendant_taxa,
                branch_length=float(child.branch_length or 0.0),
                collapsed_branch_log_likelihood=collapsed_branch_log_likelihood,
                contribution_proxy=contribution_proxy,
                warning_flags=_branch_likelihood_warning_flags(
                    branch_length=float(child.branch_length or 0.0),
                    baseline_log_likelihood=baseline_log_likelihood,
                    collapsed_branch_log_likelihood=collapsed_branch_log_likelihood,
                    contribution_proxy=contribution_proxy,
                ),
            )
        )
    return FixedTreeBranchLikelihoodDiagnosticsReport(
        model_name=model_name,
        taxa=list(taxa),
        site_count=site_count,
        pattern_count=pattern_count,
        branch_count=len(branch_diagnostics),
        tree_newick=dumps_newick(tree),
        baseline_log_likelihood=baseline_log_likelihood,
        branch_diagnostics=branch_diagnostics,
    )


def write_branch_likelihood_diagnostic_table(
    path: Path,
    report: BranchLikelihoodDiagnosticsWritableReport,
) -> Path:
    """Write one branchwise likelihood-diagnostic TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    model_name = _branch_diagnostic_model_name(report)
    baseline_log_likelihood = _branch_diagnostic_baseline_log_likelihood(report)
    branch_diagnostics = _branch_diagnostic_rows(report)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "model_name",
                "branch_id",
                "child_name",
                "descendant_taxa",
                "branch_length",
                "baseline_log_likelihood",
                "collapsed_branch_log_likelihood",
                "contribution_proxy",
                "warning_flags",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in branch_diagnostics:
            writer.writerow(
                {
                    "model_name": model_name,
                    "branch_id": row.branch_id,
                    "child_name": row.child_name or "",
                    "descendant_taxa": "|".join(row.descendant_taxa),
                    "branch_length": format(row.branch_length, ".15g"),
                    "baseline_log_likelihood": repr(baseline_log_likelihood),
                    "collapsed_branch_log_likelihood": repr(
                        row.collapsed_branch_log_likelihood
                    ),
                    "contribution_proxy": repr(row.contribution_proxy),
                    "warning_flags": "|".join(row.warning_flags),
                }
            )
    return path


def _branch_diagnostic_model_name(
    report: BranchLikelihoodDiagnosticsWritableReport,
) -> str:
    return report.model_name


def _branch_diagnostic_baseline_log_likelihood(
    report: BranchLikelihoodDiagnosticsWritableReport,
) -> float:
    if isinstance(report, FixedTreeBranchLikelihoodDiagnosticsReport):
        return report.baseline_log_likelihood
    if isinstance(report, StrictClockLikelihoodReport):
        return report.optimized_log_likelihood
    return report.optimized_log_likelihood


def _branch_diagnostic_rows(
    report: BranchLikelihoodDiagnosticsWritableReport,
) -> list[BranchLikelihoodDiagnosticRow]:
    if isinstance(report, FixedTreeBranchLikelihoodDiagnosticsReport):
        return list(report.branch_diagnostics)
    return list(report.branch_likelihood_diagnostics)


def _contribution_proxy(
    *,
    baseline_log_likelihood: float,
    collapsed_branch_log_likelihood: float,
) -> float:
    if not math.isfinite(baseline_log_likelihood):
        return float("nan")
    if not math.isfinite(collapsed_branch_log_likelihood):
        return float("inf")
    return baseline_log_likelihood - collapsed_branch_log_likelihood


def _branch_likelihood_warning_flags(
    *,
    branch_length: float,
    baseline_log_likelihood: float,
    collapsed_branch_log_likelihood: float,
    contribution_proxy: float,
) -> list[str]:
    warning_flags: list[str] = []
    if branch_length == 0.0:
        warning_flags.append("zero-branch-length")
    elif branch_length < _NEAR_ZERO_BRANCH_LENGTH:
        warning_flags.append("near-zero-branch-length")
    if not math.isfinite(collapsed_branch_log_likelihood):
        warning_flags.append("collapsed-log-likelihood-nonfinite")
        return warning_flags
    if not math.isfinite(contribution_proxy):
        warning_flags.append("contribution-proxy-nonfinite")
        return warning_flags
    if collapsed_branch_log_likelihood > (
        baseline_log_likelihood + _NEGLIGIBLE_DELTA_TOLERANCE
    ):
        warning_flags.append("collapse-improves-likelihood")
    elif abs(contribution_proxy) <= _NEGLIGIBLE_DELTA_TOLERANCE:
        warning_flags.append("collapse-has-negligible-effect")
    return warning_flags


def _empirical_branch_diagnostic_model_label(likelihood_model: str) -> str:
    return {
        "fixed-rate": "empirical protein matrix",
        "discrete-gamma": "empirical protein matrix +G",
        "invariant": "empirical protein matrix +I",
        "discrete-gamma-invariant": "empirical protein matrix +G+I",
    }[likelihood_model]


def _empirical_transition_by_node_id(
    tree: PhyloTree,
    validated_rate_matrix: numpy.ndarray,
    *,
    rate_scale: float,
) -> dict[str | None, numpy.ndarray]:
    from bijux_phylogenetics.phylo.likelihood.pruning import (
        transition_probability_matrix,
    )

    return {
        child.node_id: transition_probability_matrix(
            validated_rate_matrix,
            max(float(child.branch_length or 0.0), 0.0) * rate_scale,
        )
        for _parent, child in tree.iter_edges()
    }
