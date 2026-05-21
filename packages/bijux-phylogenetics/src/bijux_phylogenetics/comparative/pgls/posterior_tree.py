from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import statistics
import tempfile

from Bio import Phylo

from bijux_phylogenetics.comparative.pgls.design import inspect_pgls_inputs
from bijux_phylogenetics.comparative.pgls.fitting import run_pgls
from bijux_phylogenetics.comparative.pgls.models import (
    ComparativeFormulaSpecification,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    ComparativeMethodError,
    InvalidAlignmentError,
)
from bijux_phylogenetics.trees import load_tree_set


@dataclass(slots=True)
class PosteriorTreePGLSTreeFitRow:
    """One per-tree PGLS fit within a posterior-tree comparative run."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    unrooted_topology_id: str
    lambda_value: float
    log_likelihood: float


@dataclass(slots=True)
class PosteriorTreePGLSCoefficientRow:
    """One coefficient observation from one posterior-tree PGLS fit."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    term: str
    estimate: float
    p_value: float
    significant: bool
    direction: str


@dataclass(slots=True)
class PosteriorTreePGLSCoefficientSummaryRow:
    """One reviewer-facing coefficient summary across the retained tree set."""

    term: str
    tree_fit_count: int
    positive_tree_count: int
    negative_tree_count: int
    zero_tree_count: int
    dominant_direction: str
    direction_consistency: float
    significant_tree_count: int
    significance_fraction: float
    conclusion_stability: str
    mean_estimate: float
    median_estimate: float
    standard_deviation: float
    minimum_estimate: float
    maximum_estimate: float
    lower_95_empirical_estimate: float
    upper_95_empirical_estimate: float
    mean_p_value: float
    median_p_value: float
    minimum_p_value: float
    maximum_p_value: float


@dataclass(slots=True)
class PosteriorTreePGLSReport:
    """PGLS fit stability across a posterior or bootstrap tree set."""

    tree_set_path: Path
    traits_path: Path
    response: str
    formula: ComparativeFormulaSpecification
    lambda_mode: str
    burnin_fraction: float
    significance_threshold: float
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    shared_tree_taxa: list[str]
    analysis_taxa: list[str]
    rooted_topology_count: int
    unrooted_topology_count: int
    tree_rows: list[PosteriorTreePGLSTreeFitRow]
    coefficient_rows: list[PosteriorTreePGLSCoefficientRow]
    coefficient_summaries: list[PosteriorTreePGLSCoefficientSummaryRow]
    warnings: list[str]


def run_posterior_tree_pgls(
    tree_set_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
    burnin_fraction: float = 0.0,
    significance_threshold: float = 0.05,
) -> PosteriorTreePGLSReport:
    """Fit the same PGLS model across every retained tree in a posterior or bootstrap set."""
    if not 0.0 <= burnin_fraction < 1.0:
        raise ComparativeMethodError(
            "burnin fraction must be between 0 inclusive and 1 exclusive"
        )
    if not 0.0 < significance_threshold < 1.0:
        raise ComparativeMethodError("significance threshold must be between 0 and 1")
    if lambda_value != "estimate":
        lambda_value = float(lambda_value)
        if not 0.0 <= lambda_value <= 1.0:
            raise ComparativeMethodError(
                "PGLS lambda value must be between 0 and 1 inclusive"
            )

    _source_format, trees = _load_tree_set_trees(tree_set_path)
    total_tree_count = len(trees)
    burnin_tree_count = math.floor(total_tree_count * burnin_fraction)
    kept_tree_entries = [
        (source_tree_index, tree)
        for source_tree_index, tree in enumerate(trees, start=1)
    ][burnin_tree_count:]
    if not kept_tree_entries:
        raise ComparativeMethodError(
            "posterior-tree PGLS retains no trees after burn-in removal"
        )

    kept_trees = [tree for _, tree in kept_tree_entries]
    shared_tree_taxa = sorted(_shared_taxa(kept_trees))
    if len(shared_tree_taxa) < 3:
        raise ComparativeMethodError(
            "posterior-tree PGLS requires at least three shared taxa across retained trees"
        )
    warnings: list[str] = []
    if any(set(tree.tip_names) != set(shared_tree_taxa) for tree in kept_trees):
        warnings.append(
            "retained trees do not share identical tip sets and were reduced to their shared taxa"
        )

    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-posterior-tree-pgls-"
    ) as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        reference_tree_path = tmp_dir_path / "posterior-tree-pgls-reference-tree.nwk"
        reference_tree = _prune_tree_to_taxa(
            kept_tree_entries[0][1],
            shared_tree_taxa,
            scratch_path=reference_tree_path,
        )
        reference_tree_path.write_text(
            dumps_newick(reference_tree) + "\n",
            encoding="utf-8",
        )
        input_report = inspect_pgls_inputs(
            reference_tree_path,
            traits_path,
            response=response,
            predictors=predictors,
            formula=formula,
            taxon_column=taxon_column,
        )
        if not input_report.ready:
            raise ComparativeMethodError("; ".join(input_report.blockers))

        analysis_taxa = list(input_report.analysis_taxa)
        analysis_tree_set_path = tmp_dir_path / "posterior-tree-pgls-analysis-set.nwk"
        analysis_trees = [
            _prune_tree_to_taxa(
                tree,
                analysis_taxa,
                scratch_path=tmp_dir_path / f"posterior-tree-pgls-pruned-{index}.nwk",
            )
            for index, (_source_tree_index, tree) in enumerate(
                kept_tree_entries, start=1
            )
        ]
        analysis_tree_set_path.write_text(
            "".join(dumps_newick(tree) + "\n" for tree in analysis_trees),
            encoding="utf-8",
        )
        topology_summary = load_tree_set(analysis_tree_set_path)

        tree_rows: list[PosteriorTreePGLSTreeFitRow] = []
        coefficient_rows: list[PosteriorTreePGLSCoefficientRow] = []
        single_tree_path = tmp_dir_path / "posterior-tree-pgls-current-tree.nwk"
        for topology_record, (source_tree_index, tree) in zip(
            topology_summary.records,
            kept_tree_entries,
            strict=True,
        ):
            pruned_tree = _prune_tree_to_taxa(
                tree,
                analysis_taxa,
                scratch_path=single_tree_path,
            )
            single_tree_path.write_text(
                dumps_newick(pruned_tree) + "\n",
                encoding="utf-8",
            )
            fit = run_pgls(
                single_tree_path,
                traits_path,
                response=response,
                predictors=predictors,
                formula=formula,
                taxon_column=taxon_column,
                lambda_value=lambda_value,
            )
            tree_rows.append(
                PosteriorTreePGLSTreeFitRow(
                    source_tree_index=source_tree_index,
                    post_burnin_index=topology_record.index,
                    rooted_topology_id=topology_record.rooted_topology_id,
                    unrooted_topology_id=topology_record.unrooted_topology_id,
                    lambda_value=fit.lambda_value,
                    log_likelihood=fit.log_likelihood,
                )
            )
            for coefficient in fit.coefficients:
                coefficient_rows.append(
                    PosteriorTreePGLSCoefficientRow(
                        source_tree_index=source_tree_index,
                        post_burnin_index=topology_record.index,
                        rooted_topology_id=topology_record.rooted_topology_id,
                        term=coefficient.name,
                        estimate=coefficient.estimate,
                        p_value=coefficient.p_value,
                        significant=coefficient.p_value <= significance_threshold,
                        direction=_estimate_direction(coefficient.estimate),
                    )
                )

    coefficient_summaries = _summarize_coefficients(
        coefficient_rows,
        significance_threshold=significance_threshold,
        term_order=list(input_report.encoded_columns),
    )
    if any(
        row.conclusion_stability == "direction_conflict"
        for row in coefficient_summaries
    ):
        warnings.append(
            "one or more coefficients change direction across retained trees"
        )
    if any(
        row.conclusion_stability == "mixed_support" for row in coefficient_summaries
    ):
        warnings.append(
            "one or more coefficients keep one direction but not one stable support decision across retained trees"
        )

    return PosteriorTreePGLSReport(
        tree_set_path=tree_set_path,
        traits_path=traits_path,
        response=input_report.response,
        formula=input_report.formula,
        lambda_mode="estimate" if lambda_value == "estimate" else "fixed",
        burnin_fraction=burnin_fraction,
        significance_threshold=significance_threshold,
        total_tree_count=total_tree_count,
        burnin_tree_count=burnin_tree_count,
        kept_tree_count=len(kept_tree_entries),
        shared_tree_taxa=shared_tree_taxa,
        analysis_taxa=analysis_taxa,
        rooted_topology_count=topology_summary.rooted_topology_count,
        unrooted_topology_count=topology_summary.unrooted_topology_count,
        tree_rows=tree_rows,
        coefficient_rows=coefficient_rows,
        coefficient_summaries=coefficient_summaries,
        warnings=warnings,
    )


def write_posterior_tree_pgls_tree_table(
    path: Path,
    report: PosteriorTreePGLSReport,
) -> Path:
    """Write one per-tree fit ledger for posterior-tree PGLS."""
    return write_taxon_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
            "lambda_value",
            "log_likelihood",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
                "lambda_value": format(row.lambda_value, ".15g"),
                "log_likelihood": format(row.log_likelihood, ".15g"),
            }
            for row in report.tree_rows
        ],
    )


def write_posterior_tree_pgls_coefficient_table(
    path: Path,
    report: PosteriorTreePGLSReport,
) -> Path:
    """Write one per-tree coefficient ledger for posterior-tree PGLS."""
    return write_taxon_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "term",
            "estimate",
            "p_value",
            "significant",
            "direction",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "term": row.term,
                "estimate": format(row.estimate, ".15g"),
                "p_value": format(row.p_value, ".15g"),
                "significant": str(row.significant).lower(),
                "direction": row.direction,
            }
            for row in report.coefficient_rows
        ],
    )


def write_posterior_tree_pgls_summary_table(
    path: Path,
    report: PosteriorTreePGLSReport,
) -> Path:
    """Write one coefficient-distribution summary ledger for posterior-tree PGLS."""
    return write_taxon_rows(
        path,
        columns=[
            "term",
            "tree_fit_count",
            "positive_tree_count",
            "negative_tree_count",
            "zero_tree_count",
            "dominant_direction",
            "direction_consistency",
            "significant_tree_count",
            "significance_fraction",
            "conclusion_stability",
            "mean_estimate",
            "median_estimate",
            "standard_deviation",
            "minimum_estimate",
            "maximum_estimate",
            "lower_95_empirical_estimate",
            "upper_95_empirical_estimate",
            "mean_p_value",
            "median_p_value",
            "minimum_p_value",
            "maximum_p_value",
        ],
        rows=[
            {
                "term": row.term,
                "tree_fit_count": str(row.tree_fit_count),
                "positive_tree_count": str(row.positive_tree_count),
                "negative_tree_count": str(row.negative_tree_count),
                "zero_tree_count": str(row.zero_tree_count),
                "dominant_direction": row.dominant_direction,
                "direction_consistency": format(row.direction_consistency, ".15g"),
                "significant_tree_count": str(row.significant_tree_count),
                "significance_fraction": format(row.significance_fraction, ".15g"),
                "conclusion_stability": row.conclusion_stability,
                "mean_estimate": format(row.mean_estimate, ".15g"),
                "median_estimate": format(row.median_estimate, ".15g"),
                "standard_deviation": format(row.standard_deviation, ".15g"),
                "minimum_estimate": format(row.minimum_estimate, ".15g"),
                "maximum_estimate": format(row.maximum_estimate, ".15g"),
                "lower_95_empirical_estimate": format(
                    row.lower_95_empirical_estimate, ".15g"
                ),
                "upper_95_empirical_estimate": format(
                    row.upper_95_empirical_estimate, ".15g"
                ),
                "mean_p_value": format(row.mean_p_value, ".15g"),
                "median_p_value": format(row.median_p_value, ".15g"),
                "minimum_p_value": format(row.minimum_p_value, ".15g"),
                "maximum_p_value": format(row.maximum_p_value, ".15g"),
            }
            for row in report.coefficient_summaries
        ],
    )


def _load_tree_set_trees(path: Path) -> tuple[str, list[PhyloTree]]:
    if not path.exists():
        raise FileNotFoundError(f"tree-set file not found: {path}")
    source_format = detect_tree_format(path)
    bio_trees = list(Phylo.parse(path, source_format))
    if not bio_trees:
        raise InvalidAlignmentError(f"tree set contains no trees: {path}")
    return source_format, [
        tree_from_biophylo(tree, source_format=source_format) for tree in bio_trees
    ]


def _shared_taxa(trees: list[PhyloTree]) -> set[str]:
    shared = set(trees[0].tip_names)
    for tree in trees[1:]:
        shared &= set(tree.tip_names)
    return shared


def _prune_tree_to_taxa(
    tree: PhyloTree,
    requested_taxa: list[str],
    *,
    scratch_path: Path,
) -> PhyloTree:
    scratch_path.write_text(dumps_newick(tree) + "\n", encoding="utf-8")
    pruned_tree, _report = prune_tree_to_requested_taxa(scratch_path, requested_taxa)
    return pruned_tree


def _summarize_coefficients(
    rows: list[PosteriorTreePGLSCoefficientRow],
    *,
    significance_threshold: float,
    term_order: list[str],
) -> list[PosteriorTreePGLSCoefficientSummaryRow]:
    grouped: dict[str, list[PosteriorTreePGLSCoefficientRow]] = {}
    for row in rows:
        grouped.setdefault(row.term, []).append(row)
    ordered_terms = term_order + sorted(
        term for term in grouped if term not in term_order
    )
    summaries: list[PosteriorTreePGLSCoefficientSummaryRow] = []
    for term in ordered_terms:
        term_rows = grouped.get(term)
        if not term_rows:
            continue
        estimates = [row.estimate for row in term_rows]
        p_values = [row.p_value for row in term_rows]
        positive_tree_count = sum(1 for row in term_rows if row.direction == "positive")
        negative_tree_count = sum(1 for row in term_rows if row.direction == "negative")
        zero_tree_count = sum(1 for row in term_rows if row.direction == "zero")
        direction_counts = {
            "positive": positive_tree_count,
            "negative": negative_tree_count,
            "zero": zero_tree_count,
        }
        dominant_direction = max(
            direction_counts,
            key=lambda direction: (direction_counts[direction], direction),
        )
        significant_tree_count = sum(row.significant for row in term_rows)
        summaries.append(
            PosteriorTreePGLSCoefficientSummaryRow(
                term=term,
                tree_fit_count=len(term_rows),
                positive_tree_count=positive_tree_count,
                negative_tree_count=negative_tree_count,
                zero_tree_count=zero_tree_count,
                dominant_direction=dominant_direction,
                direction_consistency=direction_counts[dominant_direction]
                / len(term_rows),
                significant_tree_count=significant_tree_count,
                significance_fraction=significant_tree_count / len(term_rows),
                conclusion_stability=_conclusion_stability(
                    positive_tree_count=positive_tree_count,
                    negative_tree_count=negative_tree_count,
                    zero_tree_count=zero_tree_count,
                    significant_tree_count=significant_tree_count,
                    tree_fit_count=len(term_rows),
                    significance_threshold=significance_threshold,
                ),
                mean_estimate=statistics.fmean(estimates),
                median_estimate=statistics.median(estimates),
                standard_deviation=statistics.stdev(estimates)
                if len(estimates) > 1
                else 0.0,
                minimum_estimate=min(estimates),
                maximum_estimate=max(estimates),
                lower_95_empirical_estimate=_empirical_quantile(estimates, 0.025),
                upper_95_empirical_estimate=_empirical_quantile(estimates, 0.975),
                mean_p_value=statistics.fmean(p_values),
                median_p_value=statistics.median(p_values),
                minimum_p_value=min(p_values),
                maximum_p_value=max(p_values),
            )
        )
    return summaries


def _conclusion_stability(
    *,
    positive_tree_count: int,
    negative_tree_count: int,
    zero_tree_count: int,
    significant_tree_count: int,
    tree_fit_count: int,
    significance_threshold: float,
) -> str:
    del zero_tree_count, significance_threshold
    if positive_tree_count > 0 and negative_tree_count > 0:
        return "direction_conflict"
    if positive_tree_count == 0 and negative_tree_count == 0:
        return "zero_effect"
    if significant_tree_count == tree_fit_count:
        return "stable_supported"
    if significant_tree_count == 0:
        return "stable_unsupported"
    return "mixed_support"


def _estimate_direction(estimate: float) -> str:
    if math.isclose(estimate, 0.0, abs_tol=1e-12):
        return "zero"
    return "positive" if estimate > 0.0 else "negative"


def _empirical_quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * probability
    lower_index = math.floor(position)
    upper_index = math.ceil(position)
    if lower_index == upper_index:
        return ordered[lower_index]
    lower_value = ordered[lower_index]
    upper_value = ordered[upper_index]
    fraction = position - lower_index
    return lower_value + ((upper_value - lower_value) * fraction)
