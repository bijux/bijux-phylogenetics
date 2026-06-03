from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from statistics import median

from bijux_phylogenetics.datasets.study_inputs import (
    load_taxon_table,
    validate_traits_table,
    write_taxon_rows,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import TreeNode
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

_CONTINUOUS_EXCEPTIONAL_THRESHOLD = 0.65
_CATEGORICAL_EXCEPTIONAL_THRESHOLD = 0.33


@dataclass(slots=True)
class CladeTraitStateCount:
    """One categorical state count preserved inside a clade summary row."""

    state: str
    count: int
    fraction: float


@dataclass(slots=True)
class CladeTraitExclusion:
    """One taxon excluded before clade trait summarization."""

    taxon: str
    reason: str


@dataclass(slots=True)
class CladeTraitRow:
    """One internal clade summarized for one analyzed trait."""

    clade_id: str
    node_label: str | None
    trait_kind: str
    taxon_count: int
    taxa: list[str]
    coverage_fraction: float
    mean: float | None
    median: float | None
    minimum: float | None
    maximum: float | None
    range_width: float | None
    mean_delta_from_global: float | None
    dominant_state: str | None
    dominant_state_count: int | None
    dominant_state_fraction: float | None
    dominant_state_enrichment: float | None
    distinct_state_count: int | None
    state_counts: list[CladeTraitStateCount]
    distribution_shift: float | None
    exceptionality_score: float
    exceptional: bool
    rank: int


@dataclass(slots=True)
class CladeTraitSummaryReport:
    """Clade-level trait summary report for one analyzed trait."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    trait_kind: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    excluded_taxa: list[CladeTraitExclusion]
    minimum_clade_size: int
    clade_rows: list[CladeTraitRow]
    exceptional_clades: list[str]
    top_exceptional_clade: str | None
    top_exceptionality_score: float | None
    baseline_mean: float | None
    baseline_median: float | None
    baseline_minimum: float | None
    baseline_maximum: float | None
    baseline_range_width: float | None
    baseline_dominant_state: str | None
    baseline_dominant_state_fraction: float | None
    assumptions: list[str]
    warnings: list[str]

    @property
    def analyzed_taxon_count(self) -> int:
        return len(self.analyzed_taxa)


@dataclass(slots=True)
class _ContinuousBaseline:
    mean: float
    median: float
    minimum: float
    maximum: float
    range_width: float
    sample_standard_deviation: float


@dataclass(slots=True)
class _CategoricalBaseline:
    dominant_state: str
    dominant_state_fraction: float
    state_fractions: dict[str, float]


def summarize_clade_traits(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    minimum_clade_size: int = 2,
    trait_kind: str = "auto",
) -> CladeTraitSummaryReport:
    """Summarize one continuous or categorical trait across internal clades."""
    if minimum_clade_size < 2:
        raise ComparativeMethodError(
            "clade trait summaries require a minimum clade size of at least two taxa"
        )
    if trait_kind not in {"auto", "continuous", "categorical"}:
        raise ComparativeMethodError(
            "trait kind must be 'auto', 'continuous', or 'categorical'"
        )

    tree = load_tree(tree_path)
    if len(tree.root.children) != 2:
        raise ComparativeMethodError("clade trait summaries require a rooted tree")
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    validation = validate_traits_table(traits_path, taxon_column=taxon_column)
    column_summary = next(
        (column for column in validation.trait_columns if column.name == trait),
        None,
    )
    if column_summary is None:
        raise ComparativeMethodError(f"trait table does not contain column '{trait}'")
    if column_summary.kind == "empty":
        raise ComparativeMethodError(
            f"trait column '{trait}' has no observed values for clade summarization"
        )
    resolved_trait_kind = _resolve_trait_kind(
        requested_kind=trait_kind,
        inferred_kind=column_summary.kind,
    )

    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    analyzed_taxa: list[str] = []
    continuous_values: dict[str, float] = {}
    categorical_values: dict[str, str] = {}
    exclusions: list[CladeTraitExclusion] = []
    for taxon in tree.tip_names:
        row = rows_by_taxon.get(taxon)
        if row is None:
            exclusions.append(
                CladeTraitExclusion(taxon=taxon, reason="missing_from_trait_table")
            )
            continue
        raw_value = row.get(trait, "").strip()
        if not raw_value:
            exclusions.append(
                CladeTraitExclusion(taxon=taxon, reason="missing_trait_value")
            )
            continue
        if resolved_trait_kind == "continuous":
            try:
                parsed = float(raw_value)
            except ValueError:
                exclusions.append(
                    CladeTraitExclusion(taxon=taxon, reason="non_numeric_trait_value")
                )
                continue
            continuous_values[taxon] = parsed
        else:
            categorical_values[taxon] = raw_value
        analyzed_taxa.append(taxon)

    tree_taxa = set(tree.tip_names)
    for taxon in sorted(set(table.taxa) - tree_taxa):
        exclusions.append(CladeTraitExclusion(taxon=taxon, reason="absent_from_tree"))

    if len(analyzed_taxa) < 3:
        raise ComparativeMethodError(
            "clade trait summaries require at least three analyzed taxa"
        )

    pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, analyzed_taxa)
    if len(pruned_tree.root.children) != 2:
        raise ComparativeMethodError(
            "clade trait summaries require a rooted analyzed tree"
        )

    continuous_baseline = (
        _continuous_baseline(list(continuous_values.values()))
        if resolved_trait_kind == "continuous"
        else None
    )
    categorical_baseline = (
        _categorical_baseline(list(categorical_values.values()))
        if resolved_trait_kind == "categorical"
        else None
    )
    clade_rows = _collect_clade_rows(
        pruned_tree.root,
        trait_kind=resolved_trait_kind,
        continuous_values=continuous_values,
        continuous_baseline=continuous_baseline,
        categorical_values=categorical_values,
        categorical_baseline=categorical_baseline,
        total_taxa=len(analyzed_taxa),
        minimum_clade_size=minimum_clade_size,
    )
    if not clade_rows:
        raise ComparativeMethodError(
            "no internal clades met the requested minimum clade size"
        )

    exceptional_clades = [row.clade_id for row in clade_rows if row.exceptional]
    top_exceptional_clade = clade_rows[0].clade_id if clade_rows else None
    top_exceptionality_score = (
        clade_rows[0].exceptionality_score if clade_rows else None
    )
    warnings: list[str] = []
    exclusion_reasons = {row.reason for row in exclusions}
    if "missing_from_trait_table" in exclusion_reasons:
        warnings.append(
            "trait table is missing one or more tree taxa and those taxa were pruned"
        )
    if "absent_from_tree" in exclusion_reasons:
        warnings.append("trait table contains taxa absent from the tree")
    if "missing_trait_value" in exclusion_reasons:
        warnings.append(
            "one or more overlapping taxa have missing trait values and were pruned"
        )
    if "non_numeric_trait_value" in exclusion_reasons:
        warnings.append(
            "one or more overlapping taxa have non-numeric trait values and were pruned"
        )
    if exceptional_clades:
        warnings.append(
            "one or more internal clades show elevated trait exceptionality relative to the analyzed global distribution"
        )

    assumptions = [
        "internal non-root clades are summarized after pruning the tree to analyzed taxa",
        (
            "continuous clade exceptionality ranks weighted standardized mean shifts from the analyzed global mean"
            if resolved_trait_kind == "continuous"
            else "categorical clade exceptionality ranks weighted total-variation shifts from the analyzed global state distribution"
        ),
    ]
    return CladeTraitSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=table.taxon_column,
        trait_kind=resolved_trait_kind,
        tree_taxon_count=len(tree.tip_names),
        analyzed_taxa=analyzed_taxa,
        excluded_taxa=exclusions,
        minimum_clade_size=minimum_clade_size,
        clade_rows=clade_rows,
        exceptional_clades=exceptional_clades,
        top_exceptional_clade=top_exceptional_clade,
        top_exceptionality_score=top_exceptionality_score,
        baseline_mean=None if continuous_baseline is None else continuous_baseline.mean,
        baseline_median=(
            None if continuous_baseline is None else continuous_baseline.median
        ),
        baseline_minimum=(
            None if continuous_baseline is None else continuous_baseline.minimum
        ),
        baseline_maximum=(
            None if continuous_baseline is None else continuous_baseline.maximum
        ),
        baseline_range_width=(
            None if continuous_baseline is None else continuous_baseline.range_width
        ),
        baseline_dominant_state=(
            None
            if categorical_baseline is None
            else categorical_baseline.dominant_state
        ),
        baseline_dominant_state_fraction=(
            None
            if categorical_baseline is None
            else categorical_baseline.dominant_state_fraction
        ),
        assumptions=assumptions,
        warnings=warnings,
    )


def write_clade_trait_summary_table(
    path: Path,
    report: CladeTraitSummaryReport,
) -> Path:
    """Write one summary ledger for a clade trait report."""
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "trait_kind",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "minimum_clade_size",
            "clade_count",
            "exceptional_clade_count",
            "top_exceptional_clade",
            "top_exceptionality_score",
            "baseline_mean",
            "baseline_median",
            "baseline_minimum",
            "baseline_maximum",
            "baseline_range_width",
            "baseline_dominant_state",
            "baseline_dominant_state_fraction",
            "assumptions",
            "warnings",
        ],
        rows=[
            {
                "trait": report.trait,
                "taxon_column": report.taxon_column,
                "trait_kind": report.trait_kind,
                "tree_taxon_count": str(report.tree_taxon_count),
                "analyzed_taxon_count": str(report.analyzed_taxon_count),
                "excluded_taxon_count": str(len(report.excluded_taxa)),
                "minimum_clade_size": str(report.minimum_clade_size),
                "clade_count": str(len(report.clade_rows)),
                "exceptional_clade_count": str(len(report.exceptional_clades)),
                "top_exceptional_clade": report.top_exceptional_clade or "",
                "top_exceptionality_score": _format_optional_float(
                    report.top_exceptionality_score
                ),
                "baseline_mean": _format_optional_float(report.baseline_mean),
                "baseline_median": _format_optional_float(report.baseline_median),
                "baseline_minimum": _format_optional_float(report.baseline_minimum),
                "baseline_maximum": _format_optional_float(report.baseline_maximum),
                "baseline_range_width": _format_optional_float(
                    report.baseline_range_width
                ),
                "baseline_dominant_state": report.baseline_dominant_state or "",
                "baseline_dominant_state_fraction": _format_optional_float(
                    report.baseline_dominant_state_fraction
                ),
                "assumptions": " | ".join(report.assumptions),
                "warnings": " | ".join(report.warnings),
            }
        ],
    )


def write_clade_trait_clade_table(
    path: Path,
    report: CladeTraitSummaryReport,
) -> Path:
    """Write one internal-clade trait ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "clade_id",
            "node_label",
            "trait_kind",
            "taxon_count",
            "taxa",
            "coverage_fraction",
            "mean",
            "median",
            "minimum",
            "maximum",
            "range_width",
            "mean_delta_from_global",
            "dominant_state",
            "dominant_state_count",
            "dominant_state_fraction",
            "dominant_state_enrichment",
            "distinct_state_count",
            "state_counts",
            "distribution_shift",
            "exceptionality_score",
            "exceptional",
            "rank",
        ],
        rows=[
            {
                "clade_id": row.clade_id,
                "node_label": row.node_label or "",
                "trait_kind": row.trait_kind,
                "taxon_count": str(row.taxon_count),
                "taxa": ",".join(row.taxa),
                "coverage_fraction": format(row.coverage_fraction, ".15g"),
                "mean": _format_optional_float(row.mean),
                "median": _format_optional_float(row.median),
                "minimum": _format_optional_float(row.minimum),
                "maximum": _format_optional_float(row.maximum),
                "range_width": _format_optional_float(row.range_width),
                "mean_delta_from_global": _format_optional_float(
                    row.mean_delta_from_global
                ),
                "dominant_state": row.dominant_state or "",
                "dominant_state_count": (
                    ""
                    if row.dominant_state_count is None
                    else str(row.dominant_state_count)
                ),
                "dominant_state_fraction": _format_optional_float(
                    row.dominant_state_fraction
                ),
                "dominant_state_enrichment": _format_optional_float(
                    row.dominant_state_enrichment
                ),
                "distinct_state_count": (
                    ""
                    if row.distinct_state_count is None
                    else str(row.distinct_state_count)
                ),
                "state_counts": _format_state_counts(row.state_counts),
                "distribution_shift": _format_optional_float(row.distribution_shift),
                "exceptionality_score": format(row.exceptionality_score, ".15g"),
                "exceptional": str(row.exceptional).lower(),
                "rank": str(row.rank),
            }
            for row in report.clade_rows
        ],
    )


def write_clade_trait_exclusion_table(
    path: Path,
    report: CladeTraitSummaryReport,
) -> Path:
    """Write one excluded-taxa ledger for clade trait summarization."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {"taxon": row.taxon, "reason": row.reason} for row in report.excluded_taxa
        ],
    )


def _collect_clade_rows(
    root: TreeNode,
    *,
    trait_kind: str,
    continuous_values: dict[str, float],
    continuous_baseline: _ContinuousBaseline | None,
    categorical_values: dict[str, str],
    categorical_baseline: _CategoricalBaseline | None,
    total_taxa: int,
    minimum_clade_size: int,
) -> list[CladeTraitRow]:
    rows: list[CladeTraitRow] = []

    def visit(node: TreeNode, *, is_root: bool) -> list[str]:
        if node.is_leaf():
            return [node.name] if node.name is not None else []
        taxa: list[str] = []
        for child in node.children:
            taxa.extend(visit(child, is_root=False))
        ordered_taxa = sorted(taxa)
        if not is_root and len(ordered_taxa) >= minimum_clade_size:
            if trait_kind == "continuous":
                if continuous_baseline is None:
                    raise RuntimeError(
                        "continuous clade trait traversal requires a continuous baseline"
                    )
                rows.append(
                    _continuous_clade_row(
                        clade_id="|".join(ordered_taxa),
                        node_label=node.name,
                        taxa=ordered_taxa,
                        values_by_taxon=continuous_values,
                        baseline=continuous_baseline,
                        total_taxa=total_taxa,
                    )
                )
            else:
                if categorical_baseline is None:
                    raise RuntimeError(
                        "categorical clade trait traversal requires a categorical baseline"
                    )
                rows.append(
                    _categorical_clade_row(
                        clade_id="|".join(ordered_taxa),
                        node_label=node.name,
                        taxa=ordered_taxa,
                        values_by_taxon=categorical_values,
                        baseline=categorical_baseline,
                        total_taxa=total_taxa,
                    )
                )
        return ordered_taxa

    visit(root, is_root=True)
    rows.sort(
        key=lambda row: (-row.exceptionality_score, -row.taxon_count, row.clade_id)
    )
    for rank, row in enumerate(rows, start=1):
        row.rank = rank
    return rows


def _resolve_trait_kind(*, requested_kind: str, inferred_kind: str) -> str:
    if requested_kind != "auto":
        return requested_kind
    return "continuous" if inferred_kind == "numeric" else "categorical"


def _continuous_baseline(values: list[float]) -> _ContinuousBaseline:
    mean_value = sum(values) / len(values)
    return _ContinuousBaseline(
        mean=mean_value,
        median=float(median(values)),
        minimum=min(values),
        maximum=max(values),
        range_width=max(values) - min(values),
        sample_standard_deviation=(
            math.sqrt(
                sum((value - mean_value) ** 2 for value in values) / (len(values) - 1)
            )
            if len(values) > 1
            else 0.0
        ),
    )


def _categorical_baseline(values: list[str]) -> _CategoricalBaseline:
    counts = _state_counts(values)
    dominant_state = counts[0].state
    dominant_state_fraction = counts[0].fraction
    return _CategoricalBaseline(
        dominant_state=dominant_state,
        dominant_state_fraction=dominant_state_fraction,
        state_fractions={count.state: count.fraction for count in counts},
    )


def _continuous_clade_row(
    *,
    clade_id: str,
    node_label: str | None,
    taxa: list[str],
    values_by_taxon: dict[str, float],
    baseline: _ContinuousBaseline,
    total_taxa: int,
) -> CladeTraitRow:
    values = [values_by_taxon[taxon] for taxon in taxa]
    mean_value = sum(values) / len(values)
    mean_delta = mean_value - baseline.mean
    score = 0.0
    if baseline.sample_standard_deviation > 0.0:
        score = (
            abs(mean_delta)
            / baseline.sample_standard_deviation
            * math.sqrt(len(values) / total_taxa)
        )
    return CladeTraitRow(
        clade_id=clade_id,
        node_label=node_label,
        trait_kind="continuous",
        taxon_count=len(taxa),
        taxa=taxa,
        coverage_fraction=len(taxa) / total_taxa,
        mean=mean_value,
        median=float(median(values)),
        minimum=min(values),
        maximum=max(values),
        range_width=max(values) - min(values),
        mean_delta_from_global=mean_delta,
        dominant_state=None,
        dominant_state_count=None,
        dominant_state_fraction=None,
        dominant_state_enrichment=None,
        distinct_state_count=None,
        state_counts=[],
        distribution_shift=None,
        exceptionality_score=score,
        exceptional=score >= _CONTINUOUS_EXCEPTIONAL_THRESHOLD,
        rank=0,
    )


def _categorical_clade_row(
    *,
    clade_id: str,
    node_label: str | None,
    taxa: list[str],
    values_by_taxon: dict[str, str],
    baseline: _CategoricalBaseline,
    total_taxa: int,
) -> CladeTraitRow:
    values = [values_by_taxon[taxon] for taxon in taxa]
    counts = _state_counts(values)
    dominant = counts[0]
    state_fractions = {count.state: count.fraction for count in counts}
    distribution_shift = _total_variation_shift(
        state_fractions,
        baseline.state_fractions,
    )
    score = distribution_shift * math.sqrt(len(values) / total_taxa)
    global_dominant_fraction = baseline.state_fractions.get(dominant.state, 0.0)
    return CladeTraitRow(
        clade_id=clade_id,
        node_label=node_label,
        trait_kind="categorical",
        taxon_count=len(taxa),
        taxa=taxa,
        coverage_fraction=len(taxa) / total_taxa,
        mean=None,
        median=None,
        minimum=None,
        maximum=None,
        range_width=None,
        mean_delta_from_global=None,
        dominant_state=dominant.state,
        dominant_state_count=dominant.count,
        dominant_state_fraction=dominant.fraction,
        dominant_state_enrichment=dominant.fraction - global_dominant_fraction,
        distinct_state_count=len(counts),
        state_counts=counts,
        distribution_shift=distribution_shift,
        exceptionality_score=score,
        exceptional=score >= _CATEGORICAL_EXCEPTIONAL_THRESHOLD,
        rank=0,
    )


def _state_counts(values: list[str]) -> list[CladeTraitStateCount]:
    counts: dict[str, int] = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    total = len(values)
    return [
        CladeTraitStateCount(
            state=state,
            count=count,
            fraction=count / total,
        )
        for state, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    ]


def _total_variation_shift(
    clade_fractions: dict[str, float],
    baseline_fractions: dict[str, float],
) -> float:
    states = sorted(set(clade_fractions) | set(baseline_fractions))
    return 0.5 * sum(
        abs(clade_fractions.get(state, 0.0) - baseline_fractions.get(state, 0.0))
        for state in states
    )


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".15g")


def _format_state_counts(counts: list[CladeTraitStateCount]) -> str:
    return ";".join(f"{row.state}={row.count}" for row in counts)
