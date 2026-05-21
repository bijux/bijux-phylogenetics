from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import random

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


@dataclass(slots=True)
class PhylogeneticAnovaExclusion:
    """One taxon excluded before phylogenetic ANOVA."""

    taxon: str
    reason: str
    details: str


@dataclass(slots=True)
class PhylogeneticAnovaGroupSummaryRow:
    """One categorical-group summary in the analyzed trait table."""

    group: str
    taxon_count: int
    taxa: list[str]
    mean: float
    variance: float
    minimum: float
    maximum: float


@dataclass(slots=True)
class PhylogeneticAnovaPairwiseRow:
    """One pairwise posthoc comparison between two categorical groups."""

    left_group: str
    right_group: str
    left_taxon_count: int
    right_taxon_count: int
    observed_t_statistic: float
    uncorrected_p_value: float
    adjusted_p_value: float


@dataclass(slots=True)
class PhylogeneticAnovaNullRow:
    """One observed-or-simulated F statistic in the phylogenetic null distribution."""

    simulation_index: int
    f_statistic: float
    at_or_above_observed: bool


@dataclass(slots=True)
class PhylogeneticAnovaReport:
    """Reviewer-facing summary for one phylogenetically corrected group comparison."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    response: str
    group: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    analyzed_taxon_count: int
    excluded_taxa: list[PhylogeneticAnovaExclusion]
    group_count: int
    simulation_count: int
    seed: int
    pairwise_adjustment_method: str
    brownian_sigma_squared: float
    sum_of_squares_between: float
    sum_of_squares_within: float
    mean_square_between: float
    mean_square_within: float
    f_statistic: float
    p_value: float
    low_sample_group_count: int
    group_rows: list[PhylogeneticAnovaGroupSummaryRow]
    pairwise_rows: list[PhylogeneticAnovaPairwiseRow]
    null_rows: list[PhylogeneticAnovaNullRow]
    warnings: list[str]


def summarize_phylogenetic_anova(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str,
    group: str,
    taxon_column: str | None = None,
    simulations: int = 1000,
    seed: int = 1,
    pairwise_adjustment_method: str = "holm",
) -> PhylogeneticAnovaReport:
    """Run a simulation-based phylogenetic ANOVA over one continuous trait and one group."""
    if simulations < 2:
        raise ComparativeMethodError(
            "phylogenetic ANOVA requires at least two simulations including the observed statistic"
        )
    if pairwise_adjustment_method != "holm":
        raise ComparativeMethodError(
            "supported phylogenetic ANOVA pairwise adjustment is 'holm'"
        )
    tree = load_tree(tree_path)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    if response not in table.columns:
        raise ComparativeMethodError(
            f"trait table does not contain response column '{response}'"
        )
    if group not in table.columns:
        raise ComparativeMethodError(
            f"trait table does not contain group column '{group}'"
        )
    if len(tree.root.children) != 2:
        raise ComparativeMethodError("phylogenetic ANOVA requires a rooted tree")
    if any(
        node is not tree.root and node.branch_length is None
        for node in tree.iter_nodes()
    ):
        raise ComparativeMethodError(
            "phylogenetic ANOVA requires complete branch lengths"
        )

    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    analyzed_taxa_tree_order: list[str] = []
    analyzed_values_by_taxon: dict[str, float] = {}
    analyzed_groups_by_taxon: dict[str, str] = {}
    exclusions: list[PhylogeneticAnovaExclusion] = []
    for taxon in tree.tip_names:
        row = rows_by_taxon.get(taxon)
        if row is None:
            exclusions.append(
                PhylogeneticAnovaExclusion(
                    taxon=taxon,
                    reason="missing_from_trait_table",
                    details="taxon is present in the tree but absent from the trait table",
                )
            )
            continue
        group_value = row.get(group, "")
        response_value = row.get(response, "")
        missing_columns = [
            column
            for column, value in ((response, response_value), (group, group_value))
            if not value
        ]
        if missing_columns:
            exclusions.append(
                PhylogeneticAnovaExclusion(
                    taxon=taxon,
                    reason="missing_value",
                    details=f"taxon is missing required value(s): {', '.join(missing_columns)}",
                )
            )
            continue
        try:
            numeric_response = float(response_value)
        except ValueError as error:
            raise ComparativeMethodError(
                f"response column '{response}' must be numeric for phylogenetic ANOVA"
            ) from error
        analyzed_taxa_tree_order.append(taxon)
        analyzed_values_by_taxon[taxon] = numeric_response
        analyzed_groups_by_taxon[taxon] = group_value
    for taxon in table.taxa:
        if taxon not in tree.tip_names:
            exclusions.append(
                PhylogeneticAnovaExclusion(
                    taxon=taxon,
                    reason="absent_from_tree",
                    details="taxon is present in the trait table but absent from the tree",
                )
            )

    analyzed_taxa_input_order = [
        row[table.taxon_column]
        for row in table.rows
        if row[table.taxon_column] in analyzed_values_by_taxon
    ]
    if len(analyzed_taxa_tree_order) < 4:
        raise ComparativeMethodError(
            "phylogenetic ANOVA requires at least four taxa after pruning"
        )
    group_order = sorted(
        {analyzed_groups_by_taxon[taxon] for taxon in analyzed_taxa_tree_order}
    )
    if len(group_order) < 2:
        raise ComparativeMethodError(
            "phylogenetic ANOVA requires at least two observed groups after pruning"
        )
    taxa_by_group = {
        group_name: [
            taxon
            for taxon in analyzed_taxa_input_order
            if analyzed_groups_by_taxon[taxon] == group_name
        ]
        for group_name in group_order
    }
    group_counts = {group_name: len(taxa) for group_name, taxa in taxa_by_group.items()}
    small_groups = sorted(
        group_name for group_name, count in group_counts.items() if count < 3
    )
    singleton_groups = sorted(
        group_name for group_name, count in group_counts.items() if count < 2
    )
    if singleton_groups:
        singleton_text = ", ".join(singleton_groups)
        raise ComparativeMethodError(
            f"phylogenetic ANOVA requires at least two taxa per group after pruning; singleton group(s): {singleton_text}"
        )

    pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, analyzed_taxa_tree_order)
    brownian_sigma_squared = _brownian_rate_from_independent_contrasts(
        pruned_tree,
        analyzed_values_by_taxon,
    )
    observed_summary = _anova_summary(
        analyzed_taxa_tree_order,
        analyzed_values_by_taxon,
        analyzed_groups_by_taxon,
        group_order,
    )
    observed_pairwise = _pairwise_t_statistics(
        analyzed_values_by_taxon,
        analyzed_groups_by_taxon,
        taxa_by_group,
        group_order,
    )
    randomizer = random.Random(seed)  # nosec B311
    null_rows = [
        PhylogeneticAnovaNullRow(
            simulation_index=1,
            f_statistic=observed_summary["f_statistic"],
            at_or_above_observed=True,
        )
    ]
    pairwise_null_distributions: dict[tuple[str, str], list[float]] = {
        key: [value] for key, value in observed_pairwise.items()
    }
    exceed_count = 1
    for simulation_index in range(2, simulations + 1):
        simulated_values = _simulate_brownian_tip_values(
            pruned_tree,
            analyzed_taxa_tree_order,
            sigma_squared=brownian_sigma_squared,
            rng=randomizer,
        )
        simulated_values_by_taxon = dict(
            zip(analyzed_taxa_tree_order, simulated_values, strict=True)
        )
        simulated_summary = _anova_summary(
            analyzed_taxa_tree_order,
            simulated_values_by_taxon,
            analyzed_groups_by_taxon,
            group_order,
        )
        simulated_f = simulated_summary["f_statistic"]
        at_or_above_observed = simulated_f >= observed_summary["f_statistic"]
        if at_or_above_observed:
            exceed_count += 1
        null_rows.append(
            PhylogeneticAnovaNullRow(
                simulation_index=simulation_index,
                f_statistic=simulated_f,
                at_or_above_observed=at_or_above_observed,
            )
        )
        simulated_pairwise = _pairwise_t_statistics(
            simulated_values_by_taxon,
            analyzed_groups_by_taxon,
            taxa_by_group,
            group_order,
        )
        for key, value in simulated_pairwise.items():
            pairwise_null_distributions[key].append(value)

    pair_keys = list(observed_pairwise)
    uncorrected_pairwise_p_values = [
        sum(
            abs(value) >= abs(observed_pairwise[key])
            for value in pairwise_null_distributions[key]
        )
        / simulations
        for key in pair_keys
    ]
    adjusted_pairwise_p_values = _holm_adjustment(uncorrected_pairwise_p_values)
    pairwise_rows = [
        PhylogeneticAnovaPairwiseRow(
            left_group=left_group,
            right_group=right_group,
            left_taxon_count=group_counts[left_group],
            right_taxon_count=group_counts[right_group],
            observed_t_statistic=observed_pairwise[(left_group, right_group)],
            uncorrected_p_value=uncorrected_p_value,
            adjusted_p_value=adjusted_p_value,
        )
        for (left_group, right_group), uncorrected_p_value, adjusted_p_value in zip(
            pair_keys,
            uncorrected_pairwise_p_values,
            adjusted_pairwise_p_values,
            strict=True,
        )
    ]
    group_rows = [
        _group_summary_row(
            group_name,
            taxa_by_group[group_name],
            analyzed_values_by_taxon,
        )
        for group_name in group_order
    ]
    warnings: list[str] = []
    if any(item.reason == "missing_from_trait_table" for item in exclusions):
        warnings.append(
            "one or more tree taxa are absent from the trait table and were pruned before phylogenetic ANOVA"
        )
    if any(item.reason == "absent_from_tree" for item in exclusions):
        warnings.append("trait table contains taxa absent from the tree")
    if any(item.reason == "missing_value" for item in exclusions):
        warnings.append(
            "one or more overlapping taxa have missing response or group values and were pruned before phylogenetic ANOVA"
        )
    if small_groups:
        warnings.append(
            "one or more groups contain fewer than three taxa, so phylogenetic ANOVA power may be unstable"
        )
    return PhylogeneticAnovaReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        response=response,
        group=group,
        tree_taxon_count=tree.tip_count,
        analyzed_taxa=analyzed_taxa_input_order,
        analyzed_taxon_count=len(analyzed_taxa_input_order),
        excluded_taxa=exclusions,
        group_count=len(group_order),
        simulation_count=simulations,
        seed=seed,
        pairwise_adjustment_method=pairwise_adjustment_method,
        brownian_sigma_squared=brownian_sigma_squared,
        sum_of_squares_between=observed_summary["sum_of_squares_between"],
        sum_of_squares_within=observed_summary["sum_of_squares_within"],
        mean_square_between=observed_summary["mean_square_between"],
        mean_square_within=observed_summary["mean_square_within"],
        f_statistic=observed_summary["f_statistic"],
        p_value=exceed_count / simulations,
        low_sample_group_count=len(small_groups),
        group_rows=group_rows,
        pairwise_rows=pairwise_rows,
        null_rows=null_rows,
        warnings=list(dict.fromkeys(warnings)),
    )


def write_phylogenetic_anova_summary_table(
    path: Path,
    report: PhylogeneticAnovaReport,
) -> Path:
    """Write one summary ledger for phylogenetic ANOVA."""
    return write_taxon_rows(
        path,
        columns=[
            "response",
            "group",
            "taxon_column",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "group_count",
            "simulation_count",
            "seed",
            "pairwise_adjustment_method",
            "brownian_sigma_squared",
            "sum_of_squares_between",
            "sum_of_squares_within",
            "mean_square_between",
            "mean_square_within",
            "f_statistic",
            "p_value",
            "low_sample_group_count",
        ],
        rows=[
            {
                "response": report.response,
                "group": report.group,
                "taxon_column": report.taxon_column,
                "tree_taxon_count": str(report.tree_taxon_count),
                "analyzed_taxon_count": str(report.analyzed_taxon_count),
                "excluded_taxon_count": str(len(report.excluded_taxa)),
                "group_count": str(report.group_count),
                "simulation_count": str(report.simulation_count),
                "seed": str(report.seed),
                "pairwise_adjustment_method": report.pairwise_adjustment_method,
                "brownian_sigma_squared": format(report.brownian_sigma_squared, ".15g"),
                "sum_of_squares_between": format(report.sum_of_squares_between, ".15g"),
                "sum_of_squares_within": format(report.sum_of_squares_within, ".15g"),
                "mean_square_between": format(report.mean_square_between, ".15g"),
                "mean_square_within": format(report.mean_square_within, ".15g"),
                "f_statistic": format(report.f_statistic, ".15g"),
                "p_value": format(report.p_value, ".15g"),
                "low_sample_group_count": str(report.low_sample_group_count),
            }
        ],
    )


def write_phylogenetic_anova_group_table(
    path: Path,
    report: PhylogeneticAnovaReport,
) -> Path:
    """Write one categorical-group summary ledger for phylogenetic ANOVA."""
    return write_taxon_rows(
        path,
        columns=[
            "group",
            "taxon_count",
            "taxa",
            "mean",
            "variance",
            "minimum",
            "maximum",
        ],
        rows=[
            {
                "group": row.group,
                "taxon_count": str(row.taxon_count),
                "taxa": ",".join(row.taxa),
                "mean": format(row.mean, ".15g"),
                "variance": format(row.variance, ".15g"),
                "minimum": format(row.minimum, ".15g"),
                "maximum": format(row.maximum, ".15g"),
            }
            for row in report.group_rows
        ],
    )


def write_phylogenetic_anova_pairwise_table(
    path: Path,
    report: PhylogeneticAnovaReport,
) -> Path:
    """Write one pairwise posthoc-comparison ledger for phylogenetic ANOVA."""
    return write_taxon_rows(
        path,
        columns=[
            "left_group",
            "right_group",
            "left_taxon_count",
            "right_taxon_count",
            "observed_t_statistic",
            "uncorrected_p_value",
            "adjusted_p_value",
        ],
        rows=[
            {
                "left_group": row.left_group,
                "right_group": row.right_group,
                "left_taxon_count": str(row.left_taxon_count),
                "right_taxon_count": str(row.right_taxon_count),
                "observed_t_statistic": format(row.observed_t_statistic, ".15g"),
                "uncorrected_p_value": format(row.uncorrected_p_value, ".15g"),
                "adjusted_p_value": format(row.adjusted_p_value, ".15g"),
            }
            for row in report.pairwise_rows
        ],
    )


def write_phylogenetic_anova_simulation_table(
    path: Path,
    report: PhylogeneticAnovaReport,
) -> Path:
    """Write one observed-plus-null F-statistic ledger for phylogenetic ANOVA."""
    return write_taxon_rows(
        path,
        columns=["simulation_index", "f_statistic", "at_or_above_observed"],
        rows=[
            {
                "simulation_index": str(row.simulation_index),
                "f_statistic": format(row.f_statistic, ".15g"),
                "at_or_above_observed": (
                    "true" if row.at_or_above_observed else "false"
                ),
            }
            for row in report.null_rows
        ],
    )


def write_phylogenetic_anova_exclusion_table(
    path: Path,
    report: PhylogeneticAnovaReport,
) -> Path:
    """Write one excluded-taxa ledger for phylogenetic ANOVA."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason", "details"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
                "details": row.details,
            }
            for row in report.excluded_taxa
        ],
    )


def _group_summary_row(
    group: str,
    taxa: list[str],
    values_by_taxon: dict[str, float],
) -> PhylogeneticAnovaGroupSummaryRow:
    values = [values_by_taxon[taxon] for taxon in taxa]
    mean = sum(values) / len(values)
    variance = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
    return PhylogeneticAnovaGroupSummaryRow(
        group=group,
        taxon_count=len(taxa),
        taxa=list(taxa),
        mean=mean,
        variance=variance,
        minimum=min(values),
        maximum=max(values),
    )


def _anova_summary(
    taxa: list[str],
    values_by_taxon: dict[str, float],
    groups_by_taxon: dict[str, str],
    group_order: list[str],
) -> dict[str, float]:
    values = [values_by_taxon[taxon] for taxon in taxa]
    grand_mean = sum(values) / len(values)
    taxa_by_group = {
        group_name: [taxon for taxon in taxa if groups_by_taxon[taxon] == group_name]
        for group_name in group_order
    }
    sum_of_squares_between = sum(
        len(group_taxa)
        * (
            (sum(values_by_taxon[taxon] for taxon in group_taxa) / len(group_taxa))
            - grand_mean
        )
        ** 2
        for group_taxa in taxa_by_group.values()
    )
    sum_of_squares_within = sum(
        (
            values_by_taxon[taxon]
            - (sum(values_by_taxon[item] for item in group_taxa) / len(group_taxa))
        )
        ** 2
        for group_taxa in taxa_by_group.values()
        for taxon in group_taxa
    )
    mean_square_between = sum_of_squares_between / (len(group_order) - 1)
    mean_square_within = sum_of_squares_within / (len(taxa) - len(group_order))
    return {
        "sum_of_squares_between": sum_of_squares_between,
        "sum_of_squares_within": sum_of_squares_within,
        "mean_square_between": mean_square_between,
        "mean_square_within": mean_square_within,
        "f_statistic": mean_square_between / mean_square_within,
    }


def _pairwise_t_statistics(
    values_by_taxon: dict[str, float],
    groups_by_taxon: dict[str, str],
    taxa_by_group: dict[str, list[str]],
    group_order: list[str],
) -> dict[tuple[str, str], float]:
    means = {
        group_name: sum(values_by_taxon[taxon] for taxon in taxa) / len(taxa)
        for group_name, taxa in taxa_by_group.items()
    }
    sample_variances = {
        group_name: sum(
            (values_by_taxon[taxon] - means[group_name]) ** 2 for taxon in taxa
        )
        / (len(taxa) - 1)
        for group_name, taxa in taxa_by_group.items()
    }
    total_degrees_of_freedom = sum(len(taxa) - 1 for taxa in taxa_by_group.values())
    pooled_standard_deviation = math.sqrt(
        sum(
            sample_variances[group_name] * (len(taxa_by_group[group_name]) - 1)
            for group_name in group_order
        )
        / total_degrees_of_freedom
    )
    statistics: dict[tuple[str, str], float] = {}
    for index, left_group in enumerate(group_order):
        for right_group in group_order[index + 1 :]:
            left_taxa = taxa_by_group[left_group]
            right_taxa = taxa_by_group[right_group]
            difference = means[left_group] - means[right_group]
            standard_error = pooled_standard_deviation * math.sqrt(
                (1.0 / len(left_taxa)) + (1.0 / len(right_taxa))
            )
            statistics[(left_group, right_group)] = difference / standard_error
    return statistics


def _holm_adjustment(p_values: list[float]) -> list[float]:
    ranked = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [0.0] * len(p_values)
    running_maximum = 0.0
    family_size = len(p_values)
    for rank, (index, value) in enumerate(ranked, start=1):
        adjusted_value = min(1.0, value * (family_size - rank + 1))
        running_maximum = max(running_maximum, adjusted_value)
        adjusted[index] = running_maximum
    return adjusted


def _brownian_rate_from_independent_contrasts(
    tree: PhyloTree,
    values_by_taxon: dict[str, float],
) -> float:
    contrasts, _, _ = _compute_pic_payload(tree.root, values_by_taxon)
    if not contrasts:
        raise ComparativeMethodError(
            "phylogenetic ANOVA requires at least one independent contrast"
        )
    return sum(contrast * contrast for contrast in contrasts) / len(contrasts)


def _compute_pic_payload(
    node: TreeNode,
    values_by_taxon: dict[str, float],
) -> tuple[list[float], float, float]:
    if node.is_leaf():
        if node.name is None:
            raise ComparativeMethodError("tree contains an unnamed terminal taxon")
        if node.branch_length is None:
            raise ComparativeMethodError(
                "phylogenetic ANOVA requires complete branch lengths"
            )
        return [], values_by_taxon[node.name], node.branch_length
    if len(node.children) != 2:
        raise ComparativeMethodError(
            "phylogenetic ANOVA currently requires a strictly binary tree"
        )
    left_contrasts, left_value, left_variance = _compute_pic_payload(
        node.children[0],
        values_by_taxon,
    )
    right_contrasts, right_value, right_variance = _compute_pic_payload(
        node.children[1],
        values_by_taxon,
    )
    expected_variance = left_variance + right_variance
    contrast = (left_value - right_value) / math.sqrt(expected_variance)
    ancestral_value = (
        (left_value / left_variance) + (right_value / right_variance)
    ) / ((1.0 / left_variance) + (1.0 / right_variance))
    propagated_variance = (left_variance * right_variance) / expected_variance
    if node.branch_length is not None:
        propagated_variance += node.branch_length
    return (
        left_contrasts + right_contrasts + [contrast],
        ancestral_value,
        propagated_variance,
    )


def _simulate_brownian_tip_values(
    tree: PhyloTree,
    taxa: list[str],
    *,
    sigma_squared: float,
    rng: random.Random,
) -> list[float]:
    sigma = math.sqrt(max(sigma_squared, 0.0))
    values_by_taxon: dict[str, float] = {}

    def visit(node: TreeNode, state: float) -> None:
        current_state = state
        if node is not tree.root:
            if node.branch_length is None:
                raise ComparativeMethodError(
                    "phylogenetic ANOVA requires complete branch lengths"
                )
            current_state += rng.gauss(0.0, sigma * math.sqrt(node.branch_length))
        if node.is_leaf():
            if node.name is None:
                raise ComparativeMethodError("tree contains an unnamed terminal taxon")
            values_by_taxon[node.name] = current_state
            return
        for child in node.children:
            visit(child, current_state)

    visit(tree.root, 0.0)
    return [values_by_taxon[taxon] for taxon in taxa]
