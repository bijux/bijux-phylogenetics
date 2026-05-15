from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    AncestralContinuousDataset,
    dump_pruned_tree,
    load_continuous_dataset,
    node_descendant_taxa,
    node_signature,
    stable_value,
    write_ancestral_rows,
)
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.core.ultrametric import summarize_ultrametric_tip_depths

_NORMAL_95_CRITICAL = 1.959963984540054
_BROWNIAN_CONDITION_THRESHOLD = 1e12
_SOLVER_REGULARIZATION_EPSILON = 1e-8


@dataclass(slots=True)
class ContinuousAncestralEstimate:
    """One continuous ancestral-state estimate for a tree node."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    estimate: float
    standard_error: float
    lower_95_interval: float
    upper_95_interval: float
    uncertainty_width: float
    confidence: float
    interpretation: str
    unstable: bool
    downstream_risks: list[str]


@dataclass(slots=True)
class ContinuousAncestralBrownianFitDiagnostics:
    """Explicit Brownian covariance and solver diagnostics for one ancestral fit."""

    covariance_model: str
    tree_is_ultrametric: bool
    minimum_root_to_tip_depth: float
    maximum_root_to_tip_depth: float
    minimum_branch_length: float
    maximum_branch_length: float
    covariance_matrix_dimension: int
    covariance_matrix_rank: int
    covariance_singular: bool
    covariance_near_singular: bool
    covariance_positive_definite: bool
    covariance_condition_number: float
    covariance_log_determinant: float | None
    solver_name: str
    solver_regularized: bool
    solver_regularization_epsilon: float | None
    log_likelihood: float
    residual_sigma_squared: float


@dataclass(slots=True)
class ContinuousAncestralReport:
    """Continuous ancestral-state reconstruction report."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    alpha: float
    taxon_count: int
    analysis_tree_newick: str
    dropped_missing_taxa: list[str]
    dropped_non_numeric_taxa: list[str]
    warnings: list[str]
    unstable_nodes: list[str]
    weak_support_nodes: list[str]
    brownian_fit_diagnostics: ContinuousAncestralBrownianFitDiagnostics | None
    estimates: list[ContinuousAncestralEstimate]


@dataclass(slots=True)
class ContinuousAncestralSummary:
    """Reviewer-facing summary for one continuous ancestral reconstruction."""

    trait: str
    taxon_column: str
    model: str
    alpha: float
    analyzed_taxon_count: int
    excluded_taxon_count: int
    missing_tip_taxon_count: int
    non_numeric_tip_taxon_count: int
    internal_node_count: int
    unstable_node_count: int
    weak_support_node_count: int
    root_node: str
    root_estimate: float
    root_standard_error: float
    root_lower_95_interval: float
    root_upper_95_interval: float
    tree_is_ultrametric: bool | None
    covariance_near_singular: bool | None
    covariance_condition_number: float | None
    log_likelihood: float | None
    residual_sigma_squared: float | None
    warning_count: int


@dataclass(slots=True)
class ContinuousAncestralExclusion:
    """One excluded tip from a continuous ancestral reconstruction."""

    taxon: str
    reason: str


def reconstruct_continuous_ancestral_states(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
) -> ContinuousAncestralReport:
    """Reconstruct continuous ancestral states under a Brownian or OU-style model."""
    if model not in {"brownian", "ou"}:
        raise ValueError(f"unsupported continuous ancestral model: {model}")
    if alpha <= 0:
        raise ValueError(
            f"alpha must be positive for continuous ancestral reconstruction, got {alpha}"
        )
    dataset = load_continuous_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    brownian_fit_diagnostics = (
        _summarize_brownian_fit_diagnostics(dataset)
        if model == "brownian"
        else None
    )
    return _reconstruct_continuous_from_dataset(
        dataset,
        working_tree=dataset.tree,
        model=model,
        alpha=alpha,
        brownian_fit_diagnostics=brownian_fit_diagnostics,
    )


def _reconstruct_continuous_from_dataset(
    dataset: AncestralContinuousDataset,
    *,
    working_tree: PhyloTree,
    model: str,
    alpha: float,
    brownian_fit_diagnostics: ContinuousAncestralBrownianFitDiagnostics | None,
) -> ContinuousAncestralReport:
    global_mean = sum(dataset.values_by_taxon[taxon] for taxon in dataset.taxa) / len(
        dataset.taxa
    )
    sigma = _sample_standard_deviation(
        [dataset.values_by_taxon[taxon] for taxon in dataset.taxa]
    )
    trait_range = (
        max(dataset.values_by_taxon.values()) - min(dataset.values_by_taxon.values())
        if dataset.values_by_taxon
        else 0.0
    )
    estimates: list[ContinuousAncestralEstimate] = []

    def visit(node) -> tuple[float, float]:
        if node.is_leaf():
            estimate = dataset.values_by_taxon[node.name]
            estimates.append(
                ContinuousAncestralEstimate(
                    node=node_signature(node),
                    node_name=node.name,
                    is_tip=True,
                    descendant_taxa=node_descendant_taxa(node),
                    estimate=stable_value(estimate),
                    standard_error=0.0,
                    lower_95_interval=stable_value(estimate),
                    upper_95_interval=stable_value(estimate),
                    uncertainty_width=0.0,
                    confidence=1.0,
                    interpretation="observed tip value",
                    unstable=False,
                    downstream_risks=[],
                )
            )
            return estimate, float(node.branch_length or 0.0)

        if len(node.children) != 2:
            raise ValueError(
                "continuous ancestral reconstruction requires a fully dichotomous rooted tree"
            )

        left_child, right_child = node.children
        left_estimate, left_working_length = visit(left_child)
        right_estimate, right_working_length = visit(right_child)

        if model == "brownian":
            sum_working_lengths = max(
                left_working_length + right_working_length,
                1e-12,
            )
            estimate = (
                left_estimate * right_working_length
                + right_estimate * left_working_length
            ) / sum_working_lengths
            variance = sum_working_lengths
            returned_length = (
                float(node.branch_length or 0.0)
                + (left_working_length * right_working_length) / sum_working_lengths
            )
        else:
            child_payloads: list[tuple[float, float]] = []
            for child, child_estimate, child_variance in (
                (left_child, left_estimate, left_working_length),
                (right_child, right_estimate, right_working_length),
            ):
                branch_length = float(child.branch_length or 0.0)
                shrink = math.exp(-alpha * branch_length)
                transformed_estimate = (
                    shrink * child_estimate + (1.0 - shrink) * global_mean
                )
                stationary_variance = ((sigma**2) / (2.0 * alpha)) * (
                    1.0 - math.exp(-2.0 * alpha * branch_length)
                )
                propagated_variance = (
                    child_variance * math.exp(-2.0 * alpha * branch_length)
                ) + stationary_variance
                child_payloads.append(
                    (transformed_estimate, max(propagated_variance, 1e-12))
                )
            weight_sum = sum(
                1.0 / child_variance for _, child_variance in child_payloads
            )
            estimate = (
                sum(
                    (value / child_variance) for value, child_variance in child_payloads
                )
                / weight_sum
            )
            variance = 1.0 / weight_sum
            returned_length = variance

        standard_error = math.sqrt(max(variance, 0.0))
        lower = estimate - _NORMAL_95_CRITICAL * standard_error
        upper = estimate + _NORMAL_95_CRITICAL * standard_error
        uncertainty_width = max(0.0, upper - lower)
        confidence, unstable = _continuous_confidence(uncertainty_width, trait_range)
        interpretation = _continuous_interpretation(
            uncertainty_width, trait_range, unstable=unstable
        )
        estimates.append(
            ContinuousAncestralEstimate(
                node=node_signature(node),
                node_name=node.name,
                is_tip=False,
                descendant_taxa=node_descendant_taxa(node),
                estimate=stable_value(estimate),
                standard_error=stable_value(standard_error),
                lower_95_interval=stable_value(lower),
                upper_95_interval=stable_value(upper),
                uncertainty_width=stable_value(uncertainty_width),
                confidence=stable_value(confidence),
                interpretation=interpretation,
                unstable=unstable,
                downstream_risks=_continuous_downstream_risks(unstable),
            )
        )
        return estimate, returned_length

    visit(working_tree.root)
    ordered_estimates = _ordered_estimates(dataset, estimates)
    unstable_nodes = [
        estimate.node
        for estimate in ordered_estimates
        if not estimate.is_tip and estimate.unstable
    ]
    weak_support_nodes = [
        estimate.node
        for estimate in ordered_estimates
        if not estimate.is_tip and estimate.confidence < 0.75
    ]
    warnings = list(dataset.warnings)
    if unstable_nodes:
        warnings.append(
            "one or more continuous ancestral estimates have broad uncertainty intervals"
        )
    if weak_support_nodes:
        warnings.append(
            "low-confidence ancestral estimates should not be overinterpreted for evolutionary timing or trait polarity"
        )
    if (
        brownian_fit_diagnostics is not None
        and brownian_fit_diagnostics.covariance_near_singular
    ):
        warnings.append(
            "Brownian covariance diagnostics indicate a singular or ill-conditioned fit, so ancestral uncertainty should be interpreted cautiously"
        )
    if (
        brownian_fit_diagnostics is not None
        and brownian_fit_diagnostics.solver_regularized
    ):
        warnings.append(
            "Brownian covariance inversion required light diagonal regularization for numerical stability"
        )
    return ContinuousAncestralReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        model=model,
        alpha=stable_value(alpha),
        taxon_count=len(dataset.taxa),
        analysis_tree_newick=dump_pruned_tree(working_tree),
        dropped_missing_taxa=dataset.dropped_missing_taxa,
        dropped_non_numeric_taxa=dataset.dropped_non_numeric_taxa,
        warnings=warnings,
        unstable_nodes=unstable_nodes,
        weak_support_nodes=weak_support_nodes,
        brownian_fit_diagnostics=brownian_fit_diagnostics,
        estimates=ordered_estimates,
    )


def summarize_continuous_ancestral_report(
    report: ContinuousAncestralReport,
) -> ContinuousAncestralSummary:
    """Summarize the main review facts for one continuous ancestral report."""
    internal_estimates = [
        estimate for estimate in report.estimates if not estimate.is_tip
    ]
    if not internal_estimates:
        raise ValueError(
            "continuous ancestral summary requires at least one internal-node estimate"
        )
    root_estimate = max(
        internal_estimates,
        key=lambda estimate: (
            len(estimate.descendant_taxa),
            estimate.node,
        ),
    )
    return ContinuousAncestralSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        alpha=report.alpha,
        analyzed_taxon_count=report.taxon_count,
        excluded_taxon_count=len(report.dropped_missing_taxa)
        + len(report.dropped_non_numeric_taxa),
        missing_tip_taxon_count=len(report.dropped_missing_taxa),
        non_numeric_tip_taxon_count=len(report.dropped_non_numeric_taxa),
        internal_node_count=len(internal_estimates),
        unstable_node_count=len(report.unstable_nodes),
        weak_support_node_count=len(report.weak_support_nodes),
        root_node=root_estimate.node,
        root_estimate=root_estimate.estimate,
        root_standard_error=root_estimate.standard_error,
        root_lower_95_interval=root_estimate.lower_95_interval,
        root_upper_95_interval=root_estimate.upper_95_interval,
        tree_is_ultrametric=(
            None
            if report.brownian_fit_diagnostics is None
            else report.brownian_fit_diagnostics.tree_is_ultrametric
        ),
        covariance_near_singular=(
            None
            if report.brownian_fit_diagnostics is None
            else report.brownian_fit_diagnostics.covariance_near_singular
        ),
        covariance_condition_number=(
            None
            if report.brownian_fit_diagnostics is None
            else report.brownian_fit_diagnostics.covariance_condition_number
        ),
        log_likelihood=(
            None
            if report.brownian_fit_diagnostics is None
            else report.brownian_fit_diagnostics.log_likelihood
        ),
        residual_sigma_squared=(
            None
            if report.brownian_fit_diagnostics is None
            else report.brownian_fit_diagnostics.residual_sigma_squared
        ),
        warning_count=len(report.warnings),
    )


def continuous_ancestral_exclusions(
    report: ContinuousAncestralReport,
) -> list[ContinuousAncestralExclusion]:
    """Return one explicit exclusion row per dropped tip taxon."""
    rows = [
        ContinuousAncestralExclusion(
            taxon=taxon,
            reason="missing_trait_value",
        )
        for taxon in report.dropped_missing_taxa
    ]
    rows.extend(
        ContinuousAncestralExclusion(
            taxon=taxon,
            reason="non_numeric_trait_value",
        )
        for taxon in report.dropped_non_numeric_taxa
    )
    return rows


def write_continuous_ancestral_summary_table(
    path: Path, report: ContinuousAncestralReport
) -> Path:
    """Write one summary ledger for a continuous ancestral reconstruction."""
    summary = summarize_continuous_ancestral_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "alpha",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "missing_tip_taxon_count",
            "non_numeric_tip_taxon_count",
            "internal_node_count",
            "unstable_node_count",
            "weak_support_node_count",
            "root_node",
            "root_estimate",
            "root_standard_error",
            "root_lower_95_interval",
            "root_upper_95_interval",
            "tree_is_ultrametric",
            "covariance_near_singular",
            "covariance_condition_number",
            "log_likelihood",
            "residual_sigma_squared",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "alpha": str(summary.alpha),
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "missing_tip_taxon_count": str(summary.missing_tip_taxon_count),
                "non_numeric_tip_taxon_count": str(summary.non_numeric_tip_taxon_count),
                "internal_node_count": str(summary.internal_node_count),
                "unstable_node_count": str(summary.unstable_node_count),
                "weak_support_node_count": str(summary.weak_support_node_count),
                "root_node": summary.root_node,
                "root_estimate": str(summary.root_estimate),
                "root_standard_error": str(summary.root_standard_error),
                "root_lower_95_interval": str(summary.root_lower_95_interval),
                "root_upper_95_interval": str(summary.root_upper_95_interval),
                "tree_is_ultrametric": (
                    ""
                    if summary.tree_is_ultrametric is None
                    else str(summary.tree_is_ultrametric).lower()
                ),
                "covariance_near_singular": (
                    ""
                    if summary.covariance_near_singular is None
                    else str(summary.covariance_near_singular).lower()
                ),
                "covariance_condition_number": (
                    ""
                    if summary.covariance_condition_number is None
                    else str(summary.covariance_condition_number)
                ),
                "log_likelihood": (
                    "" if summary.log_likelihood is None else str(summary.log_likelihood)
                ),
                "residual_sigma_squared": (
                    ""
                    if summary.residual_sigma_squared is None
                    else str(summary.residual_sigma_squared)
                ),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_continuous_ancestral_uncertainty_table(
    path: Path, report: ContinuousAncestralReport
) -> Path:
    """Write one internal-node uncertainty ledger for a continuous reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "estimate",
            "standard_error",
            "lower_95_interval",
            "upper_95_interval",
            "uncertainty_width",
            "confidence",
            "interpretation",
            "unstable",
        ],
        rows=[
            {
                "node": estimate.node,
                "node_name": estimate.node_name or "",
                "descendant_taxa": ",".join(estimate.descendant_taxa),
                "estimate": str(estimate.estimate),
                "standard_error": str(estimate.standard_error),
                "lower_95_interval": str(estimate.lower_95_interval),
                "upper_95_interval": str(estimate.upper_95_interval),
                "uncertainty_width": str(estimate.uncertainty_width),
                "confidence": str(estimate.confidence),
                "interpretation": estimate.interpretation,
                "unstable": str(estimate.unstable).lower(),
            }
            for estimate in report.estimates
            if not estimate.is_tip
        ],
    )


def write_continuous_ancestral_exclusion_table(
    path: Path, report: ContinuousAncestralReport
) -> Path:
    """Write one explicit excluded-tip ledger for a continuous reconstruction."""
    exclusions = continuous_ancestral_exclusions(report)
    return write_ancestral_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
            }
            for row in exclusions
        ],
    )


def _sample_standard_deviation(values: list[float]) -> float:
    if len(values) < 2:
        return 1.0
    mean_value = sum(values) / len(values)
    variance = sum((value - mean_value) ** 2 for value in values) / (len(values) - 1)
    return math.sqrt(max(variance, 1e-12))


def _ordered_estimates(
    dataset: AncestralContinuousDataset,
    estimates: list[ContinuousAncestralEstimate],
) -> list[ContinuousAncestralEstimate]:
    node_order = {
        node_signature(node): index
        for index, node in enumerate(dataset.tree.iter_nodes())
    }
    return sorted(estimates, key=lambda estimate: node_order[estimate.node])


def _continuous_confidence(
    uncertainty_width: float, trait_range: float
) -> tuple[float, bool]:
    if uncertainty_width == 0.0:
        return 1.0, False
    scale = max(trait_range, 1e-12)
    relative_width = uncertainty_width / scale
    confidence = max(0.0, min(1.0, 1.0 - min(relative_width, 1.0)))
    return stable_value(confidence), relative_width > 0.6 or confidence < 0.55


def _continuous_interpretation(
    uncertainty_width: float, trait_range: float, *, unstable: bool
) -> str:
    if uncertainty_width == 0.0:
        return "observed tip value"
    scale = max(trait_range, 1e-12)
    relative_width = uncertainty_width / scale
    if unstable:
        return "unstable node estimate"
    if relative_width <= 0.25:
        return "narrow uncertainty"
    if relative_width <= 0.6:
        return "moderate uncertainty"
    return "broad uncertainty"


def _continuous_downstream_risks(unstable: bool) -> list[str]:
    if not unstable:
        return []
    return [
        "node ordering and trait-polarity interpretations may change across alternative trees or models",
        "publication claims about deep ancestral values should be treated as provisional",
    ]


def _summarize_brownian_fit_diagnostics(
    dataset: AncestralContinuousDataset,
) -> ContinuousAncestralBrownianFitDiagnostics:
    trait_values = [dataset.values_by_taxon[taxon] for taxon in dataset.taxa]
    covariance_matrix = _build_brownian_covariance_matrix(dataset.tree, dataset.taxa)
    covariance_matrix_rank = _matrix_rank(covariance_matrix, tolerance=1e-12)
    covariance_singular = covariance_matrix_rank < len(dataset.taxa)
    covariance_condition_number = math.inf
    if not covariance_singular:
        covariance_condition_number = _symmetric_matrix_condition_number(
            covariance_matrix
        )
    covariance_near_singular = (
        covariance_singular
        or covariance_condition_number >= _BROWNIAN_CONDITION_THRESHOLD
    )
    covariance_positive_definite, covariance_log_determinant = (
        _matrix_positive_definite_diagnostics(covariance_matrix)
    )
    solver_regularized = False
    solver_regularization_epsilon: float | None = None
    solver_covariance = covariance_matrix
    try:
        inverse_covariance = _invert_matrix(solver_covariance)
        solver_log_determinant = _log_determinant(solver_covariance)
    except ValueError:
        solver_regularized = True
        solver_regularization_epsilon = _SOLVER_REGULARIZATION_EPSILON
        solver_covariance = _stable_covariance(
            covariance_matrix,
            epsilon=_SOLVER_REGULARIZATION_EPSILON,
        )
        inverse_covariance = _invert_matrix(solver_covariance)
        solver_log_determinant = _log_determinant(solver_covariance)
    ones = [1.0] * len(trait_values)
    denominator = _quadratic_form(ones, inverse_covariance)
    root_estimate = (
        sum(
            ones[row_index]
            * sum(
                inverse_covariance[row_index][column_index]
                * trait_values[column_index]
                for column_index in range(len(trait_values))
            )
            for row_index in range(len(trait_values))
        )
        / denominator
    )
    residuals = [value - root_estimate for value in trait_values]
    residual_sigma_squared = max(
        _quadratic_form(residuals, inverse_covariance) / len(trait_values),
        1e-12,
    )
    log_likelihood = -0.5 * (
        len(trait_values) * math.log(2.0 * math.pi * residual_sigma_squared)
        + solver_log_determinant
        + len(trait_values)
    )
    root_depths = _tip_root_depths(dataset.tree, dataset.taxa)
    ultrametric_summary = summarize_ultrametric_tip_depths(root_depths, tolerance=1e-12)
    minimum_branch_length, maximum_branch_length = _branch_length_range(dataset.tree)
    return ContinuousAncestralBrownianFitDiagnostics(
        covariance_model="brownian-shared-path",
        tree_is_ultrametric=ultrametric_summary.ultrametric,
        minimum_root_to_tip_depth=stable_value(ultrametric_summary.minimum_tip_depth),
        maximum_root_to_tip_depth=stable_value(ultrametric_summary.maximum_tip_depth),
        minimum_branch_length=stable_value(minimum_branch_length),
        maximum_branch_length=stable_value(maximum_branch_length),
        covariance_matrix_dimension=len(dataset.taxa),
        covariance_matrix_rank=covariance_matrix_rank,
        covariance_singular=covariance_singular,
        covariance_near_singular=covariance_near_singular,
        covariance_positive_definite=covariance_positive_definite,
        covariance_condition_number=stable_value(covariance_condition_number),
        covariance_log_determinant=(
            None
            if covariance_log_determinant is None
            else stable_value(covariance_log_determinant)
        ),
        solver_name="gauss-jordan-inverse",
        solver_regularized=solver_regularized,
        solver_regularization_epsilon=solver_regularization_epsilon,
        log_likelihood=stable_value(log_likelihood),
        residual_sigma_squared=stable_value(residual_sigma_squared),
    )


def _matrix_positive_definite_diagnostics(
    covariance_matrix: list[list[float]],
) -> tuple[bool, float | None]:
    try:
        _invert_matrix(covariance_matrix)
        return True, _log_determinant(covariance_matrix)
    except ValueError:
        return False, None


def _branch_length_range(tree: PhyloTree) -> tuple[float, float]:
    branch_lengths = [
        float(node.branch_length)
        for node in tree.iter_nodes()
        if node is not tree.root and node.branch_length is not None
    ]
    if not branch_lengths:
        return 0.0, 0.0
    return min(branch_lengths), max(branch_lengths)


def _build_brownian_covariance_matrix(
    tree: PhyloTree, taxa: list[str]
) -> list[list[float]]:
    leaf_paths = _leaf_ancestor_depths(tree)
    matrix: list[list[float]] = []
    for left_taxon in taxa:
        left_path = leaf_paths[left_taxon]
        row: list[float] = []
        for right_taxon in taxa:
            right_path = leaf_paths[right_taxon]
            shared_ancestor_ids = set(left_path) & set(right_path)
            shared_depth = max(left_path[node_id] for node_id in shared_ancestor_ids)
            row.append(shared_depth)
        matrix.append(row)
    return matrix


def _tip_root_depths(tree: PhyloTree, taxa: list[str]) -> dict[str, float]:
    leaf_paths = _leaf_ancestor_depths(tree)
    return {taxon: max(leaf_paths[taxon].values()) for taxon in taxa}


def _leaf_ancestor_depths(tree: PhyloTree) -> dict[str, dict[str, float]]:
    depths_by_leaf: dict[str, dict[str, float]] = {}

    def visit(node: TreeNode, depth: float, path: dict[str, float]) -> None:
        branch_length = 0.0 if node is tree.root else float(node.branch_length or 0.0)
        current_depth = depth + branch_length
        current_path = dict(path)
        current_path[node.node_id or ""] = current_depth
        if node.is_leaf():
            if node.name is None:
                raise ValueError("leaf taxon name is required for ancestral reconstruction")
            depths_by_leaf[node.name] = current_path
            return
        for child in node.children:
            visit(child, current_depth, current_path)

    visit(tree.root, 0.0, {})
    return depths_by_leaf


def _quadratic_form(vector: list[float], matrix: list[list[float]]) -> float:
    return sum(
        vector[row_index]
        * sum(
            matrix[row_index][column_index] * vector[column_index]
            for column_index in range(len(vector))
        )
        for row_index in range(len(vector))
    )


def _invert_matrix(matrix: list[list[float]]) -> list[list[float]]:
    size = len(matrix)
    if size == 0:
        return []
    augmented = [
        [float(value) for value in row]
        + [1.0 if row_index == column_index else 0.0 for column_index in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for pivot_index in range(size):
        pivot_row = max(
            range(pivot_index, size),
            key=lambda row_index: abs(augmented[row_index][pivot_index]),
        )
        pivot_value = augmented[pivot_row][pivot_index]
        if math.isclose(pivot_value, 0.0, abs_tol=1e-12):
            raise ValueError("matrix is singular and cannot be inverted")
        if pivot_row != pivot_index:
            augmented[pivot_index], augmented[pivot_row] = (
                augmented[pivot_row],
                augmented[pivot_index],
            )
        pivot_value = augmented[pivot_index][pivot_index]
        augmented[pivot_index] = [
            value / pivot_value for value in augmented[pivot_index]
        ]
        for row_index in range(size):
            if row_index == pivot_index:
                continue
            factor = augmented[row_index][pivot_index]
            if math.isclose(factor, 0.0, abs_tol=1e-15):
                continue
            augmented[row_index] = [
                row_value - factor * pivot_value
                for row_value, pivot_value in zip(
                    augmented[row_index], augmented[pivot_index], strict=True
                )
            ]
    return [row[size:] for row in augmented]


def _log_determinant(matrix: list[list[float]]) -> float:
    size = len(matrix)
    if size == 0:
        return 0.0
    working = [list(row) for row in matrix]
    sign = 1.0
    log_abs_det = 0.0
    for pivot_index in range(size):
        pivot_row = max(
            range(pivot_index, size),
            key=lambda row_index: abs(working[row_index][pivot_index]),
        )
        pivot_value = working[pivot_row][pivot_index]
        if math.isclose(pivot_value, 0.0, abs_tol=1e-12):
            raise ValueError("matrix determinant is zero")
        if pivot_row != pivot_index:
            working[pivot_index], working[pivot_row] = (
                working[pivot_row],
                working[pivot_index],
            )
            sign *= -1.0
        pivot_value = working[pivot_index][pivot_index]
        if pivot_value < 0.0:
            sign *= -1.0
        log_abs_det += math.log(abs(pivot_value))
        for row_index in range(pivot_index + 1, size):
            factor = working[row_index][pivot_index] / pivot_value
            if math.isclose(factor, 0.0, abs_tol=1e-15):
                continue
            for column_index in range(pivot_index, size):
                working[row_index][column_index] -= (
                    factor * working[pivot_index][column_index]
                )
    if sign <= 0.0:
        raise ValueError("matrix determinant is not positive")
    return log_abs_det


def _stable_covariance(
    matrix: list[list[float]], *, epsilon: float
) -> list[list[float]]:
    stabilized = [list(row) for row in matrix]
    for index in range(len(stabilized)):
        stabilized[index][index] += epsilon
    return stabilized


def _symmetric_matrix_eigenvalues(
    matrix: list[list[float]],
    *,
    tolerance: float = 1e-15,
    max_iterations: int = 10_000,
) -> list[float]:
    size = len(matrix)
    if size == 0:
        return []
    if size == 1:
        return [float(matrix[0][0])]
    working = [list(row) for row in matrix]
    for _ in range(max_iterations):
        pivot_row = 0
        pivot_column = 1
        pivot_value = 0.0
        for row_index in range(size):
            for column_index in range(row_index + 1, size):
                candidate = abs(working[row_index][column_index])
                if candidate > pivot_value:
                    pivot_row = row_index
                    pivot_column = column_index
                    pivot_value = candidate
        if pivot_value <= tolerance:
            return [working[index][index] for index in range(size)]
        app = working[pivot_row][pivot_row]
        aqq = working[pivot_column][pivot_column]
        apq = working[pivot_row][pivot_column]
        tau = (aqq - app) / (2.0 * apq)
        tangent = (
            math.copysign(1.0, tau) / (abs(tau) + math.sqrt(1.0 + tau * tau))
            if not math.isclose(tau, 0.0, abs_tol=tolerance)
            else 1.0
        )
        cosine = 1.0 / math.sqrt(1.0 + tangent * tangent)
        sine = tangent * cosine
        for index in range(size):
            if index in (pivot_row, pivot_column):
                continue
            left = working[index][pivot_row]
            right = working[index][pivot_column]
            working[index][pivot_row] = working[pivot_row][index] = (
                cosine * left - sine * right
            )
            working[index][pivot_column] = working[pivot_column][index] = (
                sine * left + cosine * right
            )
        working[pivot_row][pivot_row] = (
            cosine * cosine * app
            - 2.0 * sine * cosine * apq
            + sine * sine * aqq
        )
        working[pivot_column][pivot_column] = (
            sine * sine * app
            + 2.0 * sine * cosine * apq
            + cosine * cosine * aqq
        )
        working[pivot_row][pivot_column] = 0.0
        working[pivot_column][pivot_row] = 0.0
    raise ValueError("symmetric eigenvalue iteration did not converge")


def _symmetric_matrix_condition_number(
    matrix: list[list[float]], *, tolerance: float = 1e-12
) -> float:
    singular_values = sorted(
        abs(value)
        for value in _symmetric_matrix_eigenvalues(matrix, tolerance=tolerance)
    )
    if not singular_values:
        return 0.0
    if math.isclose(singular_values[0], 0.0, abs_tol=tolerance):
        return math.inf
    return singular_values[-1] / singular_values[0]


def _matrix_rank(matrix: list[list[float]], *, tolerance: float) -> int:
    if not matrix:
        return 0
    working = [list(map(float, row)) for row in matrix]
    row_count = len(working)
    column_count = len(working[0])
    rank = 0
    pivot_row = 0
    for pivot_column in range(column_count):
        candidate_row = max(
            range(pivot_row, row_count),
            key=lambda index: abs(working[index][pivot_column]),
        )
        pivot_value = working[candidate_row][pivot_column]
        if math.isclose(pivot_value, 0.0, abs_tol=tolerance):
            continue
        working[pivot_row], working[candidate_row] = (
            working[candidate_row],
            working[pivot_row],
        )
        pivot = working[pivot_row][pivot_column]
        for column_index in range(pivot_column, column_count):
            working[pivot_row][column_index] /= pivot
        for row_index in range(row_count):
            if row_index == pivot_row:
                continue
            factor = working[row_index][pivot_column]
            if math.isclose(factor, 0.0, abs_tol=tolerance):
                continue
            for column_index in range(pivot_column, column_count):
                working[row_index][column_index] -= (
                    factor * working[pivot_row][column_index]
                )
        rank += 1
        pivot_row += 1
        if pivot_row == row_count:
            break
    return rank
