from __future__ import annotations

import math
from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
)
from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    log_determinant,
    quadratic_form,
    stable_covariance,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    run_discrete_state_transition_model,
)
from bijux_phylogenetics.comparative.signal import (
    compute_phylogenetic_independent_contrasts,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .contracts import CorrelatedTraitComparisonRow
from .contracts import CorrelatedTraitEvolutionReport
from .contracts import CorrelatedTraitExclusion
from .contracts import CorrelatedTraitObservationRow
from .preparation import _PreparedTraitRows
from .preparation import _prepare_shared_trait_rows

_FISHER_95_Z = 1.959963984540054
_LOG_2PI = math.log(2.0 * math.pi)


def summarize_correlated_trait_evolution(
    tree_path: Path,
    traits_path: Path,
    *,
    left_trait: str,
    right_trait: str,
    taxon_column: str | None = None,
    analysis_kind: str = "auto",
    binary_model: str = "all-rates-different",
) -> CorrelatedTraitEvolutionReport:
    """Summarize coupling between two traits across one phylogenetic tree."""
    if left_trait == right_trait:
        raise ComparativeMethodError(
            "correlated trait evolution requires two distinct trait columns"
        )
    if analysis_kind not in {"auto", "continuous", "binary"}:
        raise ComparativeMethodError(
            "analysis_kind must be one of: auto, continuous, binary"
        )
    if binary_model == "meristic":
        _resolve_discrete_model_name(binary_model)
    if binary_model not in {"equal-rates", "symmetric", "all-rates-different"}:
        raise ComparativeMethodError(
            "binary correlated-trait analysis requires a supported discrete model"
        )

    tree = load_tree(tree_path)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    for trait in (left_trait, right_trait):
        if trait not in table.columns:
            raise ComparativeMethodError(
                f"trait table does not contain required column '{trait}'"
            )

    prepared = _prepare_shared_trait_rows(
        tree=tree,
        table=table,
        left_trait=left_trait,
        right_trait=right_trait,
        analysis_kind=analysis_kind,
    )
    if prepared.analysis_kind == "continuous-brownian-contrasts":
        return _summarize_continuous_trait_coupling(
            tree_path=tree_path,
            traits_path=traits_path,
            tree_taxon_count=tree.tip_count,
            taxon_column=table.taxon_column,
            left_trait=left_trait,
            right_trait=right_trait,
            prepared=prepared,
        )
    return _summarize_binary_trait_coupling(
        tree_path=tree_path,
        traits_path=traits_path,
        tree_taxon_count=tree.tip_count,
        taxon_column=table.taxon_column,
        left_trait=left_trait,
        right_trait=right_trait,
        prepared=prepared,
        binary_model=binary_model,
    )


def write_correlated_trait_summary_table(
    path: Path,
    report: CorrelatedTraitEvolutionReport,
) -> Path:
    """Write one summary ledger for correlated trait evolution."""
    return write_taxon_rows(
        path,
        columns=[
            "analysis_kind",
            "left_trait",
            "right_trait",
            "taxon_column",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observation_row_count",
            "association_measure_name",
            "association_measure_value",
            "evolutionary_covariance",
            "evolutionary_correlation",
            "lower_95_confidence_interval",
            "upper_95_confidence_interval",
            "independent_parameter_count",
            "independent_log_likelihood",
            "independent_aic",
            "correlated_parameter_count",
            "correlated_log_likelihood",
            "correlated_aic",
            "better_model",
            "likelihood_ratio_statistic",
            "likelihood_ratio_degrees_of_freedom",
            "likelihood_ratio_p_value",
            "likelihood_ratio_p_value_method",
            "left_root_estimate",
            "right_root_estimate",
            "left_state_order",
            "right_state_order",
            "joint_state_count",
            "warning_count",
        ],
        rows=[
            {
                "analysis_kind": report.analysis_kind,
                "left_trait": report.left_trait,
                "right_trait": report.right_trait,
                "taxon_column": report.taxon_column,
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": len(report.analyzed_taxa),
                "excluded_taxon_count": len(report.excluded_taxa),
                "observation_row_count": len(report.observation_rows),
                "association_measure_name": report.association_measure_name,
                "association_measure_value": format(
                    report.association_measure_value, ".15g"
                ),
                "evolutionary_covariance": _format_optional(
                    report.evolutionary_covariance
                ),
                "evolutionary_correlation": _format_optional(
                    report.evolutionary_correlation
                ),
                "lower_95_confidence_interval": _format_optional(
                    report.lower_95_confidence_interval
                ),
                "upper_95_confidence_interval": _format_optional(
                    report.upper_95_confidence_interval
                ),
                "independent_parameter_count": report.independent_parameter_count,
                "independent_log_likelihood": format(
                    report.independent_log_likelihood, ".15g"
                ),
                "independent_aic": format(report.independent_aic, ".15g"),
                "correlated_parameter_count": report.correlated_parameter_count,
                "correlated_log_likelihood": format(
                    report.correlated_log_likelihood, ".15g"
                ),
                "correlated_aic": format(report.correlated_aic, ".15g"),
                "better_model": report.better_model,
                "likelihood_ratio_statistic": format(
                    report.likelihood_ratio_statistic, ".15g"
                ),
                "likelihood_ratio_degrees_of_freedom": (
                    report.likelihood_ratio_degrees_of_freedom
                ),
                "likelihood_ratio_p_value": format(
                    report.likelihood_ratio_p_value, ".15g"
                ),
                "likelihood_ratio_p_value_method": (
                    report.likelihood_ratio_p_value_method
                ),
                "left_root_estimate": _format_optional(report.left_root_estimate),
                "right_root_estimate": _format_optional(report.right_root_estimate),
                "left_state_order": ",".join(report.left_state_order),
                "right_state_order": ",".join(report.right_state_order),
                "joint_state_count": len(report.joint_state_counts),
                "warning_count": len(report.warnings),
            }
        ],
    )


def write_correlated_trait_comparison_table(
    path: Path,
    report: CorrelatedTraitEvolutionReport,
) -> Path:
    """Write one independent-versus-correlated model comparison ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "model_kind",
            "model_description",
            "parameter_count",
            "log_likelihood",
            "aic",
            "delta_aic",
            "selected",
        ],
        rows=[
            {
                "model_kind": row.model_kind,
                "model_description": row.model_description,
                "parameter_count": row.parameter_count,
                "log_likelihood": format(row.log_likelihood, ".15g"),
                "aic": format(row.aic, ".15g"),
                "delta_aic": format(row.delta_aic, ".15g"),
                "selected": str(row.selected).lower(),
            }
            for row in report.comparison_rows
        ],
    )


def write_correlated_trait_observation_table(
    path: Path,
    report: CorrelatedTraitEvolutionReport,
) -> Path:
    """Write one detailed observation ledger for a coupling analysis."""
    return write_taxon_rows(
        path,
        columns=[
            "row_kind",
            "label",
            "taxon",
            "left_taxa",
            "right_taxa",
            "left_numeric_value",
            "right_numeric_value",
            "expected_variance",
            "left_state",
            "right_state",
            "joint_state",
        ],
        rows=[
            {
                "row_kind": row.row_kind,
                "label": row.label,
                "taxon": row.taxon or "",
                "left_taxa": ",".join(row.left_taxa),
                "right_taxa": ",".join(row.right_taxa),
                "left_numeric_value": _format_optional(row.left_numeric_value),
                "right_numeric_value": _format_optional(row.right_numeric_value),
                "expected_variance": _format_optional(row.expected_variance),
                "left_state": row.left_state or "",
                "right_state": row.right_state or "",
                "joint_state": row.joint_state or "",
            }
            for row in report.observation_rows
        ],
    )


def write_correlated_trait_exclusion_table(
    path: Path,
    report: CorrelatedTraitEvolutionReport,
) -> Path:
    """Write one excluded-taxon ledger for correlated trait evolution."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason", "missing_traits"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
                "missing_traits": ",".join(row.missing_traits),
            }
            for row in report.excluded_taxa
        ],
    )


def _summarize_continuous_trait_coupling(
    *,
    tree_path: Path,
    traits_path: Path,
    tree_taxon_count: int,
    taxon_column: str,
    left_trait: str,
    right_trait: str,
    prepared: _PreparedTraitRows,
) -> CorrelatedTraitEvolutionReport:
    reduced_tree, _ = prune_tree_to_requested_taxa(tree_path, prepared.analyzed_taxa)
    reduced_rows = [
        {
            taxon_column: row[taxon_column],
            left_trait: row[left_trait],
            right_trait: row[right_trait],
        }
        for row in prepared.analyzed_rows
    ]
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-correlated-traits-"
    ) as tmp_dir:
        reduced_tree_path = Path(tmp_dir) / "correlated-traits-tree.nwk"
        reduced_table_path = Path(tmp_dir) / "correlated-traits.tsv"
        reduced_tree_path.write_text(
            dumps_newick(reduced_tree) + "\n",
            encoding="utf-8",
        )
        write_taxon_rows(
            reduced_table_path,
            columns=[taxon_column, left_trait, right_trait],
            rows=reduced_rows,
        )
        left_report = compute_phylogenetic_independent_contrasts(
            reduced_tree_path,
            reduced_table_path,
            trait=left_trait,
            taxon_column=taxon_column,
        )
        right_report = compute_phylogenetic_independent_contrasts(
            reduced_tree_path,
            reduced_table_path,
            trait=right_trait,
            taxon_column=taxon_column,
        )
    left_lookup = {row.node: row for row in left_report.contrasts}
    right_lookup = {row.node: row for row in right_report.contrasts}
    paired_nodes = [node for node in left_lookup if node in right_lookup]
    left_values = [left_lookup[node].contrast for node in paired_nodes]
    right_values = [right_lookup[node].contrast for node in paired_nodes]
    covariance_matrix = _estimate_trait_covariance(left_values, right_values)
    stabilized = False
    try:
        correlated_log_likelihood = _multivariate_normal_log_likelihood(
            observations=[
                [left, right]
                for left, right in zip(left_values, right_values, strict=True)
            ],
            covariance_matrix=covariance_matrix,
        )
    except ValueError:
        covariance_matrix = stable_covariance(covariance_matrix)
        correlated_log_likelihood = _multivariate_normal_log_likelihood(
            observations=[
                [left, right]
                for left, right in zip(left_values, right_values, strict=True)
            ],
            covariance_matrix=covariance_matrix,
        )
        stabilized = True
    independent_covariance = [
        [covariance_matrix[0][0], 0.0],
        [0.0, covariance_matrix[1][1]],
    ]
    independent_log_likelihood = _multivariate_normal_log_likelihood(
        observations=[
            [left, right] for left, right in zip(left_values, right_values, strict=True)
        ],
        covariance_matrix=independent_covariance,
    )
    evolutionary_covariance = covariance_matrix[0][1]
    evolutionary_correlation = _correlation(
        covariance_matrix[0][0],
        covariance_matrix[1][1],
        covariance_matrix[0][1],
    )
    lower_95, upper_95 = _fisher_interval(
        evolutionary_correlation,
        len(left_values),
    )
    comparison_rows = _comparison_rows(
        independent_log_likelihood=independent_log_likelihood,
        independent_parameter_count=2,
        correlated_log_likelihood=correlated_log_likelihood,
        correlated_parameter_count=3,
        independent_description="brownian-diagonal-contrast-covariance",
        correlated_description="brownian-full-contrast-covariance",
    )
    warnings = list(prepared.warnings)
    if stabilized:
        warnings.append(
            "continuous trait coupling required light covariance stabilization because the fitted contrast covariance was nearly singular"
        )
    observation_rows = [
        CorrelatedTraitObservationRow(
            row_kind="contrast",
            label=node,
            taxon=None,
            left_taxa=list(left_lookup[node].left_taxa),
            right_taxa=list(left_lookup[node].right_taxa),
            left_numeric_value=left_lookup[node].contrast,
            right_numeric_value=right_lookup[node].contrast,
            expected_variance=left_lookup[node].expected_variance,
            left_state=None,
            right_state=None,
            joint_state=None,
        )
        for node in paired_nodes
    ]
    better_model = next(row.model_kind for row in comparison_rows if row.selected)
    likelihood_ratio_statistic = max(
        0.0,
        2.0 * (correlated_log_likelihood - independent_log_likelihood),
    )
    return CorrelatedTraitEvolutionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        left_trait=left_trait,
        right_trait=right_trait,
        taxon_column=taxon_column,
        analysis_kind="continuous-brownian-contrasts",
        tree_taxon_count=tree_taxon_count,
        analyzed_taxa=list(prepared.analyzed_taxa),
        excluded_taxa=list(prepared.excluded_taxa),
        observation_rows=observation_rows,
        comparison_rows=comparison_rows,
        association_measure_name="evolutionary_correlation",
        association_measure_value=evolutionary_correlation,
        evolutionary_covariance=evolutionary_covariance,
        evolutionary_correlation=evolutionary_correlation,
        lower_95_confidence_interval=lower_95,
        upper_95_confidence_interval=upper_95,
        independent_parameter_count=2,
        independent_log_likelihood=independent_log_likelihood,
        independent_aic=_aic(2, independent_log_likelihood),
        correlated_parameter_count=3,
        correlated_log_likelihood=correlated_log_likelihood,
        correlated_aic=_aic(3, correlated_log_likelihood),
        better_model=better_model,
        likelihood_ratio_statistic=likelihood_ratio_statistic,
        likelihood_ratio_degrees_of_freedom=1,
        likelihood_ratio_p_value=_chi_square_survival(likelihood_ratio_statistic, 1),
        likelihood_ratio_p_value_method="chi-square-approximation",
        left_root_estimate=left_report.root_estimate,
        right_root_estimate=right_report.root_estimate,
        left_state_order=[],
        right_state_order=[],
        joint_state_counts={},
        warnings=warnings,
    )


def _summarize_binary_trait_coupling(
    *,
    tree_path: Path,
    traits_path: Path,
    tree_taxon_count: int,
    taxon_column: str,
    left_trait: str,
    right_trait: str,
    prepared: _PreparedTraitRows,
    binary_model: str,
) -> CorrelatedTraitEvolutionReport:
    reduced_tree, _ = prune_tree_to_requested_taxa(tree_path, prepared.analyzed_taxa)
    left_codes = {state: index for index, state in enumerate(prepared.left_state_order)}
    right_codes = {
        state: index for index, state in enumerate(prepared.right_state_order)
    }
    reduced_rows = [
        {
            taxon_column: row[taxon_column],
            left_trait: row[left_trait],
            right_trait: row[right_trait],
            "_joint_state": (
                f"{left_codes[row[left_trait]]}{right_codes[row[right_trait]]}"
            ),
        }
        for row in prepared.analyzed_rows
    ]
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-correlated-traits-"
    ) as tmp_dir:
        reduced_tree_path = Path(tmp_dir) / "correlated-traits-tree.nwk"
        reduced_table_path = Path(tmp_dir) / "correlated-traits.tsv"
        reduced_tree_path.write_text(
            dumps_newick(reduced_tree) + "\n",
            encoding="utf-8",
        )
        write_taxon_rows(
            reduced_table_path,
            columns=[taxon_column, left_trait, right_trait, "_joint_state"],
            rows=reduced_rows,
        )
        correlated_report = run_discrete_state_transition_model(
            reduced_tree_path,
            reduced_table_path,
            trait="_joint_state",
            taxon_column=taxon_column,
            model=binary_model,
            allowed_states=["00", "01", "10", "11"],
        )
        left_report = run_discrete_state_transition_model(
            reduced_tree_path,
            reduced_table_path,
            trait=left_trait,
            taxon_column=taxon_column,
            model=binary_model,
            allowed_states=prepared.left_state_order,
        )
        right_report = run_discrete_state_transition_model(
            reduced_tree_path,
            reduced_table_path,
            trait=right_trait,
            taxon_column=taxon_column,
            model=binary_model,
            allowed_states=prepared.right_state_order,
        )
    joint_state_counts: dict[str, int] = {}
    left_numeric: list[float] = []
    right_numeric: list[float] = []
    observation_rows: list[CorrelatedTraitObservationRow] = []
    for row in reduced_rows:
        joint_state = row["_joint_state"]
        joint_state_counts[joint_state] = joint_state_counts.get(joint_state, 0) + 1
        left_value = float(left_codes[row[left_trait]])
        right_value = float(right_codes[row[right_trait]])
        left_numeric.append(left_value)
        right_numeric.append(right_value)
        observation_rows.append(
            CorrelatedTraitObservationRow(
                row_kind="tip-state",
                label=row[taxon_column],
                taxon=row[taxon_column],
                left_taxa=[],
                right_taxa=[],
                left_numeric_value=left_value,
                right_numeric_value=right_value,
                expected_variance=None,
                left_state=row[left_trait],
                right_state=row[right_trait],
                joint_state=joint_state,
            )
        )
    observed_covariance, observed_correlation = _sample_covariance_and_correlation(
        left_numeric,
        right_numeric,
    )
    independent_log_likelihood = (
        left_report.transition_model.pseudo_log_likelihood
        + right_report.transition_model.pseudo_log_likelihood
    )
    independent_parameter_count = (
        left_report.transition_model.parameter_count
        + right_report.transition_model.parameter_count
    )
    correlated_log_likelihood = correlated_report.transition_model.pseudo_log_likelihood
    correlated_parameter_count = correlated_report.transition_model.parameter_count
    comparison_rows = _comparison_rows(
        independent_log_likelihood=independent_log_likelihood,
        independent_parameter_count=independent_parameter_count,
        correlated_log_likelihood=correlated_log_likelihood,
        correlated_parameter_count=correlated_parameter_count,
        independent_description="binary-independent-discrete-pseudo-likelihood",
        correlated_description="binary-correlated-joint-state-pseudo-likelihood",
    )
    likelihood_ratio_statistic = max(
        0.0,
        2.0 * (correlated_log_likelihood - independent_log_likelihood),
    )
    lower_95, upper_95 = _fisher_interval(
        observed_correlation,
        len(left_numeric),
    )
    warnings = list(
        dict.fromkeys(
            [
                *prepared.warnings,
                *correlated_report.warnings,
                *left_report.warnings,
                *right_report.warnings,
            ]
        )
    )
    better_model = next(row.model_kind for row in comparison_rows if row.selected)
    return CorrelatedTraitEvolutionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        left_trait=left_trait,
        right_trait=right_trait,
        taxon_column=taxon_column,
        analysis_kind="binary-joint-state",
        tree_taxon_count=tree_taxon_count,
        analyzed_taxa=list(prepared.analyzed_taxa),
        excluded_taxa=list(prepared.excluded_taxa),
        observation_rows=observation_rows,
        comparison_rows=comparison_rows,
        association_measure_name="phi_correlation",
        association_measure_value=observed_correlation,
        evolutionary_covariance=observed_covariance,
        evolutionary_correlation=observed_correlation,
        lower_95_confidence_interval=lower_95,
        upper_95_confidence_interval=upper_95,
        independent_parameter_count=independent_parameter_count,
        independent_log_likelihood=independent_log_likelihood,
        independent_aic=_aic(independent_parameter_count, independent_log_likelihood),
        correlated_parameter_count=correlated_parameter_count,
        correlated_log_likelihood=correlated_log_likelihood,
        correlated_aic=_aic(correlated_parameter_count, correlated_log_likelihood),
        better_model=better_model,
        likelihood_ratio_statistic=likelihood_ratio_statistic,
        likelihood_ratio_degrees_of_freedom=(
            correlated_parameter_count - independent_parameter_count
        ),
        likelihood_ratio_p_value=_chi_square_survival(
            likelihood_ratio_statistic,
            correlated_parameter_count - independent_parameter_count,
        ),
        likelihood_ratio_p_value_method="chi-square-approximation",
        left_root_estimate=None,
        right_root_estimate=None,
        left_state_order=list(prepared.left_state_order),
        right_state_order=list(prepared.right_state_order),
        joint_state_counts=dict(sorted(joint_state_counts.items())),
        warnings=warnings,
    )


def _estimate_trait_covariance(
    left_values: list[float],
    right_values: list[float],
) -> list[list[float]]:
    count = len(left_values)
    left_variance = sum(value * value for value in left_values) / count
    right_variance = sum(value * value for value in right_values) / count
    covariance = (
        sum(left * right for left, right in zip(left_values, right_values, strict=True))
        / count
    )
    return [
        [left_variance, covariance],
        [covariance, right_variance],
    ]


def _multivariate_normal_log_likelihood(
    *,
    observations: list[list[float]],
    covariance_matrix: list[list[float]],
) -> float:
    inverse_covariance = invert_matrix(covariance_matrix)
    log_det = log_determinant(covariance_matrix)
    dimension = len(covariance_matrix)
    total = 0.0
    for observation in observations:
        total += -0.5 * (
            (dimension * _LOG_2PI)
            + log_det
            + quadratic_form(observation, inverse_covariance)
        )
    return total


def _comparison_rows(
    *,
    independent_log_likelihood: float,
    independent_parameter_count: int,
    correlated_log_likelihood: float,
    correlated_parameter_count: int,
    independent_description: str,
    correlated_description: str,
) -> list[CorrelatedTraitComparisonRow]:
    independent_aic = _aic(independent_parameter_count, independent_log_likelihood)
    correlated_aic = _aic(correlated_parameter_count, correlated_log_likelihood)
    best_aic = min(independent_aic, correlated_aic)
    return [
        CorrelatedTraitComparisonRow(
            model_kind="independent",
            model_description=independent_description,
            parameter_count=independent_parameter_count,
            log_likelihood=independent_log_likelihood,
            aic=independent_aic,
            delta_aic=independent_aic - best_aic,
            selected=math.isclose(independent_aic, best_aic, abs_tol=1e-12),
        ),
        CorrelatedTraitComparisonRow(
            model_kind="correlated",
            model_description=correlated_description,
            parameter_count=correlated_parameter_count,
            log_likelihood=correlated_log_likelihood,
            aic=correlated_aic,
            delta_aic=correlated_aic - best_aic,
            selected=math.isclose(correlated_aic, best_aic, abs_tol=1e-12),
        ),
    ]


def _aic(parameter_count: int, log_likelihood: float) -> float:
    return (2.0 * parameter_count) - (2.0 * log_likelihood)


def _sample_covariance_and_correlation(
    left_values: list[float],
    right_values: list[float],
) -> tuple[float, float]:
    count = len(left_values)
    left_mean = sum(left_values) / count
    right_mean = sum(right_values) / count
    covariance = (
        sum(
            (left - left_mean) * (right - right_mean)
            for left, right in zip(left_values, right_values, strict=True)
        )
        / count
    )
    left_variance = sum((value - left_mean) ** 2 for value in left_values) / count
    right_variance = sum((value - right_mean) ** 2 for value in right_values) / count
    return covariance, _correlation(left_variance, right_variance, covariance)


def _correlation(
    left_variance: float,
    right_variance: float,
    covariance: float,
) -> float:
    denominator = math.sqrt(left_variance * right_variance)
    if math.isclose(denominator, 0.0, abs_tol=1e-12):
        return 0.0
    return covariance / denominator


def _fisher_interval(
    correlation: float,
    count: int,
) -> tuple[float | None, float | None]:
    if count <= 3 or abs(correlation) >= 1.0:
        return None, None
    fisher_z = math.atanh(correlation)
    standard_error = 1.0 / math.sqrt(count - 3)
    return (
        math.tanh(fisher_z - (_FISHER_95_Z * standard_error)),
        math.tanh(fisher_z + (_FISHER_95_Z * standard_error)),
    )


def _chi_square_survival(statistic: float, degrees_of_freedom: int) -> float:
    if statistic <= 0.0 or degrees_of_freedom <= 0:
        return 1.0
    if degrees_of_freedom == 1:
        return math.erfc(math.sqrt(statistic / 2.0))
    z_score = (
        ((statistic / degrees_of_freedom) ** (1.0 / 3.0))
        - (1.0 - (2.0 / (9.0 * degrees_of_freedom)))
    ) / math.sqrt(2.0 / (9.0 * degrees_of_freedom))
    return 0.5 * math.erfc(z_score / math.sqrt(2.0))


def _format_optional(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".15g")
