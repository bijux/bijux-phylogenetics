from __future__ import annotations

from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative._math import stable_covariance
from bijux_phylogenetics.comparative.signal import (
    compute_phylogenetic_independent_contrasts,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa

from .contracts import CorrelatedTraitEvolutionReport, CorrelatedTraitObservationRow
from .preparation import _PreparedTraitRows
from .statistics import (
    _aic,
    _chi_square_survival,
    _comparison_rows,
    _correlation,
    _estimate_trait_covariance,
    _fisher_interval,
    _multivariate_normal_log_likelihood,
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
