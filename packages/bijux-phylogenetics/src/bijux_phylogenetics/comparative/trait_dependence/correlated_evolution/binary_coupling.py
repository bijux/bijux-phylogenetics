from __future__ import annotations

from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative.discrete_evolution import (
    run_discrete_state_transition_model,
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
    _fisher_interval,
    _sample_covariance_and_correlation,
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
