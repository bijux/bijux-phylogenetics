from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.tree_set import DiscreteAncestralTreeSetReport
from bijux_phylogenetics.comparative.pgls import PGLSResult, run_pgls
from bijux_phylogenetics.comparative.pgls.posterior_tree import (
    PosteriorTreePGLSReport,
)
from bijux_phylogenetics.comparative.reporting import ComparativeMethodReport
from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    RabiesMethodSensitivityPanelWorkflowReport,
)
from bijux_phylogenetics.diagnostics.conclusion_stability import (
    ConclusionStabilityReport,
    build_ancestral_state_stability_rows,
    build_comparative_coefficient_stability_rows,
    build_conclusion_stability_report,
    build_key_clade_stability_rows,
    build_support_value_stability_rows,
)
from bijux_phylogenetics.trees import (
    compute_clade_frequency_table,
    extract_tree_clades,
)

from .tree_transforms import _write_comparative_tree


def _build_conclusion_stability_report(
    *,
    rooted_tree_path: Path,
    rooted_bootstrap_tree_set_path: Path,
    comparative_traits_path: Path,
    comparative_tree_path: Path,
    comparative_report: ComparativeMethodReport,
    comparative_posterior_tree_report: PosteriorTreePGLSReport,
    method_sensitivity_report: RabiesMethodSensitivityPanelWorkflowReport,
    metadata_path: Path,
    host_trait: str,
    host_model: str,
    geography_trait: str,
    geography_model: str,
    comparative_formula: str,
    comparative_branch_length_floor: float,
    host_ancestral_report: DiscreteAncestralReport,
    host_ancestral_tree_set_report: DiscreteAncestralTreeSetReport,
    geography_ancestral_report: DiscreteAncestralReport,
    geography_ancestral_tree_set_report: DiscreteAncestralTreeSetReport,
) -> ConclusionStabilityReport:
    baseline_clades = extract_tree_clades(rooted_tree_path)
    bootstrap_frequencies = compute_clade_frequency_table(
        rooted_bootstrap_tree_set_path
    )
    method_tree_paths = [
        variant.rooted_iqtree_path for variant in method_sensitivity_report.variant_runs
    ] + [
        variant.rooted_fasttree_path
        for variant in method_sensitivity_report.variant_runs
    ]
    method_clade_reports = [extract_tree_clades(path) for path in method_tree_paths]
    key_clade_rows = build_key_clade_stability_rows(
        baseline_clades=baseline_clades,
        bootstrap_frequencies=bootstrap_frequencies,
        method_clade_reports=method_clade_reports,
    )
    support_value_rows = build_support_value_stability_rows(
        baseline_clades=baseline_clades,
        bootstrap_frequencies=bootstrap_frequencies,
        method_clade_reports=method_clade_reports,
    )
    host_method_reports = [
        reconstruct_discrete_ancestral_states(
            path,
            metadata_path,
            trait=host_trait,
            taxon_column="taxon",
            model=host_model,
        )
        for path in method_tree_paths
    ]
    geography_method_reports = [
        reconstruct_discrete_ancestral_states(
            path,
            metadata_path,
            trait=geography_trait,
            taxon_column="taxon",
            model=geography_model,
        )
        for path in method_tree_paths
    ]
    ancestral_state_rows = build_ancestral_state_stability_rows(
        baseline_report=host_ancestral_report,
        bootstrap_report=host_ancestral_tree_set_report,
        method_reports=host_method_reports,
    ) + build_ancestral_state_stability_rows(
        baseline_report=geography_ancestral_report,
        bootstrap_report=geography_ancestral_tree_set_report,
        method_reports=geography_method_reports,
    )
    comparative_method_results = [
        _run_comparative_pgls_on_tree(
            tree_path=path,
            comparative_traits_path=comparative_traits_path,
            formula=comparative_formula,
            branch_length_floor=comparative_branch_length_floor,
        )
        for path in method_tree_paths
    ]
    comparative_coefficient_rows = build_comparative_coefficient_stability_rows(
        baseline_result=comparative_report.snapshot.pgls_model,
        bootstrap_report=comparative_posterior_tree_report,
        method_results=comparative_method_results,
    )
    return build_conclusion_stability_report(
        key_clade_rows=key_clade_rows,
        support_value_rows=support_value_rows,
        ancestral_state_rows=ancestral_state_rows,
        comparative_coefficient_rows=comparative_coefficient_rows,
    )


def _canonicalize_discrete_tree_set_model(model: str) -> str:
    normalized = model.strip().lower()
    alias_map = {
        "er": "equal-rates",
        "equal-rates": "equal-rates",
        "sym": "symmetric",
        "symmetric": "symmetric",
        "ard": "all-rates-different",
        "all-rates-different": "all-rates-different",
        "fitch": "fitch",
    }
    return alias_map.get(normalized, model)


def _run_comparative_pgls_on_tree(
    *,
    tree_path: Path,
    comparative_traits_path: Path,
    formula: str,
    branch_length_floor: float,
) -> PGLSResult:
    with TemporaryDirectory(prefix="bijux-rabies-comparative-tree-") as temporary_root:
        adjusted_tree_path, _repairs = _write_comparative_tree(
            tree_path,
            out_path=Path(temporary_root) / "comparative-tree.nwk",
            branch_length_floor=branch_length_floor,
        )
        return run_pgls(
            adjusted_tree_path,
            comparative_traits_path,
            formula=formula,
            taxon_column="taxon",
            lambda_value="estimate",
        )
