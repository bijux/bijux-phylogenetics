# ruff: noqa: F401, F403, F405
from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from hashlib import sha256
from html import escape
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from Bio import Phylo

from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.tree_set import (
    DiscreteAncestralTreeSetReport,
    summarize_discrete_ancestral_tree_set,
)
from bijux_phylogenetics.biogeography import (
    BiogeographyReportPackageResult,
    build_biogeography_report_package,
)
from bijux_phylogenetics.trees import (
    CladeTableReport,
    CladeTableRow,
    extract_tree_clades,
    write_clade_table,
)
from bijux_phylogenetics.comparative.pgls import (
    PGLSResult,
    run_pgls,
    write_pgls_model_matrix_table,
)
from bijux_phylogenetics.comparative.pgls_categorical_contrasts import (
    PGLSCategoricalContrastReport,
    summarize_pgls_categorical_contrasts,
    write_pgls_categorical_contrast_table,
)
from bijux_phylogenetics.comparative.pgls_lambda_fit import (
    write_pgls_lambda_profile_table,
)
from bijux_phylogenetics.comparative.posterior_tree_pgls import (
    PosteriorTreePGLSReport,
    run_posterior_tree_pgls,
)
from bijux_phylogenetics.comparative.report_package import (
    ComparativeAnalysisSummaryRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeResidualTableRow,
    ComparativeSignalTableRow,
    summarize_comparative_analysis,
    summarize_comparative_audit,
    summarize_comparative_coefficients,
    summarize_comparative_interpretation,
    summarize_comparative_residuals,
    summarize_comparative_signal,
    write_comparative_audit_table,
    write_comparative_coefficient_table,
    write_comparative_contrast_table,
    write_comparative_interpretation_table,
    write_comparative_model_comparison_table,
    write_comparative_residual_table,
    write_comparative_signal_table,
    write_comparative_summary_table,
)
from bijux_phylogenetics.comparative.reporting import (
    ComparativeMethodReport,
    build_comparative_method_report,
)
from bijux_phylogenetics.compare.reports import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.compare.topology import write_tree_comparison_table
from bijux_phylogenetics.core.alignment import (
    AlignmentQualityReport,
    SequenceQualityRankingReport,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.topology import (
    TreeRootingReport,
    root_tree_on_outgroup,
    write_tree_rooting_report,
)
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    RabiesMethodSensitivityPanelWorkflowReport,
    run_rabies_method_sensitivity_panel_workflow,
)
from bijux_phylogenetics.diagnostics.conclusion_stability import (
    ConclusionStabilityReport,
    build_ancestral_state_stability_rows,
    build_comparative_coefficient_stability_rows,
    build_conclusion_stability_report,
    build_key_clade_stability_rows,
    build_support_value_stability_rows,
    write_ancestral_state_stability_table,
    write_comparative_coefficient_stability_table,
    write_conclusion_stability_report_html,
    write_conclusion_stability_summary_table,
    write_key_clade_stability_table,
    write_support_value_stability_table,
)
from bijux_phylogenetics.engines.fasta_to_tree import (
    FastaToTreeWorkflowReport,
    run_fasta_to_tree_workflow,
)
from bijux_phylogenetics.host_association import (
    HostSwitchingReport,
    summarize_host_switching,
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_count_table,
    write_host_switch_exclusion_table,
    write_host_switch_fit_table,
    write_host_switch_summary_table,
    write_unsupported_host_switch_claim_table,
)
from bijux_phylogenetics.io.biopython import tree_from_biophylo, tree_to_biophylo
from bijux_phylogenetics.io.fasta._shared import load_permissive_fasta_records
from bijux_phylogenetics.io.fasta.quality import (
    build_alignment_quality_report,
    build_sequence_quality_ranking,
)
from bijux_phylogenetics.io.fasta.records import validate_fasta_input
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.simulation import write_tree_set
from bijux_phylogenetics.trees import (
    BootstrapTreeSetArtifactReport,
    BootstrapTreeSetSummaryReport,
    compute_clade_frequency_table,
    write_bootstrap_tree_set_artifacts,
)

from .config import (
    _build_workflow_config_audit_rows,
    _load_workflow_config,
    _raise_for_failed_config_audit,
    load_rabies_cross_host_geography_panel_dataset,
)
from .models import *

def run_rabies_cross_host_geography_panel_workflow(
    out_dir: Path,
    *,
    config_path: Path | None = None,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    fasttree_executable: str | Path = "FastTree",
    iqtree_seed: int | None = None,
    iqtree_threads: int | None = None,
    bootstrap_replicates: int | None = None,
) -> RabiesCrossHostGeographyPanelWorkflowReport:
    """Run the full integrated rabies workflow from sequences and metadata."""
    dataset = load_rabies_cross_host_geography_panel_dataset(config_path)
    config = _load_workflow_config(config_path)
    config_audit_rows = _build_workflow_config_audit_rows(config)
    _raise_for_failed_config_audit(config_audit_rows)
    workflow = run_fasta_to_tree_workflow(
        dataset.sequences_path,
        out_dir=out_dir,
        prefix=dataset.workflow_prefix,
        sequence_type=dataset.sequence_type,
        mafft_executable=mafft_executable,
        alignment_mode=config.alignment_mode,
        trimal_executable=trimal_executable,
        trimming_mode=config.trimming_mode,
        iqtree_executable=iqtree_executable,
        iqtree_seed=config.iqtree_seed if iqtree_seed is None else iqtree_seed,
        iqtree_threads=(
            config.iqtree_threads if iqtree_threads is None else iqtree_threads
        ),
        trim_gap_threshold=config.trim_gap_threshold,
        bootstrap_replicates=(
            config.bootstrap_replicates
            if bootstrap_replicates is None
            else bootstrap_replicates
        ),
        timeout_seconds=config.timeout_seconds,
    )
    rooted_tree, rooting_report = root_tree_on_outgroup(
        workflow.output_paths["tree"],
        outgroup_taxa=list(dataset.outgroup_taxa),
    )
    rooted_tree_path = out_dir / f"{dataset.workflow_prefix}.rooted.tree"
    write_newick(rooted_tree_path, rooted_tree)
    host_switching = summarize_host_switching(
        rooted_tree_path,
        dataset.metadata_path,
        trait=dataset.host_trait,
        taxon_column="taxon",
        model=dataset.host_model,
    )
    biogeography_report = build_biogeography_report_package(
        tree_path=rooted_tree_path,
        traits_path=dataset.metadata_path,
        centroids_path=dataset.centroids_path,
        trait=dataset.geography_trait,
        out_dir=out_dir / "biogeography-report",
        taxon_column="taxon",
        model=dataset.geography_model,
        region_column="region",
        latitude_column="latitude",
        longitude_column="longitude",
    )
    aligned_quality = build_alignment_quality_report(workflow.output_paths["alignment"])
    trimmed_quality = build_alignment_quality_report(
        workflow.output_paths["trimmed_alignment"]
    )
    trimmed_sequence_ranking = build_sequence_quality_ranking(
        workflow.output_paths["trimmed_alignment"]
    )
    comparative_traits_rows = _build_comparative_trait_rows(
        metadata_path=dataset.metadata_path,
        centroids_path=dataset.centroids_path,
        host_trait=dataset.host_trait,
        geography_trait=dataset.geography_trait,
    )
    comparative_traits_path = (
        out_dir / f"{dataset.workflow_prefix}.comparative-traits.tsv"
    )
    write_taxon_rows(
        comparative_traits_path,
        columns=list(comparative_traits_rows[0].keys()),
        rows=comparative_traits_rows,
    )
    comparative_tree_path, comparative_branch_repairs = _write_comparative_tree(
        rooted_tree_path,
        out_path=out_dir / f"{dataset.workflow_prefix}.comparative.tree",
        branch_length_floor=config.comparative_branch_length_floor,
    )
    comparative_report = build_comparative_method_report(
        comparative_tree_path,
        comparative_traits_path,
        formula=config.comparative_formula,
        taxon_column="taxon",
        lambda_value="estimate",
    )
    comparative_categorical_contrasts = summarize_pgls_categorical_contrasts(
        comparative_tree_path,
        comparative_traits_path,
        formula=config.comparative_formula,
        taxon_column="taxon",
        lambda_value="estimate",
    )
    method_sensitivity_report = run_rabies_method_sensitivity_panel_workflow(
        out_dir / "method-sensitivity-review",
        mafft_executable=mafft_executable,
        trimal_executable=trimal_executable,
        iqtree_executable=iqtree_executable,
        fasttree_executable=fasttree_executable,
        iqtree_seed=config.iqtree_seed if iqtree_seed is None else iqtree_seed,
        iqtree_threads=(
            config.iqtree_threads if iqtree_threads is None else iqtree_threads
        ),
        bootstrap_replicates=(
            config.bootstrap_replicates
            if bootstrap_replicates is None
            else bootstrap_replicates
        ),
    )
    rooted_bootstrap_tree_set_path = _write_rooted_tree_set_on_outgroup(
        workflow.bootstrap_workflow.output_paths["bootstrap_trees"],
        out_path=out_dir / f"{dataset.workflow_prefix}.rooted-bootstrap.trees",
        outgroup_taxa=list(dataset.outgroup_taxa),
    )
    comparative_bootstrap_tree_set_path = _write_comparative_tree_set(
        rooted_bootstrap_tree_set_path,
        out_path=out_dir / f"{dataset.workflow_prefix}.comparative-bootstrap.trees",
        reference_tree_path=comparative_tree_path,
        branch_length_floor=config.comparative_branch_length_floor,
    )
    host_ancestral_report = reconstruct_discrete_ancestral_states(
        rooted_tree_path,
        dataset.metadata_path,
        trait=dataset.host_trait,
        taxon_column="taxon",
        model=dataset.host_model,
    )
    host_ancestral_tree_set_report = summarize_discrete_ancestral_tree_set(
        rooted_bootstrap_tree_set_path,
        dataset.metadata_path,
        trait=dataset.host_trait,
        taxon_column="taxon",
        model=_canonicalize_discrete_tree_set_model(dataset.host_model),
    )
    geography_ancestral_report = reconstruct_discrete_ancestral_states(
        rooted_tree_path,
        dataset.metadata_path,
        trait=dataset.geography_trait,
        taxon_column="taxon",
        model=dataset.geography_model,
    )
    geography_ancestral_tree_set_report = summarize_discrete_ancestral_tree_set(
        rooted_bootstrap_tree_set_path,
        dataset.metadata_path,
        trait=dataset.geography_trait,
        taxon_column="taxon",
        model=_canonicalize_discrete_tree_set_model(dataset.geography_model),
    )
    comparative_posterior_tree_report = run_posterior_tree_pgls(
        comparative_bootstrap_tree_set_path,
        comparative_traits_path,
        formula=config.comparative_formula,
        taxon_column="taxon",
        lambda_value="estimate",
    )
    conclusion_stability_report = _build_conclusion_stability_report(
        rooted_tree_path=rooted_tree_path,
        rooted_bootstrap_tree_set_path=rooted_bootstrap_tree_set_path,
        comparative_traits_path=comparative_traits_path,
        comparative_tree_path=comparative_tree_path,
        comparative_report=comparative_report,
        comparative_posterior_tree_report=comparative_posterior_tree_report,
        method_sensitivity_report=method_sensitivity_report,
        metadata_path=dataset.metadata_path,
        host_trait=dataset.host_trait,
        host_model=dataset.host_model,
        geography_trait=dataset.geography_trait,
        geography_model=dataset.geography_model,
        comparative_formula=config.comparative_formula,
        comparative_branch_length_floor=config.comparative_branch_length_floor,
        host_ancestral_report=host_ancestral_report,
        host_ancestral_tree_set_report=host_ancestral_tree_set_report,
        geography_ancestral_report=geography_ancestral_report,
        geography_ancestral_tree_set_report=geography_ancestral_tree_set_report,
    )
    return RabiesCrossHostGeographyPanelWorkflowReport(
        dataset=dataset,
        config=config,
        config_audit_rows=config_audit_rows,
        fasta_to_tree=workflow,
        rooted_tree_path=rooted_tree_path,
        rooting_report=rooting_report,
        host_switching=host_switching,
        biogeography_report=biogeography_report,
        aligned_quality=aligned_quality,
        trimmed_quality=trimmed_quality,
        trimmed_sequence_ranking=trimmed_sequence_ranking,
        comparative_traits_rows=comparative_traits_rows,
        comparative_tree_path=comparative_tree_path,
        comparative_branch_repairs=comparative_branch_repairs,
        comparative_report=comparative_report,
        comparative_categorical_contrasts=comparative_categorical_contrasts,
        method_sensitivity_report=method_sensitivity_report,
        rooted_bootstrap_tree_set_path=rooted_bootstrap_tree_set_path,
        comparative_bootstrap_tree_set_path=comparative_bootstrap_tree_set_path,
        host_ancestral_report=host_ancestral_report,
        host_ancestral_tree_set_report=host_ancestral_tree_set_report,
        geography_ancestral_report=geography_ancestral_report,
        geography_ancestral_tree_set_report=geography_ancestral_tree_set_report,
        comparative_posterior_tree_report=comparative_posterior_tree_report,
        conclusion_stability_report=conclusion_stability_report,
    )



def _build_comparative_trait_rows(
    *,
    metadata_path: Path,
    centroids_path: Path,
    host_trait: str,
    geography_trait: str,
) -> list[dict[str, str]]:
    centroids_by_region: dict[str, dict[str, str]] = {}
    with centroids_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            centroids_by_region[row["region"].strip()] = row
    rows: list[dict[str, str]] = []
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            region = row[geography_trait].strip()
            centroid = centroids_by_region[region]
            rows.append(
                {
                    "taxon": row["taxon"].strip(),
                    "host_group": row[host_trait].strip(),
                    "region_group": region,
                    "region_latitude": centroid["latitude"].strip(),
                    "region_longitude": centroid["longitude"].strip(),
                }
            )
    return rows


def _write_comparative_tree(
    rooted_tree_path: Path,
    *,
    out_path: Path,
    branch_length_floor: float,
) -> tuple[Path, list[RabiesComparativeBranchRepair]]:
    tree = load_tree(rooted_tree_path)
    repairs = _apply_branch_length_floor(tree.root, floor=branch_length_floor)
    write_newick(out_path, tree)
    return out_path, repairs


def _write_comparative_tree_set(
    rooted_tree_set_path: Path,
    *,
    out_path: Path,
    reference_tree_path: Path,
    branch_length_floor: float,
) -> Path:
    reference_length_lookup = _build_branch_length_lookup(
        load_tree(reference_tree_path)
    )
    rooted_trees = _load_tree_set_trees(rooted_tree_set_path)
    adjusted_trees = []
    for tree in rooted_trees:
        _overlay_branch_lengths_from_reference(
            tree.root,
            reference_length_lookup=reference_length_lookup,
            floor=branch_length_floor,
        )
        _apply_branch_length_floor(tree.root, floor=branch_length_floor)
        adjusted_trees.append(tree)
    return write_tree_set(out_path, adjusted_trees)


def _write_rooted_tree_set_on_outgroup(
    tree_set_path: Path,
    *,
    out_path: Path,
    outgroup_taxa: list[str],
) -> Path:
    rooted_trees = [
        _root_tree_on_outgroup_from_tree(tree, outgroup_taxa=outgroup_taxa)
        for tree in _load_tree_set_trees(tree_set_path)
    ]
    return write_tree_set(out_path, rooted_trees)


def _load_tree_set_trees(path: Path) -> list[PhyloTree]:
    source_format = "newick"
    return [
        tree_from_biophylo(tree, source_format=source_format)
        for tree in Phylo.parse(path, source_format)
    ]


def _root_tree_on_outgroup_from_tree(
    tree: PhyloTree, *, outgroup_taxa: list[str]
) -> PhyloTree:
    biophylo_tree = tree_to_biophylo(tree)
    matched = [
        next(biophylo_tree.find_clades(name=taxon), None) for taxon in outgroup_taxa
    ]
    if not any(clade is not None for clade in matched):
        raise ValueError(
            "none of the requested outgroup taxa were found while rooting a tree set"
        )
    biophylo_tree.root_with_outgroup(*[clade for clade in matched if clade is not None])
    return tree_from_biophylo(biophylo_tree, source_format="newick")


def _build_branch_length_lookup(tree: PhyloTree) -> dict[str, float]:
    lookup: dict[str, float] = {}
    for node in tree.iter_nodes():
        if node is tree.root or node.branch_length is None:
            continue
        lookup[_node_signature(node)] = node.branch_length
    return lookup


def _overlay_branch_lengths_from_reference(
    node: TreeNode,
    *,
    reference_length_lookup: dict[str, float],
    floor: float,
) -> None:
    for child in node.children:
        signature = _node_signature(child)
        child.branch_length = reference_length_lookup.get(signature, floor)
        _overlay_branch_lengths_from_reference(
            child,
            reference_length_lookup=reference_length_lookup,
            floor=floor,
        )


def _node_signature(node: TreeNode) -> str:
    if node.is_leaf():
        return node.name or "<unnamed>"
    return "|".join(sorted(_descendant_taxa(node)))


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return taxa


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


def _apply_branch_length_floor(
    root: TreeNode,
    *,
    floor: float,
) -> list[RabiesComparativeBranchRepair]:
    repairs: list[RabiesComparativeBranchRepair] = []

    def visit(node: TreeNode, *, is_root: bool) -> None:
        if not is_root and node.branch_length is not None and node.branch_length <= 0.0:
            repairs.append(
                RabiesComparativeBranchRepair(
                    node_label=node.name or "<internal>",
                    original_branch_length=float(node.branch_length),
                    repaired_branch_length=floor,
                    reason=(
                        "comparative branch-length methods require strictly positive "
                        "nonroot branch lengths"
                    ),
                )
            )
            node.branch_length = floor
        for child in node.children:
            visit(child, is_root=False)

    visit(root, is_root=True)
    return repairs


def _stabilize_clade_report(
    report: CladeTableReport,
    *,
    stable_source_path: Path,
) -> CladeTableReport:
    return CladeTableReport(
        path=stable_source_path,
        source_format=report.source_format,
        tree_count=report.tree_count,
        metadata_path=report.metadata_path,
        taxon_column=report.taxon_column,
        metadata_columns=list(report.metadata_columns),
        rows=[
            CladeTableRow(
                source_path=stable_source_path,
                tree_index=row.tree_index,
                node_kind=row.node_kind,
                clade_id=row.clade_id,
                node_label=row.node_label,
                taxon_count=row.taxon_count,
                taxa=list(row.taxa),
                support=row.support,
                support_fraction=row.support_fraction,
                branch_length=row.branch_length,
                root_depth=row.root_depth,
                descendant_tip_depth_min=row.descendant_tip_depth_min,
                descendant_tip_depth_max=row.descendant_tip_depth_max,
                node_age=row.node_age,
                metadata=list(row.metadata),
            )
            for row in report.rows
        ],
    )

