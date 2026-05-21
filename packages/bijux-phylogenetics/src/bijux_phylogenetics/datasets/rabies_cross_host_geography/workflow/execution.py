# ruff: noqa: F401, F403, F405
from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.ancestral.tree_set import (
    summarize_discrete_ancestral_tree_set,
)
from bijux_phylogenetics.biogeography import build_biogeography_report_package
from bijux_phylogenetics.comparative.pgls.categorical_contrasts import (
    summarize_pgls_categorical_contrasts,
)
from bijux_phylogenetics.comparative.pgls.posterior_tree import (
    run_posterior_tree_pgls,
)
from bijux_phylogenetics.comparative.reporting import build_comparative_method_report
from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    run_rabies_method_sensitivity_panel_workflow,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.ecology import summarize_host_switching
from bijux_phylogenetics.engines.inference import run_fasta_to_tree_workflow
from bijux_phylogenetics.io.fasta.quality import (
    build_alignment_quality_report,
    build_sequence_quality_ranking,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.phylo.topology import root_tree_on_outgroup

from ..config import (
    _build_workflow_config_audit_rows,
    _load_workflow_config,
    _raise_for_failed_config_audit,
    load_rabies_cross_host_geography_panel_dataset,
)
from ..models import *
from .comparative_inputs import _build_comparative_trait_rows
from .conclusion_stability import (
    _build_conclusion_stability_report,
    _canonicalize_discrete_tree_set_model,
)
from .tree_transforms import (
    _apply_branch_length_floor,
    _write_comparative_tree,
    _write_comparative_tree_set,
    _write_rooted_tree_set_on_outgroup,
)


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
