from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from hashlib import sha256
from html import escape
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from bijux_phylogenetics.biogeography import (
    BiogeographyReportPackageResult,
    build_biogeography_report_package,
)
from bijux_phylogenetics.core.topology import (
    TreeRootingReport,
    root_tree_on_outgroup,
    write_tree_rooting_report,
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
from bijux_phylogenetics.io.fasta import validate_fasta_input
from bijux_phylogenetics.io.newick import write_newick

_DATASET_ID = "rabies_cross_host_geography_panel"
_DATASET_LABEL = "Rabies cross-host geography panel"
_SEQUENCE_TYPE = "dna"
_WORKFLOW_PREFIX = "rabies-cross-host-geography-panel"
_HOST_TRAIT = "host_group"
_GEOGRAPHY_TRAIT = "region_group"
_HOST_MODEL = "ard"
_GEOGRAPHY_MODEL = "ard"
_IQTREE_SEED = 1
_IQTREE_THREADS = 1
_BOOTSTRAP_REPLICATES = 1000
_OUTGROUP_TAXA = ("bat_chile_rv108",)
_SOURCE_ACCESSIONS = (
    "MG458305",
    "MG458304",
    "PV641713",
    "PX845689",
    "OQ693985",
    "PX845683",
    "PX845681",
    "PX845678",
    "PX845676",
)


@dataclass(slots=True)
class RabiesCrossHostGeographyPanelDataset:
    """Packaged rabies panel for one complete host and geography workflow."""

    dataset_id: str
    label: str
    dataset_root: Path
    sequences_path: Path
    metadata_path: Path
    centroids_path: Path
    reference_output_root: Path
    sequence_count: int
    sequence_type: str
    workflow_prefix: str
    host_trait: str
    geography_trait: str
    host_model: str
    geography_model: str
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    outgroup_taxa: tuple[str, ...]
    observed_host_group_count: int
    observed_region_group_count: int
    source_accessions: tuple[str, ...]
    source_summary: str


@dataclass(slots=True)
class RabiesCrossHostGeographyPanelExportResult:
    """Materialized copy of the packaged rabies integrated dataset."""

    output_root: Path
    readme_path: Path
    sequences_path: Path
    metadata_path: Path
    centroids_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class RabiesCrossHostGeographyPanelWorkflowReport:
    """One full raw-sequence-to-result workflow run over the packaged rabies panel."""

    dataset: RabiesCrossHostGeographyPanelDataset
    fasta_to_tree: FastaToTreeWorkflowReport
    rooted_tree_path: Path
    rooting_report: TreeRootingReport
    host_switching: HostSwitchingReport
    biogeography_report: BiogeographyReportPackageResult


@dataclass(slots=True)
class RabiesCrossHostGeographyPanelWorkflowBundle:
    """Written integrated workflow outputs for the packaged rabies panel."""

    output_root: Path
    selected_model: str
    minimum_support: float | None
    maximum_support: float | None
    median_support: float | None
    weakly_supported_clade_count: int
    rooted_outgroup_taxa: tuple[str, ...]
    root_host: str
    root_host_confidence: float
    host_switch_count: int
    certain_host_switch_count: int
    uncertain_host_switch_count: int
    root_region: str
    root_region_probability: float
    changed_region_branch_count: int
    migration_event_count: int
    strongly_supported_migration_event_count: int
    workflow_summary_path: Path
    alignment_path: Path
    trimmed_alignment_path: Path
    tree_path: Path
    rooting_report_path: Path
    model_table_path: Path
    support_table_path: Path
    log_path: Path
    manifest_path: Path
    engine_artifact_root: Path
    host_switch_summary_path: Path
    host_state_nodes_path: Path
    host_switch_branches_path: Path
    host_switch_counts_path: Path
    host_switch_fits_path: Path
    host_switch_unsupported_path: Path
    host_switch_exclusions_path: Path
    biogeography_output_root: Path
    biogeography_report_path: Path
    biogeography_tree_figure_path: Path
    biogeography_map_path: Path
    final_report_path: Path
    final_manifest_path: Path


@dataclass(slots=True)
class RabiesCrossHostGeographyPanelDemoResult:
    """Dataset export plus integrated workflow outputs for the public rabies demo."""

    output_root: Path
    dataset: RabiesCrossHostGeographyPanelDataset
    dataset_export: RabiesCrossHostGeographyPanelExportResult
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle
    overview_path: Path


def load_rabies_cross_host_geography_panel_dataset() -> (
    RabiesCrossHostGeographyPanelDataset
):
    """Expose the packaged rabies host-and-geography panel as one owned surface."""
    dataset_root = _resource_root()
    sequences_path = dataset_root / "sequences.fasta"
    metadata_path = dataset_root / "metadata.csv"
    centroids_path = dataset_root / "region-centroids.csv"
    validation = validate_fasta_input(sequences_path, sequence_type=_SEQUENCE_TYPE)
    observed_host_groups, observed_region_groups = _read_observed_groups(metadata_path)
    return RabiesCrossHostGeographyPanelDataset(
        dataset_id=_DATASET_ID,
        label=_DATASET_LABEL,
        dataset_root=dataset_root,
        sequences_path=sequences_path,
        metadata_path=metadata_path,
        centroids_path=centroids_path,
        reference_output_root=dataset_root / "expected",
        sequence_count=validation.summary.sequence_count,
        sequence_type=_SEQUENCE_TYPE,
        workflow_prefix=_WORKFLOW_PREFIX,
        host_trait=_HOST_TRAIT,
        geography_trait=_GEOGRAPHY_TRAIT,
        host_model=_HOST_MODEL,
        geography_model=_GEOGRAPHY_MODEL,
        iqtree_seed=_IQTREE_SEED,
        iqtree_threads=_IQTREE_THREADS,
        bootstrap_replicates=_BOOTSTRAP_REPLICATES,
        outgroup_taxa=_OUTGROUP_TAXA,
        observed_host_group_count=len(observed_host_groups),
        observed_region_group_count=len(observed_region_groups),
        source_accessions=_SOURCE_ACCESSIONS,
        source_summary=(
            "Real rabies virus nucleoprotein sequences paired with grouped host "
            "and macroregion metadata so one governed workflow can rerun tree "
            "inference, host switching, and geographic transition review from "
            "raw sequence inputs."
        ),
    )


def export_rabies_cross_host_geography_panel_dataset(
    destination: Path,
) -> RabiesCrossHostGeographyPanelExportResult:
    """Copy the packaged integrated rabies dataset and stable expected outputs."""
    dataset = load_rabies_cross_host_geography_panel_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(
        dataset.dataset_root / "README.md", destination / "README.md"
    )
    sequences_path = shutil.copy2(
        dataset.sequences_path, destination / "sequences.fasta"
    )
    metadata_path = shutil.copy2(dataset.metadata_path, destination / "metadata.csv")
    centroids_path = shutil.copy2(
        dataset.centroids_path, destination / "region-centroids.csv"
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return RabiesCrossHostGeographyPanelExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        sequences_path=Path(sequences_path),
        metadata_path=Path(metadata_path),
        centroids_path=Path(centroids_path),
        expected_output_root=expected_output_root,
    )


def run_rabies_cross_host_geography_panel_workflow(
    out_dir: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = _IQTREE_SEED,
    iqtree_threads: int = _IQTREE_THREADS,
    bootstrap_replicates: int = _BOOTSTRAP_REPLICATES,
) -> RabiesCrossHostGeographyPanelWorkflowReport:
    """Run the full integrated rabies workflow from sequences and metadata."""
    dataset = load_rabies_cross_host_geography_panel_dataset()
    workflow = run_fasta_to_tree_workflow(
        dataset.sequences_path,
        out_dir=out_dir,
        prefix=dataset.workflow_prefix,
        sequence_type=dataset.sequence_type,
        mafft_executable=mafft_executable,
        trimal_executable=trimal_executable,
        iqtree_executable=iqtree_executable,
        iqtree_seed=iqtree_seed,
        iqtree_threads=iqtree_threads,
        bootstrap_replicates=bootstrap_replicates,
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
    return RabiesCrossHostGeographyPanelWorkflowReport(
        dataset=dataset,
        fasta_to_tree=workflow,
        rooted_tree_path=rooted_tree_path,
        rooting_report=rooting_report,
        host_switching=host_switching,
        biogeography_report=biogeography_report,
    )


def write_rabies_cross_host_geography_panel_workflow_bundle(
    output_root: Path,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
) -> RabiesCrossHostGeographyPanelWorkflowBundle:
    """Write the complete integrated workflow bundle for the packaged rabies panel."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    workflow = report.fasta_to_tree
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary

    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report,
    )
    alignment_path = _copy_output(
        workflow.output_paths["alignment"],
        output_root / workflow.output_paths["alignment"].name,
    )
    trimmed_alignment_path = _copy_output(
        workflow.output_paths["trimmed_alignment"],
        output_root / workflow.output_paths["trimmed_alignment"].name,
    )
    tree_path = _copy_output(
        report.rooted_tree_path,
        output_root / report.rooted_tree_path.name,
    )
    stable_rooting_report = replace(
        report.rooting_report, tree_path=Path(tree_path.name)
    )
    rooting_report_path = write_tree_rooting_report(
        output_root / f"{report.dataset.workflow_prefix}.rooting.tsv",
        stable_rooting_report,
    )
    model_table_path = _copy_output(
        workflow.output_paths["model_table"],
        output_root / workflow.output_paths["model_table"].name,
    )
    support_table_path = _copy_output(
        workflow.output_paths["support_table"],
        output_root / workflow.output_paths["support_table"].name,
    )
    log_path = _copy_output(
        workflow.output_paths["log"],
        output_root / workflow.output_paths["log"].name,
    )
    manifest_path = _copy_output(
        workflow.manifest_path,
        output_root / workflow.manifest_path.name,
    )
    engine_artifact_root = (
        output_root / "engine-artifacts" / report.dataset.workflow_prefix
    )
    shutil.copytree(workflow.engine_artifact_dir, engine_artifact_root)

    host_switch_summary_path = write_host_switch_summary_table(
        output_root / "host-switch-summary.tsv",
        report.host_switching,
    )
    host_state_nodes_path = write_host_state_node_table(
        output_root / "host-state-nodes.tsv",
        report.host_switching,
    )
    host_switch_branches_path = write_host_switch_branch_table(
        output_root / "host-switch-branches.tsv",
        report.host_switching,
    )
    host_switch_counts_path = write_host_switch_count_table(
        output_root / "host-switch-counts.tsv",
        report.host_switching,
    )
    host_switch_fits_path = write_host_switch_fit_table(
        output_root / "host-switch-fits.tsv",
        report.host_switching,
    )
    host_switch_unsupported_path = write_unsupported_host_switch_claim_table(
        output_root / "host-switch-unsupported.tsv",
        report.host_switching,
    )
    host_switch_exclusions_path = write_host_switch_exclusion_table(
        output_root / "host-switch-exclusions.tsv",
        report.host_switching,
    )

    biogeography_output_root = output_root / "biogeography"
    shutil.copytree(report.biogeography_report.output_dir, biogeography_output_root)
    biogeography_report_path = biogeography_output_root / "biogeography-report.html"
    biogeography_tree_figure_path = (
        biogeography_output_root / "ancestral-region-tree.svg"
    )
    biogeography_map_path = biogeography_output_root / "geographic-region-map.html"

    final_report_path = _write_integrated_report(
        output_root / "rabies-cross-host-geography-report.html",
        report,
    )
    final_manifest_path = _write_manifest(
        output_root / "rabies-cross-host-geography.manifest.json",
        report,
        bundle_paths={
            "workflow_summary": workflow_summary_path,
            "alignment": alignment_path,
            "trimmed_alignment": trimmed_alignment_path,
            "rooted_tree": tree_path,
            "rooting_report": rooting_report_path,
            "model_table": model_table_path,
            "support_table": support_table_path,
            "host_switch_summary": host_switch_summary_path,
            "host_state_nodes": host_state_nodes_path,
            "host_switch_branches": host_switch_branches_path,
            "host_switch_counts": host_switch_counts_path,
            "host_switch_fits": host_switch_fits_path,
            "host_switch_unsupported": host_switch_unsupported_path,
            "host_switch_exclusions": host_switch_exclusions_path,
            "biogeography_report": biogeography_report_path,
            "biogeography_tree_figure": biogeography_tree_figure_path,
            "biogeography_map": biogeography_map_path,
            "final_report": final_report_path,
        },
    )

    return RabiesCrossHostGeographyPanelWorkflowBundle(
        output_root=output_root,
        selected_model=workflow.selected_model,
        minimum_support=workflow.support_summary.minimum_support,
        maximum_support=workflow.support_summary.maximum_support,
        median_support=workflow.support_summary.median_support,
        weakly_supported_clade_count=workflow.support_summary.weakly_supported_clade_count,
        rooted_outgroup_taxa=tuple(report.rooting_report.rooted_outgroup_taxa),
        root_host=host_summary.root_host,
        root_host_confidence=host_summary.root_confidence,
        host_switch_count=host_summary.host_switch_count,
        certain_host_switch_count=host_summary.certain_host_switch_count,
        uncertain_host_switch_count=host_summary.uncertain_host_switch_count,
        root_region=geography_summary.root_region,
        root_region_probability=geography_summary.root_region_probability,
        changed_region_branch_count=geography_summary.changed_branch_count,
        migration_event_count=migration_summary.event_count,
        strongly_supported_migration_event_count=migration_summary.strongly_supported_event_count,
        workflow_summary_path=workflow_summary_path,
        alignment_path=alignment_path,
        trimmed_alignment_path=trimmed_alignment_path,
        tree_path=tree_path,
        rooting_report_path=rooting_report_path,
        model_table_path=model_table_path,
        support_table_path=support_table_path,
        log_path=log_path,
        manifest_path=manifest_path,
        engine_artifact_root=engine_artifact_root,
        host_switch_summary_path=host_switch_summary_path,
        host_state_nodes_path=host_state_nodes_path,
        host_switch_branches_path=host_switch_branches_path,
        host_switch_counts_path=host_switch_counts_path,
        host_switch_fits_path=host_switch_fits_path,
        host_switch_unsupported_path=host_switch_unsupported_path,
        host_switch_exclusions_path=host_switch_exclusions_path,
        biogeography_output_root=biogeography_output_root,
        biogeography_report_path=biogeography_report_path,
        biogeography_tree_figure_path=biogeography_tree_figure_path,
        biogeography_map_path=biogeography_map_path,
        final_report_path=final_report_path,
        final_manifest_path=final_manifest_path,
    )


def run_rabies_cross_host_geography_panel_demo(
    output_root: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = _IQTREE_SEED,
    iqtree_threads: int = _IQTREE_THREADS,
    bootstrap_replicates: int = _BOOTSTRAP_REPLICATES,
) -> RabiesCrossHostGeographyPanelDemoResult:
    """Materialize the packaged integrated rabies dataset and rerun the full workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_rabies_cross_host_geography_panel_dataset()
    dataset_export = export_rabies_cross_host_geography_panel_dataset(
        output_root / "dataset"
    )
    with TemporaryDirectory(prefix="rabies-cross-host-geography-") as temporary_root:
        workflow_report = run_rabies_cross_host_geography_panel_workflow(
            Path(temporary_root),
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
        )
        workflow_bundle = write_rabies_cross_host_geography_panel_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    overview_path = _write_overview(
        output_root / "overview.md", dataset, workflow_bundle
    )
    return RabiesCrossHostGeographyPanelDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "resources"
        / "datasets"
        / "pathogens"
        / _DATASET_ID
    )


def _read_observed_groups(metadata_path: Path) -> tuple[set[str], set[str]]:
    host_groups: set[str] = set()
    region_groups: set[str] = set()
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            host_group = row.get(_HOST_TRAIT, "").strip()
            region_group = row.get(_GEOGRAPHY_TRAIT, "").strip()
            if host_group:
                host_groups.add(host_group)
            if region_group:
                region_groups.add(region_group)
    return host_groups, region_groups


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _write_workflow_summary_table(
    path: Path,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
) -> Path:
    support = report.fasta_to_tree.support_summary
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary
    lines = [
        (
            "dataset_id\tsequence_count\tselected_model\tminimum_support\tmaximum_support\t"
            "median_support\tweakly_supported_clade_count\toutgroup_taxa\troot_host\t"
            "root_host_confidence\thost_switch_count\tcertain_host_switch_count\t"
            "uncertain_host_switch_count\troot_region\troot_region_probability\t"
            "changed_region_branch_count\tmigration_event_count\t"
            "strongly_supported_migration_event_count"
        ),
        "\t".join(
            [
                report.dataset.dataset_id,
                str(report.dataset.sequence_count),
                report.fasta_to_tree.selected_model,
                _format_number(support.minimum_support),
                _format_number(support.maximum_support),
                _format_number(support.median_support),
                str(support.weakly_supported_clade_count),
                ",".join(report.dataset.outgroup_taxa),
                host_summary.root_host,
                _format_number(host_summary.root_confidence),
                str(host_summary.host_switch_count),
                str(host_summary.certain_host_switch_count),
                str(host_summary.uncertain_host_switch_count),
                geography_summary.root_region,
                _format_number(geography_summary.root_region_probability),
                str(geography_summary.changed_branch_count),
                str(migration_summary.event_count),
                str(migration_summary.strongly_supported_event_count),
            ]
        ),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_overview(
    path: Path,
    dataset: RabiesCrossHostGeographyPanelDataset,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
) -> Path:
    lines = [
        "# Rabies Cross-Host Geography Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- sequence count: `{dataset.sequence_count}`",
        f"- host workflow trait: `{dataset.host_trait}`",
        f"- geography workflow trait: `{dataset.geography_trait}`",
        "",
        "Generated outputs:",
        "",
        f"- integrated workflow summary: `{workflow_bundle.workflow_summary_path.name}`",
        f"- rooted tree: `{workflow_bundle.tree_path.name}`",
        f"- host-switch summary: `{workflow_bundle.host_switch_summary_path.name}`",
        f"- biogeography report: `biogeography/{workflow_bundle.biogeography_report_path.name}`",
        f"- final report: `{workflow_bundle.final_report_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_manifest(
    path: Path,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    *,
    bundle_paths: dict[str, Path],
) -> Path:
    manifest = {
        "report_kind": "rabies_cross_host_geography_workflow_bundle",
        "dataset_id": report.dataset.dataset_id,
        "input_checksums": {
            "sequences.fasta": _checksum(report.dataset.sequences_path),
            "metadata.csv": _checksum(report.dataset.metadata_path),
            "region-centroids.csv": _checksum(report.dataset.centroids_path),
        },
        "output_checksums": {
            key: _checksum(value) for key, value in bundle_paths.items()
        },
        "metrics": {
            "sequence_count": report.dataset.sequence_count,
            "selected_model": report.fasta_to_tree.selected_model,
            "minimum_support": report.fasta_to_tree.support_summary.minimum_support,
            "maximum_support": report.fasta_to_tree.support_summary.maximum_support,
            "host_switch_count": report.host_switching.summary.host_switch_count,
            "migration_event_count": report.biogeography_report.event_report.summary.event_count,
            "root_host": report.host_switching.summary.root_host,
            "root_region": report.biogeography_report.state_report.summary.root_region,
        },
    }
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _write_integrated_report(
    path: Path,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
) -> Path:
    support_summary = report.fasta_to_tree.support_summary
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Bijux Rabies Host and Geography Workflow</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f4f1ea 0%, #e7efe7 100%); color: #163222; }",
            "    main { max-width: 1360px; margin: 0 auto; padding: 24px; }",
            "    h1 { margin: 0 0 8px; font-size: 34px; }",
            "    h2 { margin: 0 0 10px; font-size: 22px; }",
            "    p { line-height: 1.55; }",
            "    .cards { display: grid; grid-template-columns: repeat(5, minmax(0, 1fr)); gap: 14px; margin: 18px 0 24px; }",
            "    .card, .panel { background: rgba(255,255,255,0.86); border: 1px solid rgba(22,50,34,0.12); border-radius: 18px; padding: 18px; box-shadow: 0 16px 42px rgba(22,50,34,0.08); }",
            "    .label { color: #5b7466; font-size: 13px; text-transform: uppercase; letter-spacing: 0.04em; }",
            "    .card strong { display: block; font-size: 21px; margin-top: 6px; }",
            "    .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }",
            "    .full { grid-column: 1 / -1; }",
            "    .figure-shell { overflow: auto; }",
            "    .figure-shell img { width: 100%; height: auto; display: block; }",
            "    table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 14px; }",
            "    th, td { border-bottom: 1px solid rgba(22,50,34,0.10); padding: 8px 10px; text-align: left; vertical-align: top; }",
            "    th { color: #365443; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    a { color: #16543a; }",
            "    iframe { width: 100%; min-height: 760px; border: 1px solid rgba(22,50,34,0.12); border-radius: 14px; background: white; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Bijux Rabies Host and Geography Workflow</h1>",
            "  <p>Complete end-to-end review for one real rabies nucleoprotein panel. The workflow starts from raw sequences plus combined host and geography metadata, infers a maximum-likelihood tree with bootstrap support, roots that tree on one explicit outgroup, then carries host-switching and biogeographic reconstruction into one reviewer-facing handoff.</p>",
            '  <section class="cards">',
            f'    <div class="card"><span class="label">sequences</span><strong>{report.dataset.sequence_count}</strong></div>',
            f'    <div class="card"><span class="label">selected model</span><strong>{escape(report.fasta_to_tree.selected_model)}</strong></div>',
            f'    <div class="card"><span class="label">support range</span><strong>{_support_range_text(support_summary.minimum_support, support_summary.maximum_support)}</strong></div>',
            f'    <div class="card"><span class="label">root host</span><strong>{escape(host_summary.root_host)}</strong></div>',
            f'    <div class="card"><span class="label">root region</span><strong>{escape(geography_summary.root_region)}</strong></div>',
            "  </section>",
            '  <section class="panel">',
            "    <h2>Reviewer Summary</h2>",
            _html_list(
                [
                    f"The inferred rooted tree used substitution model {report.fasta_to_tree.selected_model} with bootstrap support spanning {_support_range_text(support_summary.minimum_support, support_summary.maximum_support)}.",
                    f"The outgroup root is anchored on {', '.join(report.dataset.outgroup_taxa)}.",
                    f"Host reconstruction over {report.dataset.host_trait} inferred {host_summary.host_switch_count} host-switch branches, with {host_summary.certain_host_switch_count} certain and {host_summary.uncertain_host_switch_count} uncertain changes.",
                    f"Geographic reconstruction over {report.dataset.geography_trait} inferred {migration_summary.event_count} migration events across {geography_summary.changed_branch_count} changed branches.",
                    "The workflow keeps the alignment, trimmed alignment, rooted tree, support table, host tables, geographic tables, tree figure, map, and integrated report together in one durable bundle.",
                ]
            ),
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel">',
            "      <h2>Sequence-to-Tree Outputs</h2>",
            _html_list(
                [
                    "alignment: rabies-cross-host-geography-panel.aln",
                    "trimmed alignment: rabies-cross-host-geography-panel.trimmed.aln",
                    "rooted tree: rabies-cross-host-geography-panel.rooted.tree",
                    "support table: rabies-cross-host-geography-panel.support.tsv",
                    "rooting report: rabies-cross-host-geography-panel.rooting.tsv",
                ]
            ),
            _support_table(report.fasta_to_tree),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Host Switching</h2>",
            _html_list(
                [
                    f"workflow trait: {report.dataset.host_trait}",
                    f"root host confidence: {_format_number(host_summary.root_confidence)}",
                    f"host-switch rows: {len(report.host_switching.count_rows)}",
                    "see host-switch-summary.tsv, host-state-nodes.tsv, host-switch-branches.tsv, and host-switch-counts.tsv",
                ]
            ),
            _host_count_table(report.host_switching),
            "    </section>",
            '    <section class="panel full">',
            "      <h2>Biogeography</h2>",
            '      <p>The bundle includes the detailed biogeography package at <a href="biogeography/biogeography-report.html">biogeography/biogeography-report.html</a> together with the ancestral-region tree SVG and the self-contained geographic map.</p>',
            '      <div class="grid">',
            '        <div class="panel">',
            "          <h2>Ancestral-Region Tree</h2>",
            '          <div class="figure-shell">',
            '            <img src="biogeography/ancestral-region-tree.svg" alt="Ancestral region tree">',
            "          </div>",
            "        </div>",
            '        <div class="panel">',
            "          <h2>Geographic Map</h2>",
            '          <iframe src="biogeography/geographic-region-map.html" title="Geographic region map"></iframe>',
            "        </div>",
            "      </div>",
            _migration_event_table(report.biogeography_report),
            "    </section>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Key Files</h2>",
            _html_list(
                [
                    '<a href="workflow-summary.tsv">workflow-summary.tsv</a>',
                    '<a href="host-switch-summary.tsv">host-switch-summary.tsv</a>',
                    '<a href="host-switch-counts.tsv">host-switch-counts.tsv</a>',
                    '<a href="biogeography/summary.tsv">biogeography/summary.tsv</a>',
                    '<a href="biogeography/event-table.tsv">biogeography/event-table.tsv</a>',
                    '<a href="rabies-cross-host-geography.manifest.json">rabies-cross-host-geography.manifest.json</a>',
                ]
            ),
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )
    path.write_text(html + "\n", encoding="utf-8")
    return path


def _support_table(report: FastaToTreeWorkflowReport) -> str:
    return _table(
        headers=["node", "descendant_taxa", "support", "support_fraction"],
        rows=[
            [
                row.node,
                ", ".join(row.descendant_taxa),
                _format_number(row.support),
                _format_number(row.support_fraction),
            ]
            for row in report.support_rows
        ],
    )


def _host_count_table(report: HostSwitchingReport) -> str:
    return _table(
        headers=[
            "transition",
            "certain_switch_count",
            "uncertain_switch_count",
            "total_switch_count",
        ],
        rows=[
            [
                row.transition,
                str(row.certain_switch_count),
                str(row.uncertain_switch_count),
                str(row.total_switch_count),
            ]
            for row in report.count_rows
        ],
    )


def _migration_event_table(report: BiogeographyReportPackageResult) -> str:
    return _table(
        headers=[
            "branch_id",
            "source_region",
            "target_region",
            "support",
            "midpoint_depth",
        ],
        rows=[
            [
                row.branch_id,
                row.source_region,
                row.target_region,
                _format_number(row.support),
                _format_number(row.midpoint_depth),
            ]
            for row in report.event_report.event_rows
        ],
    )


def _html_list(items: list[str]) -> str:
    return "<ul>" + "".join(f"<li>{item}</li>" for item in items) + "</ul>"


def _table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{escape(header)}</th>" for header in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{escape(cell)}</td>" for cell in row) + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _support_range_text(
    minimum_support: float | None,
    maximum_support: float | None,
) -> str:
    if minimum_support is None or maximum_support is None:
        return "not available"
    return f"{_format_number(minimum_support)}-{_format_number(maximum_support)}"


def _format_number(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".12g")


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()
