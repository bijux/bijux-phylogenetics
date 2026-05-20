from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import shutil

from bijux_phylogenetics.ecology import (
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
from bijux_phylogenetics.io.fasta.records import validate_fasta_input

_DATASET_ID = "rabies_cross_host_panel"
_DATASET_LABEL = "Rabies cross-host nucleoprotein panel"
_SEQUENCE_TYPE = "dna"
_WORKFLOW_TRAIT = "host_group"
_WORKFLOW_MODEL = "ard"
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
class RabiesCrossHostPanelDataset:
    """Packaged pathogen dataset for host-switching workflow review."""

    dataset_id: str
    label: str
    dataset_root: Path
    sequences_path: Path
    tree_path: Path
    hosts_path: Path
    reference_output_root: Path
    taxon_count: int
    sequence_type: str
    workflow_trait: str
    workflow_model: str
    observed_host_group_count: int
    source_accessions: tuple[str, ...]
    source_summary: str


@dataclass(slots=True)
class RabiesCrossHostPanelExportResult:
    """Materialized copy of the packaged pathogen dataset."""

    output_root: Path
    readme_path: Path
    sequences_path: Path
    tree_path: Path
    hosts_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class RabiesCrossHostPanelWorkflowReport:
    """One governed host-switching workflow run over the packaged rabies panel."""

    dataset: RabiesCrossHostPanelDataset
    host_switching: HostSwitchingReport


@dataclass(slots=True)
class RabiesCrossHostPanelWorkflowBundle:
    """Written host-switching outputs for the packaged rabies panel."""

    output_root: Path
    analysis_constraint_mode: str
    root_host: str
    root_confidence: float
    host_switch_count: int
    certain_host_switch_count: int
    uncertain_host_switch_count: int
    workflow_summary_path: Path
    host_switch_summary_path: Path
    host_state_nodes_path: Path
    host_switch_branches_path: Path
    host_switch_counts_path: Path
    host_switch_fits_path: Path
    host_switch_unsupported_path: Path
    host_switch_exclusions_path: Path


@dataclass(slots=True)
class RabiesCrossHostPanelDemoResult:
    """Dataset export plus workflow outputs for the public pathogen demo."""

    output_root: Path
    dataset: RabiesCrossHostPanelDataset
    dataset_export: RabiesCrossHostPanelExportResult
    workflow_bundle: RabiesCrossHostPanelWorkflowBundle
    overview_path: Path


def load_rabies_cross_host_panel_dataset() -> RabiesCrossHostPanelDataset:
    """Expose the packaged rabies host-switching panel as a first-class surface."""
    dataset_root = _resource_root()
    sequences_path = dataset_root / "sequences.fasta"
    tree_path = dataset_root / "tree.nwk"
    hosts_path = dataset_root / "hosts.csv"
    validation = validate_fasta_input(sequences_path, sequence_type=_SEQUENCE_TYPE)
    observed_host_groups = _read_observed_host_groups(hosts_path)
    return RabiesCrossHostPanelDataset(
        dataset_id=_DATASET_ID,
        label=_DATASET_LABEL,
        dataset_root=dataset_root,
        sequences_path=sequences_path,
        tree_path=tree_path,
        hosts_path=hosts_path,
        reference_output_root=dataset_root / "expected",
        taxon_count=validation.summary.sequence_count,
        sequence_type=_SEQUENCE_TYPE,
        workflow_trait=_WORKFLOW_TRAIT,
        workflow_model=_WORKFLOW_MODEL,
        observed_host_group_count=len(observed_host_groups),
        source_accessions=_SOURCE_ACCESSIONS,
        source_summary=(
            "Real rabies virus nucleoprotein sequences spanning bat, canid, and "
            "livestock hosts, packaged with grouped host metadata and a rooted "
            "tree for host-switching review."
        ),
    )


def export_rabies_cross_host_panel_dataset(
    destination: Path,
) -> RabiesCrossHostPanelExportResult:
    """Copy the packaged rabies host-switching dataset and reference outputs."""
    dataset = load_rabies_cross_host_panel_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(
        dataset.dataset_root / "README.md", destination / "README.md"
    )
    sequences_path = shutil.copy2(
        dataset.sequences_path, destination / "sequences.fasta"
    )
    tree_path = shutil.copy2(dataset.tree_path, destination / "tree.nwk")
    hosts_path = shutil.copy2(dataset.hosts_path, destination / "hosts.csv")
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return RabiesCrossHostPanelExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        sequences_path=Path(sequences_path),
        tree_path=Path(tree_path),
        hosts_path=Path(hosts_path),
        expected_output_root=expected_output_root,
    )


def run_rabies_cross_host_panel_workflow() -> RabiesCrossHostPanelWorkflowReport:
    """Run the owned host-switching workflow over the packaged rabies panel."""
    dataset = load_rabies_cross_host_panel_dataset()
    host_switching = summarize_host_switching(
        dataset.tree_path,
        dataset.hosts_path,
        trait=dataset.workflow_trait,
        taxon_column="taxon",
        model=dataset.workflow_model,
    )
    return RabiesCrossHostPanelWorkflowReport(
        dataset=dataset,
        host_switching=host_switching,
    )


def write_rabies_cross_host_panel_workflow_bundle(
    output_root: Path,
    report: RabiesCrossHostPanelWorkflowReport,
) -> RabiesCrossHostPanelWorkflowBundle:
    """Write the governed host-switching outputs for the packaged rabies panel."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    summary = report.host_switching.summary
    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report,
    )
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
    return RabiesCrossHostPanelWorkflowBundle(
        output_root=output_root,
        analysis_constraint_mode=summary.analysis_constraint_mode,
        root_host=summary.root_host,
        root_confidence=summary.root_confidence,
        host_switch_count=summary.host_switch_count,
        certain_host_switch_count=summary.certain_host_switch_count,
        uncertain_host_switch_count=summary.uncertain_host_switch_count,
        workflow_summary_path=workflow_summary_path,
        host_switch_summary_path=host_switch_summary_path,
        host_state_nodes_path=host_state_nodes_path,
        host_switch_branches_path=host_switch_branches_path,
        host_switch_counts_path=host_switch_counts_path,
        host_switch_fits_path=host_switch_fits_path,
        host_switch_unsupported_path=host_switch_unsupported_path,
        host_switch_exclusions_path=host_switch_exclusions_path,
    )


def run_rabies_cross_host_panel_demo(
    output_root: Path,
) -> RabiesCrossHostPanelDemoResult:
    """Materialize the packaged pathogen dataset and rerun host-switching review."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    report = run_rabies_cross_host_panel_workflow()
    dataset_export = export_rabies_cross_host_panel_dataset(output_root / "dataset")
    workflow_bundle = write_rabies_cross_host_panel_workflow_bundle(
        output_root / "workflow",
        report,
    )
    overview_path = _write_overview(
        output_root / "overview.md", report, workflow_bundle
    )
    return RabiesCrossHostPanelDemoResult(
        output_root=output_root,
        dataset=report.dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _write_workflow_summary_table(
    path: Path,
    report: RabiesCrossHostPanelWorkflowReport,
) -> Path:
    summary = report.host_switching.summary
    rows = [
        "\t".join(
            [
                "dataset_id",
                "taxon_count",
                "workflow_trait",
                "workflow_model",
                "observed_host_group_count",
                "root_host",
                "root_confidence",
                "host_switch_count",
                "certain_host_switch_count",
                "uncertain_host_switch_count",
                "ambiguous_internal_node_count",
                "preferred_constraint",
                "warning_count",
            ]
        ),
        "\t".join(
            [
                report.dataset.dataset_id,
                str(report.dataset.taxon_count),
                report.dataset.workflow_trait,
                report.dataset.workflow_model,
                str(report.dataset.observed_host_group_count),
                summary.root_host,
                _format_number(summary.root_confidence),
                str(summary.host_switch_count),
                str(summary.certain_host_switch_count),
                str(summary.uncertain_host_switch_count),
                str(summary.ambiguous_internal_node_count),
                summary.preferred_constraint,
                str(summary.warning_count),
            ]
        ),
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _write_overview(
    path: Path,
    report: RabiesCrossHostPanelWorkflowReport,
    bundle: RabiesCrossHostPanelWorkflowBundle,
) -> Path:
    summary = report.host_switching.summary
    lines = [
        "# Rabies Cross-Host Demo",
        "",
        f"- dataset id: `{report.dataset.dataset_id}`",
        f"- taxon count: `{report.dataset.taxon_count}`",
        f"- workflow trait: `{report.dataset.workflow_trait}`",
        f"- observed host groups: `{report.dataset.observed_host_group_count}`",
        f"- root host: `{summary.root_host}`",
        f"- host switch count: `{summary.host_switch_count}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{bundle.workflow_summary_path.name}`",
        f"- host-switch summary: `{bundle.host_switch_summary_path.name}`",
        f"- host-state nodes: `{bundle.host_state_nodes_path.name}`",
        f"- host-switch branches: `{bundle.host_switch_branches_path.name}`",
        f"- host-switch counts: `{bundle.host_switch_counts_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _read_observed_host_groups(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        observed = {
            row["host_group"].strip()
            for row in reader
            if row.get("host_group") and row["host_group"].strip()
        }
    return sorted(observed)


def _format_number(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".12g")


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "resources"
        / "datasets"
        / "pathogens"
        / _DATASET_ID
    )
