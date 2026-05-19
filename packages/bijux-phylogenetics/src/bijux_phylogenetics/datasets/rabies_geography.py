from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import shutil

from bijux_phylogenetics.biogeography import (
    GeographicMigrationEventReport,
    GeographicStateModelReport,
    summarize_geographic_migration_events,
    summarize_geographic_state_model,
    write_geographic_exclusion_table,
    write_geographic_migration_event_summary_table,
    write_geographic_migration_event_table,
    write_geographic_migration_exclusion_table,
    write_geographic_region_probability_table,
    write_geographic_state_summary_table,
    write_geographic_transition_event_table,
    write_geographic_transition_rate_table,
)
from bijux_phylogenetics.io.fasta.records import validate_fasta_input

_DATASET_ID = "rabies_geographic_transition_panel"
_DATASET_LABEL = "Rabies geographic transition panel"
_SEQUENCE_TYPE = "dna"
_WORKFLOW_TRAIT = "region_group"
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
class RabiesGeographicTransitionPanelDataset:
    """Packaged pathogen geography dataset for biogeographic transition review."""

    dataset_id: str
    label: str
    dataset_root: Path
    sequences_path: Path
    tree_path: Path
    regions_path: Path
    reference_output_root: Path
    taxon_count: int
    sequence_type: str
    workflow_trait: str
    workflow_model: str
    observed_region_group_count: int
    source_accessions: tuple[str, ...]
    source_summary: str


@dataclass(slots=True)
class RabiesGeographicTransitionPanelExportResult:
    """Materialized copy of the packaged pathogen geography dataset."""

    output_root: Path
    readme_path: Path
    sequences_path: Path
    tree_path: Path
    regions_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class RabiesGeographicTransitionPanelWorkflowReport:
    """One governed geography-transition workflow run over the packaged rabies panel."""

    dataset: RabiesGeographicTransitionPanelDataset
    geographic_state_model: GeographicStateModelReport
    migration_events: GeographicMigrationEventReport


@dataclass(slots=True)
class RabiesGeographicTransitionPanelWorkflowBundle:
    """Written geography-transition outputs for the packaged rabies panel."""

    output_root: Path
    root_region: str
    root_region_probability: float
    changed_branch_count: int
    strongly_supported_transition_count: int
    migration_event_count: int
    strongly_supported_migration_event_count: int
    workflow_summary_path: Path
    geographic_state_summary_path: Path
    geographic_region_probability_path: Path
    geographic_transition_rate_path: Path
    geographic_transition_event_path: Path
    geographic_state_exclusion_path: Path
    geographic_migration_summary_path: Path
    geographic_migration_event_path: Path
    geographic_migration_exclusion_path: Path


@dataclass(slots=True)
class RabiesGeographicTransitionPanelDemoResult:
    """Dataset export plus geography-transition outputs for the public rabies demo."""

    output_root: Path
    dataset: RabiesGeographicTransitionPanelDataset
    dataset_export: RabiesGeographicTransitionPanelExportResult
    workflow_bundle: RabiesGeographicTransitionPanelWorkflowBundle
    overview_path: Path


def load_rabies_geographic_transition_panel_dataset() -> (
    RabiesGeographicTransitionPanelDataset
):
    """Expose the packaged rabies geography panel as a first-class surface."""
    dataset_root = _resource_root()
    sequences_path = dataset_root / "sequences.fasta"
    tree_path = dataset_root / "tree.nwk"
    regions_path = dataset_root / "regions.csv"
    validation = validate_fasta_input(sequences_path, sequence_type=_SEQUENCE_TYPE)
    observed_region_groups = _read_observed_region_groups(regions_path)
    return RabiesGeographicTransitionPanelDataset(
        dataset_id=_DATASET_ID,
        label=_DATASET_LABEL,
        dataset_root=dataset_root,
        sequences_path=sequences_path,
        tree_path=tree_path,
        regions_path=regions_path,
        reference_output_root=dataset_root / "expected",
        taxon_count=validation.summary.sequence_count,
        sequence_type=_SEQUENCE_TYPE,
        workflow_trait=_WORKFLOW_TRAIT,
        workflow_model=_WORKFLOW_MODEL,
        observed_region_group_count=len(observed_region_groups),
        source_accessions=_SOURCE_ACCESSIONS,
        source_summary=(
            "Real rabies virus nucleoprotein sequences spanning South America, "
            "North America, Europe, South Asia, and northern Asia, packaged "
            "with grouped region metadata and a rooted tree for geographic "
            "transition review."
        ),
    )


def export_rabies_geographic_transition_panel_dataset(
    destination: Path,
) -> RabiesGeographicTransitionPanelExportResult:
    """Copy the packaged rabies geography dataset and reference outputs."""
    dataset = load_rabies_geographic_transition_panel_dataset()
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
    regions_path = shutil.copy2(dataset.regions_path, destination / "regions.csv")
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return RabiesGeographicTransitionPanelExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        sequences_path=Path(sequences_path),
        tree_path=Path(tree_path),
        regions_path=Path(regions_path),
        expected_output_root=expected_output_root,
    )


def run_rabies_geographic_transition_panel_workflow() -> (
    RabiesGeographicTransitionPanelWorkflowReport
):
    """Run the owned biogeography workflow over the packaged rabies panel."""
    dataset = load_rabies_geographic_transition_panel_dataset()
    geographic_state_model = summarize_geographic_state_model(
        dataset.tree_path,
        dataset.regions_path,
        trait=dataset.workflow_trait,
        taxon_column="taxon",
        model=dataset.workflow_model,
    )
    migration_events = summarize_geographic_migration_events(
        dataset.tree_path,
        dataset.regions_path,
        trait=dataset.workflow_trait,
        taxon_column="taxon",
        model=dataset.workflow_model,
    )
    return RabiesGeographicTransitionPanelWorkflowReport(
        dataset=dataset,
        geographic_state_model=geographic_state_model,
        migration_events=migration_events,
    )


def write_rabies_geographic_transition_panel_workflow_bundle(
    output_root: Path,
    report: RabiesGeographicTransitionPanelWorkflowReport,
) -> RabiesGeographicTransitionPanelWorkflowBundle:
    """Write the governed geography-transition outputs for the rabies panel."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    state_summary = report.geographic_state_model.summary
    migration_summary = report.migration_events.summary
    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report,
    )
    geographic_state_summary_path = write_geographic_state_summary_table(
        output_root / "geographic-state-summary.tsv",
        report.geographic_state_model,
    )
    geographic_region_probability_path = write_geographic_region_probability_table(
        output_root / "geographic-region-probabilities.tsv",
        report.geographic_state_model,
    )
    geographic_transition_rate_path = write_geographic_transition_rate_table(
        output_root / "geographic-transition-rates.tsv",
        report.geographic_state_model,
    )
    geographic_transition_event_path = write_geographic_transition_event_table(
        output_root / "geographic-transition-events.tsv",
        report.geographic_state_model,
    )
    geographic_state_exclusion_path = write_geographic_exclusion_table(
        output_root / "geographic-state-exclusions.tsv",
        report.geographic_state_model,
    )
    geographic_migration_summary_path = write_geographic_migration_event_summary_table(
        output_root / "geographic-migration-summary.tsv",
        report.migration_events,
    )
    geographic_migration_event_path = write_geographic_migration_event_table(
        output_root / "geographic-migration-events.tsv",
        report.migration_events,
    )
    geographic_migration_exclusion_path = write_geographic_migration_exclusion_table(
        output_root / "geographic-migration-exclusions.tsv",
        report.migration_events,
    )
    return RabiesGeographicTransitionPanelWorkflowBundle(
        output_root=output_root,
        root_region=state_summary.root_region,
        root_region_probability=state_summary.root_region_probability,
        changed_branch_count=state_summary.changed_branch_count,
        strongly_supported_transition_count=(
            state_summary.strongly_supported_transition_count
        ),
        migration_event_count=migration_summary.event_count,
        strongly_supported_migration_event_count=(
            migration_summary.strongly_supported_event_count
        ),
        workflow_summary_path=workflow_summary_path,
        geographic_state_summary_path=geographic_state_summary_path,
        geographic_region_probability_path=geographic_region_probability_path,
        geographic_transition_rate_path=geographic_transition_rate_path,
        geographic_transition_event_path=geographic_transition_event_path,
        geographic_state_exclusion_path=geographic_state_exclusion_path,
        geographic_migration_summary_path=geographic_migration_summary_path,
        geographic_migration_event_path=geographic_migration_event_path,
        geographic_migration_exclusion_path=geographic_migration_exclusion_path,
    )


def run_rabies_geographic_transition_panel_demo(
    output_root: Path,
) -> RabiesGeographicTransitionPanelDemoResult:
    """Materialize the packaged rabies geography dataset and rerun the workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    report = run_rabies_geographic_transition_panel_workflow()
    dataset_export = export_rabies_geographic_transition_panel_dataset(
        output_root / "dataset"
    )
    workflow_bundle = write_rabies_geographic_transition_panel_workflow_bundle(
        output_root / "workflow",
        report,
    )
    overview_path = _write_overview(
        output_root / "overview.md", report, workflow_bundle
    )
    return RabiesGeographicTransitionPanelDemoResult(
        output_root=output_root,
        dataset=report.dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _write_workflow_summary_table(
    path: Path,
    report: RabiesGeographicTransitionPanelWorkflowReport,
) -> Path:
    state_summary = report.geographic_state_model.summary
    migration_summary = report.migration_events.summary
    rows = [
        "\t".join(
            [
                "dataset_id",
                "taxon_count",
                "workflow_trait",
                "workflow_model",
                "observed_region_group_count",
                "root_region",
                "root_region_probability",
                "changed_branch_count",
                "strongly_supported_transition_count",
                "migration_event_count",
                "strongly_supported_migration_event_count",
                "warning_count",
            ]
        ),
        "\t".join(
            [
                report.dataset.dataset_id,
                str(report.dataset.taxon_count),
                report.dataset.workflow_trait,
                report.dataset.workflow_model,
                str(report.dataset.observed_region_group_count),
                state_summary.root_region,
                _format_number(state_summary.root_region_probability),
                str(state_summary.changed_branch_count),
                str(state_summary.strongly_supported_transition_count),
                str(migration_summary.event_count),
                str(migration_summary.strongly_supported_event_count),
                str(state_summary.warning_count),
            ]
        ),
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _write_overview(
    path: Path,
    report: RabiesGeographicTransitionPanelWorkflowReport,
    bundle: RabiesGeographicTransitionPanelWorkflowBundle,
) -> Path:
    state_summary = report.geographic_state_model.summary
    migration_summary = report.migration_events.summary
    lines = [
        "# Rabies Geographic Transition Demo",
        "",
        f"- dataset id: `{report.dataset.dataset_id}`",
        f"- taxon count: `{report.dataset.taxon_count}`",
        f"- workflow trait: `{report.dataset.workflow_trait}`",
        f"- observed region groups: `{report.dataset.observed_region_group_count}`",
        f"- root region: `{state_summary.root_region}`",
        f"- changed branches: `{state_summary.changed_branch_count}`",
        f"- migration event count: `{migration_summary.event_count}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{bundle.workflow_summary_path.name}`",
        f"- geographic-state summary: `{bundle.geographic_state_summary_path.name}`",
        f"- geographic-region probabilities: `{bundle.geographic_region_probability_path.name}`",
        f"- geographic-transition rates: `{bundle.geographic_transition_rate_path.name}`",
        f"- geographic-transition events: `{bundle.geographic_transition_event_path.name}`",
        f"- geographic-migration summary: `{bundle.geographic_migration_summary_path.name}`",
        f"- geographic-migration events: `{bundle.geographic_migration_event_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _read_observed_region_groups(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        observed = {
            row["region_group"].strip()
            for row in reader
            if row.get("region_group") and row["region_group"].strip()
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
