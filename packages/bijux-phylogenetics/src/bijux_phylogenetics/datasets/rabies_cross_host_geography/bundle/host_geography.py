from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil

from bijux_phylogenetics.datasets.rabies_cross_host_geography.models import (
    RabiesCrossHostGeographyPanelWorkflowReport,
)
from bijux_phylogenetics.ecology import (
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_count_table,
    write_host_switch_exclusion_table,
    write_host_switch_fit_table,
    write_host_switch_summary_table,
    write_unsupported_host_switch_claim_table,
)


@dataclass(frozen=True)
class HostGeographyArtifacts:
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


def _write_host_geography_artifacts(
    output_root: Path,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
) -> HostGeographyArtifacts:
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
    return HostGeographyArtifacts(
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
    )
