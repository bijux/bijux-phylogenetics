from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.rabies_cross_host_geography.models import (
    RabiesCrossHostGeographyPanelWorkflowReport,
)
from bijux_phylogenetics.diagnostics.conclusion_stability import (
    write_ancestral_state_stability_table,
    write_comparative_coefficient_stability_table,
    write_conclusion_stability_report_html,
    write_conclusion_stability_summary_table,
    write_key_clade_stability_table,
    write_support_value_stability_table,
)


@dataclass(frozen=True)
class ConclusionStabilityArtifacts:
    output_root: Path
    summary_path: Path
    key_clade_stability_path: Path
    support_value_stability_path: Path
    ancestral_state_stability_path: Path
    comparative_coefficient_stability_path: Path
    report_path: Path


def _write_conclusion_stability_artifacts(
    output_root: Path,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
) -> ConclusionStabilityArtifacts:
    conclusion_stability_output_root = output_root / "conclusion-stability"
    conclusion_stability_output_root.mkdir(parents=True, exist_ok=True)
    conclusion_stability_summary_path = write_conclusion_stability_summary_table(
        conclusion_stability_output_root / "conclusion-stability-summary.tsv",
        report.conclusion_stability_report,
    )
    key_clade_stability_path = write_key_clade_stability_table(
        conclusion_stability_output_root / "key-clade-stability.tsv",
        report.conclusion_stability_report.key_clade_rows,
    )
    support_value_stability_path = write_support_value_stability_table(
        conclusion_stability_output_root / "support-value-stability.tsv",
        report.conclusion_stability_report.support_value_rows,
    )
    ancestral_state_stability_path = write_ancestral_state_stability_table(
        conclusion_stability_output_root / "ancestral-state-stability.tsv",
        report.conclusion_stability_report.ancestral_state_rows,
    )
    comparative_coefficient_stability_path = (
        write_comparative_coefficient_stability_table(
            conclusion_stability_output_root / "comparative-coefficient-stability.tsv",
            report.conclusion_stability_report.comparative_coefficient_rows,
        )
    )
    conclusion_stability_report_path = write_conclusion_stability_report_html(
        conclusion_stability_output_root / "conclusion-stability-report.html",
        report.conclusion_stability_report,
    )
    return ConclusionStabilityArtifacts(
        output_root=conclusion_stability_output_root,
        summary_path=conclusion_stability_summary_path,
        key_clade_stability_path=key_clade_stability_path,
        support_value_stability_path=support_value_stability_path,
        ancestral_state_stability_path=ancestral_state_stability_path,
        comparative_coefficient_stability_path=comparative_coefficient_stability_path,
        report_path=conclusion_stability_report_path,
    )
