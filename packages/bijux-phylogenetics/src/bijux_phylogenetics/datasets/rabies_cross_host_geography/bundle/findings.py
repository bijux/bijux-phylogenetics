from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.reporting.analysis_package import (
    ComparativeAnalysisSummaryRow,
    ComparativeInterpretationRow,
)
from bijux_phylogenetics.compare.presentation import ComparisonReportBuildResult
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..models import (
    RabiesCrossHostGeographyPanelWorkflowReport,
    RabiesScientificFindingRow,
)
from ..shared import _format_number


def _build_scientific_finding_rows(
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    comparative_interpretation_rows: list[ComparativeInterpretationRow],
) -> list[RabiesScientificFindingRow]:
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary
    bootstrap_question = (
        "Does the bootstrap consensus preserve the rooted ML conclusion?"
    )
    if bootstrap_tree_comparison_report.topology.topology_equal:
        bootstrap_claim = "The bootstrap consensus preserves the rooted ML topology on the shared taxon set."
    else:
        bootstrap_claim = "The bootstrap consensus differs from the rooted ML topology after support-driven summarization."
    comparative_claim = next(
        (
            row.claim
            for row in comparative_interpretation_rows
            if row.topic == "coefficient" and row.claim
        ),
        "The comparative layer did not expose one stable host-associated longitude shift.",
    )
    return [
        RabiesScientificFindingRow(
            finding_id="root_host_state",
            question="What host state anchors the rooted rabies panel?",
            claim=f"The rooted tree places the ancestral host state in {host_summary.root_host}.",
            evidence=(
                f"root host confidence {_format_number(host_summary.root_confidence)} "
                f"with outgroup {','.join(report.dataset.outgroup_taxa)}"
            ),
            caution=(
                "The panel is compact and grouped by broad host classes rather than species-level host states."
            ),
            source_artifact=report.dataset.workflow_prefix + ".rooting.tsv",
        ),
        RabiesScientificFindingRow(
            finding_id="root_region_state",
            question="What geographic regime anchors the rooted rabies panel?",
            claim=f"The rooted tree places the ancestral region in {geography_summary.root_region}.",
            evidence=(
                f"root region probability {_format_number(geography_summary.root_region_probability)} "
                f"across {geography_summary.changed_branch_count} changed branches"
            ),
            caution=(
                "Grouped macroregions simplify the raw locality labels so the result should be treated as regional rather than site-level history."
            ),
            source_artifact="biogeography/summary.tsv",
        ),
        RabiesScientificFindingRow(
            finding_id="host_switching",
            question="How much host-switching signal appears in the rooted tree?",
            claim=(
                f"The host reconstruction inferred {host_summary.host_switch_count} host-switch branch changes."
            ),
            evidence=(
                f"certain changes {host_summary.certain_host_switch_count}; "
                f"uncertain changes {host_summary.uncertain_host_switch_count}"
            ),
            caution=(
                "Branch-wise host changes depend on the grouped host coding and should not be over-read as one exhaustive host-jump catalogue."
            ),
            source_artifact="host-switch-summary.tsv",
        ),
        RabiesScientificFindingRow(
            finding_id="bootstrap_consensus",
            question=bootstrap_question,
            claim=bootstrap_claim,
            evidence=(
                f"rooted RF distance {bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance}; "
                f"high-support conflicts "
                f"{len([row for row in bootstrap_tree_comparison_report.support.conflicting_clades if row.conflict_classification == 'high_support_conflict'])}"
            ),
            caution=(
                "Consensus trees can collapse low-support branches, so exact rooted agreement is stricter than shared major clades."
            ),
            source_artifact=(
                "bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv"
            ),
        ),
        RabiesScientificFindingRow(
            finding_id="comparative_longitude",
            question=(
                "Do host-associated rabies lineages occupy one distinct longitudinal regime in this panel?"
            ),
            claim=comparative_claim,
            evidence=(
                f"selected model {comparative_summary_row.selected_model}; "
                f"PGLS lambda {_format_number(comparative_summary_row.pgls_lambda)}; "
                f"r-squared {_format_number(comparative_summary_row.pgls_r_squared)}"
            ),
            caution=(
                "The comparative claim is associational, uses only nine taxa, and retains residual-diagnostic cautions."
            ),
            source_artifact="comparative/interpretation-table.tsv",
        ),
        RabiesScientificFindingRow(
            finding_id="migration_events",
            question="How much regional movement is implied by the geographic reconstruction?",
            claim=(
                f"The biogeography layer inferred {migration_summary.event_count} migration events across the rooted tree."
            ),
            evidence=(
                f"strongly supported migration events {migration_summary.strongly_supported_event_count}"
            ),
            caution=(
                "Event counts summarize transitions over grouped regions and do not replace one dated dispersal analysis."
            ),
            source_artifact="biogeography/event-table.tsv",
        ),
    ]


def _write_scientific_findings_table(
    path: Path,
    rows: list[RabiesScientificFindingRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "finding_id",
            "question",
            "claim",
            "evidence",
            "caution",
            "source_artifact",
        ],
        rows=[
            {
                "finding_id": row.finding_id,
                "question": row.question,
                "claim": row.claim,
                "evidence": row.evidence,
                "caution": row.caution,
                "source_artifact": row.source_artifact,
            }
            for row in rows
        ],
    )
