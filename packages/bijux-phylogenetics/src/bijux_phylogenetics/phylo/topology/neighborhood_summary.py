from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from bijux_phylogenetics.phylo.topology.models import (
    RootedNniNeighborhoodReport,
    RootedSprNeighborhoodReport,
    RootedTbrNeighborhoodReport,
    TopologyNeighborhoodSummaryReport,
)


def summarize_topology_neighborhood(
    report: (
        RootedNniNeighborhoodReport
        | RootedSprNeighborhoodReport
        | RootedTbrNeighborhoodReport
    ),
) -> TopologyNeighborhoodSummaryReport:
    """Project one rooted topology-neighborhood report into a comparable summary row."""
    if isinstance(report, RootedNniNeighborhoodReport):
        candidate_count = report.expected_neighbor_count
        valid_count = report.generated_neighbor_count
        duplicate_count = max(0, candidate_count - valid_count)
        return TopologyNeighborhoodSummaryReport(
            neighborhood_family="rooted-nni",
            algorithm=report.algorithm,
            candidate_count=candidate_count,
            valid_count=valid_count,
            duplicate_count=duplicate_count,
            skipped_count=0,
            skipped_reason="none",
            budget_reason="unbounded",
        )
    if isinstance(report, RootedSprNeighborhoodReport):
        candidate_count = (
            report.generated_move_candidate_count
            + report.skipped_budget_move_candidate_count
        )
        skipped_count = (
            report.identity_move_candidate_count
            + report.self_regraft_candidate_count
            + report.skipped_budget_move_candidate_count
        )
        duplicate_count = max(
            0,
            candidate_count - report.generated_neighbor_count - skipped_count,
        )
        return TopologyNeighborhoodSummaryReport(
            neighborhood_family="rooted-spr",
            algorithm=report.algorithm,
            candidate_count=candidate_count,
            valid_count=report.generated_neighbor_count,
            duplicate_count=duplicate_count,
            skipped_count=skipped_count,
            skipped_reason=_summarize_rooted_spr_skipped_reason(report),
            budget_reason=_summarize_rooted_spr_budget_reason(report),
        )
    candidate_count = report.generated_reconnection_count
    skipped_count = report.identity_reconnection_count
    duplicate_count = max(
        0,
        candidate_count - report.generated_neighbor_count - skipped_count,
    )
    return TopologyNeighborhoodSummaryReport(
        neighborhood_family="rooted-tbr",
        algorithm=report.algorithm,
        candidate_count=candidate_count,
        valid_count=report.generated_neighbor_count,
        duplicate_count=duplicate_count,
        skipped_count=skipped_count,
        skipped_reason=_summarize_rooted_tbr_skipped_reason(report),
        budget_reason="unbounded",
    )


def write_topology_neighborhood_summary_table(
    path: Path,
    reports: Sequence[TopologyNeighborhoodSummaryReport],
) -> Path:
    """Write comparable topology-neighborhood summary rows as a TSV table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "neighborhood_family",
                "algorithm",
                "candidate_count",
                "valid_count",
                "duplicate_count",
                "skipped_count",
                "skipped_reason",
                "budget_reason",
            ]
        )
    ]
    for report in reports:
        lines.append(
            "\t".join(
                [
                    report.neighborhood_family,
                    report.algorithm,
                    str(report.candidate_count),
                    str(report.valid_count),
                    str(report.duplicate_count),
                    str(report.skipped_count),
                    report.skipped_reason,
                    report.budget_reason,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _summarize_rooted_spr_skipped_reason(
    report: RootedSprNeighborhoodReport,
) -> str:
    reasons: list[str] = []
    if report.identity_move_candidate_count:
        reasons.append(f"identity-candidate:{report.identity_move_candidate_count}")
    if report.self_regraft_candidate_count:
        reasons.append(f"self-regraft-candidate:{report.self_regraft_candidate_count}")
    if report.skipped_budget_move_candidate_count:
        reasons.append(f"budget-candidate:{report.skipped_budget_move_candidate_count}")
    return "none" if not reasons else ";".join(reasons)


def _summarize_rooted_spr_budget_reason(
    report: RootedSprNeighborhoodReport,
) -> str:
    reasons: list[str] = []
    if report.max_pruned_clade_count is not None:
        reasons.append(f"max-pruned-clades={report.max_pruned_clade_count}")
    if report.max_regraft_target_count_per_pruned_clade is not None:
        reasons.append(
            "max-regraft-targets-per-pruned-clade="
            f"{report.max_regraft_target_count_per_pruned_clade}"
        )
    return "unbounded" if not reasons else ";".join(reasons)


def _summarize_rooted_tbr_skipped_reason(
    report: RootedTbrNeighborhoodReport,
) -> str:
    if report.identity_reconnection_count == 0:
        return "none"
    return f"identity-candidate:{report.identity_reconnection_count}"
