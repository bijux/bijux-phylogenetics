from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .contracts import NicheTransitionReport
from .shared import format_optional_float


def write_niche_transition_summary_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one overall ecological niche transition summary ledger."""
    summary = report.summary
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "internal_model",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observed_niche_count",
            "internal_node_count",
            "ambiguous_internal_node_count",
            "log_likelihood",
            "parameter_count",
            "aic",
            "transition_rate_row_count",
            "changed_branch_count",
            "certain_transition_count",
            "uncertain_transition_count",
            "strongly_supported_transition_count",
            "clade_shift_row_count",
            "repeated_shift_clade_count",
            "root_niche",
            "root_confidence",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "internal_model": summary.internal_model,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "observed_niche_count": str(summary.observed_niche_count),
                "internal_node_count": str(summary.internal_node_count),
                "ambiguous_internal_node_count": str(
                    summary.ambiguous_internal_node_count
                ),
                "log_likelihood": str(summary.log_likelihood),
                "parameter_count": str(summary.parameter_count),
                "aic": str(summary.aic),
                "transition_rate_row_count": str(summary.transition_rate_row_count),
                "changed_branch_count": str(summary.changed_branch_count),
                "certain_transition_count": str(summary.certain_transition_count),
                "uncertain_transition_count": str(summary.uncertain_transition_count),
                "strongly_supported_transition_count": str(
                    summary.strongly_supported_transition_count
                ),
                "clade_shift_row_count": str(summary.clade_shift_row_count),
                "repeated_shift_clade_count": str(summary.repeated_shift_clade_count),
                "root_niche": summary.root_niche,
                "root_confidence": str(summary.root_confidence),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_niche_state_node_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one internal-node ecological niche probability ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "most_likely_niche",
            "niche_probabilities",
            "confidence",
            "ambiguous",
            "is_root",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "descendant_taxa": ",".join(row.descendant_taxa),
                "most_likely_niche": row.most_likely_niche,
                "niche_probabilities": json.dumps(
                    row.niche_probabilities,
                    sort_keys=True,
                ),
                "confidence": str(row.confidence),
                "ambiguous": str(row.ambiguous).lower(),
                "is_root": str(row.is_root).lower(),
            }
            for row in report.node_rows
        ],
    )


def write_niche_transition_rate_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one ecological niche transition-rate ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "source_niche",
            "target_niche",
            "transition_allowed",
            "step_distance",
            "rate",
        ],
        rows=[
            {
                "source_niche": row.source_niche,
                "target_niche": row.target_niche,
                "transition_allowed": str(row.transition_allowed).lower(),
                "step_distance": str(row.step_distance),
                "rate": str(row.rate),
            }
            for row in report.rate_rows
        ],
    )


def write_niche_transition_branch_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one branchwise ecological niche transition ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "branch_length",
            "parent_most_likely_niche",
            "child_most_likely_niche",
            "parent_niche_set",
            "child_niche_set",
            "overlapping_niches",
            "changed",
            "transition",
            "certainty_class",
            "support",
            "strongly_supported",
            "parent_confidence",
            "child_confidence",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "branch_length": format_optional_float(row.branch_length),
                "parent_most_likely_niche": row.parent_most_likely_niche,
                "child_most_likely_niche": row.child_most_likely_niche,
                "parent_niche_set": ",".join(row.parent_niche_set),
                "child_niche_set": ",".join(row.child_niche_set),
                "overlapping_niches": ",".join(row.overlapping_niches),
                "changed": str(row.changed).lower(),
                "transition": row.transition,
                "certainty_class": row.certainty_class,
                "support": str(row.support),
                "strongly_supported": str(row.strongly_supported).lower(),
                "parent_confidence": str(row.parent_confidence),
                "child_confidence": str(row.child_confidence),
            }
            for row in report.branch_rows
        ],
    )


def write_niche_transition_count_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one aggregated ecological niche transition-count ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "transition",
            "source_niche",
            "target_niche",
            "certain_transition_count",
            "uncertain_transition_count",
            "total_transition_count",
            "strongly_supported_transition_count",
        ],
        rows=[
            {
                "transition": row.transition,
                "source_niche": row.source_niche,
                "target_niche": row.target_niche,
                "certain_transition_count": str(row.certain_transition_count),
                "uncertain_transition_count": str(row.uncertain_transition_count),
                "total_transition_count": str(row.total_transition_count),
                "strongly_supported_transition_count": str(
                    row.strongly_supported_transition_count
                ),
            }
            for row in report.count_rows
        ],
    )


def write_niche_transition_clade_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one internal-clade ecological niche shift ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "descendant_taxon_count",
            "descendant_internal_node_count",
            "changed_branch_count",
            "certain_transition_count",
            "uncertain_transition_count",
            "strongly_supported_transition_count",
            "transition_diversity",
            "dominant_transition",
            "dominant_transition_count",
            "shift_burden_score",
            "contains_repeated_shifts",
            "rank",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "descendant_taxa": ",".join(row.descendant_taxa),
                "descendant_taxon_count": str(row.descendant_taxon_count),
                "descendant_internal_node_count": str(
                    row.descendant_internal_node_count
                ),
                "changed_branch_count": str(row.changed_branch_count),
                "certain_transition_count": str(row.certain_transition_count),
                "uncertain_transition_count": str(row.uncertain_transition_count),
                "strongly_supported_transition_count": str(
                    row.strongly_supported_transition_count
                ),
                "transition_diversity": str(row.transition_diversity),
                "dominant_transition": row.dominant_transition,
                "dominant_transition_count": str(row.dominant_transition_count),
                "shift_burden_score": str(row.shift_burden_score),
                "contains_repeated_shifts": str(row.contains_repeated_shifts).lower(),
                "rank": str(row.rank),
            }
            for row in report.clade_rows
        ],
    )


def write_niche_transition_exclusion_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one excluded ecological niche row ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "raw_niche",
            "normalized_niche",
            "reason",
            "note",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "raw_niche": row.raw_niche,
                "normalized_niche": row.normalized_niche or "",
                "reason": row.reason,
                "note": row.note,
            }
            for row in report.exclusion_rows
        ],
    )
