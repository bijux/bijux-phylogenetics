from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .contracts import HostSwitchingReport
from .shared import format_optional_float


def write_host_switch_summary_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one overall summary ledger for host-switching analysis."""
    summary = report.summary
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "internal_model",
            "analysis_constraint_mode",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observed_host_count",
            "internal_node_count",
            "ambiguous_internal_node_count",
            "host_switch_count",
            "certain_host_switch_count",
            "uncertain_host_switch_count",
            "allowed_transition_count",
            "forbidden_transition_count",
            "constrained_log_likelihood",
            "unconstrained_log_likelihood",
            "constrained_aic",
            "unconstrained_aic",
            "preferred_constraint",
            "unsupported_switch_claim_count",
            "root_host",
            "root_confidence",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "internal_model": summary.internal_model,
                "analysis_constraint_mode": summary.analysis_constraint_mode,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "observed_host_count": str(summary.observed_host_count),
                "internal_node_count": str(summary.internal_node_count),
                "ambiguous_internal_node_count": str(
                    summary.ambiguous_internal_node_count
                ),
                "host_switch_count": str(summary.host_switch_count),
                "certain_host_switch_count": str(summary.certain_host_switch_count),
                "uncertain_host_switch_count": str(summary.uncertain_host_switch_count),
                "allowed_transition_count": str(summary.allowed_transition_count),
                "forbidden_transition_count": str(summary.forbidden_transition_count),
                "constrained_log_likelihood": format_optional_float(
                    summary.constrained_log_likelihood
                ),
                "unconstrained_log_likelihood": str(
                    summary.unconstrained_log_likelihood
                ),
                "constrained_aic": format_optional_float(summary.constrained_aic),
                "unconstrained_aic": str(summary.unconstrained_aic),
                "preferred_constraint": summary.preferred_constraint,
                "unsupported_switch_claim_count": str(
                    summary.unsupported_switch_claim_count
                ),
                "root_host": summary.root_host,
                "root_confidence": str(summary.root_confidence),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_host_state_node_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one internal-node host-state ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "most_likely_host",
            "host_probabilities",
            "confidence",
            "ambiguous",
            "is_root",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "descendant_taxa": ",".join(row.descendant_taxa),
                "most_likely_host": row.most_likely_host,
                "host_probabilities": json.dumps(
                    row.host_probabilities,
                    sort_keys=True,
                ),
                "confidence": str(row.confidence),
                "ambiguous": str(row.ambiguous).lower(),
                "is_root": str(row.is_root).lower(),
            }
            for row in report.node_rows
        ],
    )


def write_host_switch_branch_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one branchwise host-switch ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "branch_length",
            "parent_most_likely_host",
            "child_most_likely_host",
            "parent_host_set",
            "child_host_set",
            "overlapping_hosts",
            "changed",
            "transition",
            "certainty_class",
            "parent_confidence",
            "child_confidence",
            "transition_allowed",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "branch_length": format_optional_float(row.branch_length),
                "parent_most_likely_host": row.parent_most_likely_host,
                "child_most_likely_host": row.child_most_likely_host,
                "parent_host_set": ",".join(row.parent_host_set),
                "child_host_set": ",".join(row.child_host_set),
                "overlapping_hosts": ",".join(row.overlapping_hosts),
                "changed": str(row.changed).lower(),
                "transition": row.transition,
                "certainty_class": row.certainty_class,
                "parent_confidence": str(row.parent_confidence),
                "child_confidence": str(row.child_confidence),
                "transition_allowed": str(row.transition_allowed).lower(),
            }
            for row in report.branch_rows
        ],
    )


def write_host_switch_count_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one aggregated host-switch count ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "transition",
            "source_host",
            "target_host",
            "transition_allowed",
            "certain_switch_count",
            "uncertain_switch_count",
            "total_switch_count",
        ],
        rows=[
            {
                "transition": row.transition,
                "source_host": row.source_host,
                "target_host": row.target_host,
                "transition_allowed": str(row.transition_allowed).lower(),
                "certain_switch_count": str(row.certain_switch_count),
                "uncertain_switch_count": str(row.uncertain_switch_count),
                "total_switch_count": str(row.total_switch_count),
            }
            for row in report.count_rows
        ],
    )


def write_host_switch_fit_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one fit-comparison ledger for unconstrained and constrained host models."""
    return write_taxon_rows(
        path,
        columns=[
            "constraint_mode",
            "model",
            "analyzed_taxon_count",
            "log_likelihood",
            "parameter_count",
            "aic",
            "root_host",
            "root_confidence",
        ],
        rows=[
            {
                "constraint_mode": row.constraint_mode,
                "model": row.model,
                "analyzed_taxon_count": str(row.analyzed_taxon_count),
                "log_likelihood": str(row.log_likelihood),
                "parameter_count": str(row.parameter_count),
                "aic": str(row.aic),
                "root_host": row.root_host,
                "root_confidence": str(row.root_confidence),
            }
            for row in report.fit_rows
        ],
    )


def write_unsupported_host_switch_claim_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one forbidden unconstrained host-switch claim ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "unconstrained_source_host",
            "unconstrained_target_host",
            "unconstrained_certainty_class",
            "constrained_source_host",
            "constrained_target_host",
            "constrained_certainty_class",
            "claim_resolved",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "unconstrained_source_host": row.unconstrained_source_host,
                "unconstrained_target_host": row.unconstrained_target_host,
                "unconstrained_certainty_class": row.unconstrained_certainty_class,
                "constrained_source_host": row.constrained_source_host,
                "constrained_target_host": row.constrained_target_host,
                "constrained_certainty_class": row.constrained_certainty_class,
                "claim_resolved": str(row.claim_resolved).lower(),
            }
            for row in report.unsupported_claim_rows
        ],
    )


def write_host_switch_exclusion_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one excluded-row ledger for host-switching analysis."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "raw_host",
            "normalized_host",
            "reason",
            "note",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "raw_host": row.raw_host,
                "normalized_host": row.normalized_host or "",
                "reason": row.reason,
                "note": row.note,
            }
            for row in report.exclusion_rows
        ],
    )
