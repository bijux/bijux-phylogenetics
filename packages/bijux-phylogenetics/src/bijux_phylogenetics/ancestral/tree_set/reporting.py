from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.ancestral.common import write_ancestral_rows

from .models import (
    ContinuousAncestralTreeSetReport,
    ContinuousAncestralTreeSetSummary,
    DiscreteAncestralTreeSetReport,
    DiscreteAncestralTreeSetSummary,
)


def summarize_continuous_ancestral_tree_set_report(
    report: ContinuousAncestralTreeSetReport,
) -> ContinuousAncestralTreeSetSummary:
    """Summarize the main review facts for one continuous ancestral tree-set report."""
    unstable_clades = [
        row for row in report.clade_summaries if row.stability_class != "stable"
    ]
    top_unstable = (
        max(
            unstable_clades,
            key=lambda row: (row.instability_score, row.clade_id),
        ).clade_id
        if unstable_clades
        else None
    )
    return ContinuousAncestralTreeSetSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        alpha=report.alpha,
        total_tree_count=report.total_tree_count,
        burnin_tree_count=report.burnin_tree_count,
        kept_tree_count=report.kept_tree_count,
        shared_tree_taxon_count=len(report.shared_tree_taxa),
        analysis_taxon_count=len(report.analysis_taxa),
        rooted_topology_count=report.rooted_topology_count,
        unrooted_topology_count=report.unrooted_topology_count,
        clade_summary_count=len(report.clade_summaries),
        unstable_clade_count=len(unstable_clades),
        top_unstable_clade=top_unstable,
        warning_count=len(report.warnings),
    )


def summarize_discrete_ancestral_tree_set_report(
    report: DiscreteAncestralTreeSetReport,
) -> DiscreteAncestralTreeSetSummary:
    """Summarize the main review facts for one discrete ancestral tree-set report."""
    unstable_clades = [
        row for row in report.clade_summaries if row.stability_class != "stable"
    ]
    top_unstable = (
        max(
            unstable_clades,
            key=lambda row: (row.instability_score, row.clade_id),
        ).clade_id
        if unstable_clades
        else None
    )
    return DiscreteAncestralTreeSetSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        state_ordering=report.state_ordering,
        total_tree_count=report.total_tree_count,
        burnin_tree_count=report.burnin_tree_count,
        kept_tree_count=report.kept_tree_count,
        shared_tree_taxon_count=len(report.shared_tree_taxa),
        analysis_taxon_count=len(report.analysis_taxa),
        rooted_topology_count=report.rooted_topology_count,
        unrooted_topology_count=report.unrooted_topology_count,
        clade_summary_count=len(report.clade_summaries),
        unstable_clade_count=len(unstable_clades),
        top_unstable_clade=top_unstable,
        warning_count=len(report.warnings),
    )


def write_ancestral_tree_set_tree_table(
    path: Path,
    report: ContinuousAncestralTreeSetReport | DiscreteAncestralTreeSetReport,
) -> Path:
    """Write one retained-tree ledger for an ancestral tree-set analysis."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
            "internal_clade_count",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
                "internal_clade_count": str(row.internal_clade_count),
            }
            for row in report.tree_rows
        ],
    )


def write_ancestral_tree_set_exclusion_table(
    path: Path,
    report: ContinuousAncestralTreeSetReport | DiscreteAncestralTreeSetReport,
) -> Path:
    """Write one excluded-taxa ledger for an ancestral tree-set analysis."""
    return write_ancestral_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
            }
            for row in report.exclusions
        ],
    )


def write_continuous_ancestral_tree_set_summary_table(
    path: Path,
    report: ContinuousAncestralTreeSetReport,
) -> Path:
    """Write one overall summary ledger for continuous ancestral tree-set analysis."""
    summary = summarize_continuous_ancestral_tree_set_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "alpha",
            "total_tree_count",
            "burnin_tree_count",
            "kept_tree_count",
            "shared_tree_taxon_count",
            "analysis_taxon_count",
            "rooted_topology_count",
            "unrooted_topology_count",
            "clade_summary_count",
            "unstable_clade_count",
            "top_unstable_clade",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "alpha": str(summary.alpha),
                "total_tree_count": str(summary.total_tree_count),
                "burnin_tree_count": str(summary.burnin_tree_count),
                "kept_tree_count": str(summary.kept_tree_count),
                "shared_tree_taxon_count": str(summary.shared_tree_taxon_count),
                "analysis_taxon_count": str(summary.analysis_taxon_count),
                "rooted_topology_count": str(summary.rooted_topology_count),
                "unrooted_topology_count": str(summary.unrooted_topology_count),
                "clade_summary_count": str(summary.clade_summary_count),
                "unstable_clade_count": str(summary.unstable_clade_count),
                "top_unstable_clade": summary.top_unstable_clade or "",
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_continuous_ancestral_tree_set_node_table(
    path: Path,
    report: ContinuousAncestralTreeSetReport,
) -> Path:
    """Write one per-tree internal-node ledger for continuous ancestral tree-set analysis."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
            "clade_id",
            "clade_taxa",
            "estimate",
            "standard_error",
            "lower_95_interval",
            "upper_95_interval",
            "confidence",
            "unstable",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
                "clade_id": row.clade_id,
                "clade_taxa": ",".join(row.clade_taxa),
                "estimate": str(row.estimate),
                "standard_error": str(row.standard_error),
                "lower_95_interval": str(row.lower_95_interval),
                "upper_95_interval": str(row.upper_95_interval),
                "confidence": str(row.confidence),
                "unstable": str(row.unstable).lower(),
            }
            for row in report.node_rows
        ],
    )


def write_continuous_ancestral_tree_set_clade_table(
    path: Path,
    report: ContinuousAncestralTreeSetReport,
) -> Path:
    """Write one comparable-clade summary ledger for continuous ancestral tree-set analysis."""
    return write_ancestral_rows(
        path,
        columns=[
            "clade_id",
            "clade_taxa",
            "tree_presence_count",
            "tree_presence_fraction",
            "mean_estimate",
            "median_estimate",
            "standard_deviation",
            "minimum_estimate",
            "maximum_estimate",
            "lower_95_empirical_estimate",
            "upper_95_empirical_estimate",
            "empirical_interval_width",
            "mean_standard_error",
            "unstable_tree_count",
            "unstable_tree_fraction",
            "instability_score",
            "stability_class",
        ],
        rows=[
            {
                "clade_id": row.clade_id,
                "clade_taxa": ",".join(row.clade_taxa),
                "tree_presence_count": str(row.tree_presence_count),
                "tree_presence_fraction": str(row.tree_presence_fraction),
                "mean_estimate": str(row.mean_estimate),
                "median_estimate": str(row.median_estimate),
                "standard_deviation": str(row.standard_deviation),
                "minimum_estimate": str(row.minimum_estimate),
                "maximum_estimate": str(row.maximum_estimate),
                "lower_95_empirical_estimate": str(row.lower_95_empirical_estimate),
                "upper_95_empirical_estimate": str(row.upper_95_empirical_estimate),
                "empirical_interval_width": str(row.empirical_interval_width),
                "mean_standard_error": str(row.mean_standard_error),
                "unstable_tree_count": str(row.unstable_tree_count),
                "unstable_tree_fraction": str(row.unstable_tree_fraction),
                "instability_score": str(row.instability_score),
                "stability_class": row.stability_class,
            }
            for row in report.clade_summaries
        ],
    )


def write_discrete_ancestral_tree_set_summary_table(
    path: Path,
    report: DiscreteAncestralTreeSetReport,
) -> Path:
    """Write one overall summary ledger for discrete ancestral tree-set analysis."""
    summary = summarize_discrete_ancestral_tree_set_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "state_ordering",
            "total_tree_count",
            "burnin_tree_count",
            "kept_tree_count",
            "shared_tree_taxon_count",
            "analysis_taxon_count",
            "rooted_topology_count",
            "unrooted_topology_count",
            "clade_summary_count",
            "unstable_clade_count",
            "top_unstable_clade",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "state_ordering": summary.state_ordering,
                "total_tree_count": str(summary.total_tree_count),
                "burnin_tree_count": str(summary.burnin_tree_count),
                "kept_tree_count": str(summary.kept_tree_count),
                "shared_tree_taxon_count": str(summary.shared_tree_taxon_count),
                "analysis_taxon_count": str(summary.analysis_taxon_count),
                "rooted_topology_count": str(summary.rooted_topology_count),
                "unrooted_topology_count": str(summary.unrooted_topology_count),
                "clade_summary_count": str(summary.clade_summary_count),
                "unstable_clade_count": str(summary.unstable_clade_count),
                "top_unstable_clade": summary.top_unstable_clade or "",
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_discrete_ancestral_tree_set_node_table(
    path: Path,
    report: DiscreteAncestralTreeSetReport,
) -> Path:
    """Write one per-tree internal-node ledger for discrete ancestral tree-set analysis."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
            "clade_id",
            "clade_taxa",
            "most_likely_state",
            "state_set",
            "confidence",
            "ambiguous",
            "unstable",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
                "clade_id": row.clade_id,
                "clade_taxa": ",".join(row.clade_taxa),
                "most_likely_state": row.most_likely_state,
                "state_set": ",".join(row.state_set),
                "confidence": str(row.confidence),
                "ambiguous": str(row.ambiguous).lower(),
                "unstable": str(row.unstable).lower(),
            }
            for row in report.node_rows
        ],
    )


def write_discrete_ancestral_tree_set_clade_table(
    path: Path,
    report: DiscreteAncestralTreeSetReport,
) -> Path:
    """Write one comparable-clade summary ledger for discrete ancestral tree-set analysis."""
    return write_ancestral_rows(
        path,
        columns=[
            "clade_id",
            "clade_taxa",
            "tree_presence_count",
            "tree_presence_fraction",
            "unique_state_count",
            "dominant_state",
            "dominant_state_tree_count",
            "dominant_state_fraction",
            "ambiguous_tree_count",
            "ambiguous_tree_fraction",
            "unstable_tree_count",
            "unstable_tree_fraction",
            "state_distribution",
            "instability_score",
            "stability_class",
        ],
        rows=[
            {
                "clade_id": row.clade_id,
                "clade_taxa": ",".join(row.clade_taxa),
                "tree_presence_count": str(row.tree_presence_count),
                "tree_presence_fraction": str(row.tree_presence_fraction),
                "unique_state_count": str(row.unique_state_count),
                "dominant_state": row.dominant_state,
                "dominant_state_tree_count": str(row.dominant_state_tree_count),
                "dominant_state_fraction": str(row.dominant_state_fraction),
                "ambiguous_tree_count": str(row.ambiguous_tree_count),
                "ambiguous_tree_fraction": str(row.ambiguous_tree_fraction),
                "unstable_tree_count": str(row.unstable_tree_count),
                "unstable_tree_fraction": str(row.unstable_tree_fraction),
                "state_distribution": json.dumps(
                    row.state_distribution,
                    sort_keys=True,
                ),
                "instability_score": str(row.instability_score),
                "stability_class": row.stability_class,
            }
            for row in report.clade_summaries
        ],
    )
