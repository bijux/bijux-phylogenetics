from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .contracts import GeographicSamplingBiasReport


def write_geographic_sampling_bias_summary_table(
    path: Path,
    report: GeographicSamplingBiasReport,
) -> Path:
    """Write one summary ledger for geographic sampling-bias review."""
    summary = report.summary
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "internal_model",
            "weighting_mode",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observed_region_count",
            "region_dominated",
            "dominant_region",
            "dominant_region_fraction",
            "weighted_region_dominated",
            "weighted_dominant_region",
            "weighted_dominant_region_fraction",
            "root_region_unweighted",
            "root_region_weighted",
            "root_region_changed",
            "compared_internal_node_count",
            "changed_internal_node_count",
            "compared_transition_count",
            "changed_transition_count",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "internal_model": summary.internal_model,
                "weighting_mode": summary.weighting_mode,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "observed_region_count": str(summary.observed_region_count),
                "region_dominated": str(summary.region_dominated).lower(),
                "dominant_region": summary.dominant_region,
                "dominant_region_fraction": str(summary.dominant_region_fraction),
                "weighted_region_dominated": str(
                    summary.weighted_region_dominated
                ).lower(),
                "weighted_dominant_region": summary.weighted_dominant_region,
                "weighted_dominant_region_fraction": str(
                    summary.weighted_dominant_region_fraction
                ),
                "root_region_unweighted": summary.root_region_unweighted,
                "root_region_weighted": summary.root_region_weighted,
                "root_region_changed": str(summary.root_region_changed).lower(),
                "compared_internal_node_count": str(
                    summary.compared_internal_node_count
                ),
                "changed_internal_node_count": str(summary.changed_internal_node_count),
                "compared_transition_count": str(summary.compared_transition_count),
                "changed_transition_count": str(summary.changed_transition_count),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_geographic_sampling_count_table(
    path: Path,
    report: GeographicSamplingBiasReport,
) -> Path:
    """Write one region sample-count and weighting ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "region",
            "sample_count",
            "sample_fraction",
            "applied_weight",
            "weighted_sample_count",
            "weighted_sample_fraction",
            "dominant_unweighted",
            "dominant_weighted",
        ],
        rows=[
            {
                "region": row.region,
                "sample_count": str(row.sample_count),
                "sample_fraction": str(row.sample_fraction),
                "applied_weight": str(row.applied_weight),
                "weighted_sample_count": str(row.weighted_sample_count),
                "weighted_sample_fraction": str(row.weighted_sample_fraction),
                "dominant_unweighted": str(row.dominant_unweighted).lower(),
                "dominant_weighted": str(row.dominant_weighted).lower(),
            }
            for row in report.count_rows
        ],
    )


def write_geographic_sampling_bias_node_table(
    path: Path,
    report: GeographicSamplingBiasReport,
) -> Path:
    """Write one weighted-versus-unweighted node comparison ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "is_root",
            "unweighted_region",
            "weighted_region",
            "unweighted_confidence",
            "weighted_confidence",
            "confidence_delta",
            "changed",
            "unweighted_region_probabilities",
            "weighted_region_probabilities",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "descendant_taxa": ",".join(row.descendant_taxa),
                "is_root": str(row.is_root).lower(),
                "unweighted_region": row.unweighted_region,
                "weighted_region": row.weighted_region,
                "unweighted_confidence": str(row.unweighted_confidence),
                "weighted_confidence": str(row.weighted_confidence),
                "confidence_delta": str(row.confidence_delta),
                "changed": str(row.changed).lower(),
                "unweighted_region_probabilities": json.dumps(
                    row.unweighted_region_probabilities,
                    sort_keys=True,
                ),
                "weighted_region_probabilities": json.dumps(
                    row.weighted_region_probabilities,
                    sort_keys=True,
                ),
            }
            for row in report.node_rows
        ],
    )


def write_geographic_sampling_bias_transition_table(
    path: Path,
    report: GeographicSamplingBiasReport,
) -> Path:
    """Write one weighted-versus-unweighted transition comparison ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "unweighted_source_region",
            "unweighted_target_region",
            "weighted_source_region",
            "weighted_target_region",
            "unweighted_transition",
            "weighted_transition",
            "unweighted_changed",
            "weighted_changed",
            "changed_by_weighting",
            "unweighted_support",
            "weighted_support",
        ],
        rows=[
            {
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "unweighted_source_region": row.unweighted_source_region,
                "unweighted_target_region": row.unweighted_target_region,
                "weighted_source_region": row.weighted_source_region,
                "weighted_target_region": row.weighted_target_region,
                "unweighted_transition": row.unweighted_transition,
                "weighted_transition": row.weighted_transition,
                "unweighted_changed": str(row.unweighted_changed).lower(),
                "weighted_changed": str(row.weighted_changed).lower(),
                "changed_by_weighting": str(row.changed_by_weighting).lower(),
                "unweighted_support": str(row.unweighted_support),
                "weighted_support": str(row.weighted_support),
            }
            for row in report.transition_rows
        ],
    )


def write_geographic_sampling_bias_exclusion_table(
    path: Path,
    report: GeographicSamplingBiasReport,
) -> Path:
    """Write one excluded-taxa ledger for geographic sampling-bias review."""
    return write_taxon_rows(
        path,
        columns=["taxon", "raw_region", "normalized_region", "reason", "note"],
        rows=[
            {
                "taxon": row.taxon,
                "raw_region": row.raw_region,
                "normalized_region": row.normalized_region or "",
                "reason": row.reason,
                "note": row.note,
            }
            for row in report.exclusion_rows
        ],
    )
