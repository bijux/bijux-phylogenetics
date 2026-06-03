from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

from ....artifacts import section
from .shared import artifact_link, preview_payload


def build_topology_sections(
    *,
    clade_frequencies: Any,
    diversity: Any,
    clusters: Any,
    multimodality: Any | None,
    scaled_report_note: dict[str, object],
    artifact_paths: dict[str, Path],
    out_path: Path,
    preview_limit: int,
    clade_frequency_rows: list[dict[str, Any]],
    clade_frequency_truncated: int,
    rf_rows: list[dict[str, Any]],
    rf_truncated: int,
    cluster_rows: list[dict[str, Any]],
    cluster_truncated: int,
) -> list[tuple[str, object]]:
    """Build topology-focused sections for the tree uncertainty report."""
    return [
        section(
            "clade-frequencies",
            {
                "tree_count": clade_frequencies.tree_count,
                "shared_taxa": clade_frequencies.shared_taxa,
                **preview_payload(
                    rows=clade_frequency_rows,
                    row_count=len(clade_frequencies.clade_frequencies),
                    truncated_row_count=clade_frequency_truncated,
                    preview_limit=preview_limit,
                ),
                "artifact_path": artifact_link(
                    artifact_paths=artifact_paths,
                    artifact_key="clade_frequencies",
                    out_path=out_path,
                ),
            },
        ),
        section(
            "rf-distance-distribution",
            {
                "tree_count": diversity.tree_count,
                "pair_count": diversity.pair_count,
                **preview_payload(
                    rows=rf_rows,
                    row_count=len(diversity.rf_distribution),
                    truncated_row_count=rf_truncated,
                    preview_limit=preview_limit,
                ),
                "artifact_path": artifact_link(
                    artifact_paths=artifact_paths,
                    artifact_key="rf_distance_distribution",
                    out_path=out_path,
                ),
            },
        ),
        section(
            "topology-clusters",
            {
                "tree_count": clusters.tree_count,
                "rooted_topology_count": clusters.rooted_topology_count,
                **preview_payload(
                    rows=cluster_rows,
                    row_count=len(clusters.clusters),
                    truncated_row_count=cluster_truncated,
                    preview_limit=preview_limit,
                ),
                "artifact_path": artifact_link(
                    artifact_paths=artifact_paths,
                    artifact_key="topology_clusters",
                    out_path=out_path,
                ),
            },
        ),
        section(
            "topological-diversity",
            {
                **asdict(diversity),
                "artifact_path": artifact_link(
                    artifact_paths=artifact_paths,
                    artifact_key="topological_diversity",
                    out_path=out_path,
                ),
                "rf_distribution": (
                    f"{len(diversity.rf_distribution)} rows written to linked TSV"
                ),
            },
        ),
        section(
            "topology-multimodality",
            (
                {
                    **asdict(multimodality),
                    "artifact_path": artifact_link(
                        artifact_paths=artifact_paths,
                        artifact_key="topology_multimodality",
                        out_path=out_path,
                    ),
                }
                if multimodality is not None
                else {
                    **scaled_report_note,
                    "artifact_path": artifact_link(
                        artifact_paths=artifact_paths,
                        artifact_key="topology_multimodality",
                        out_path=out_path,
                    ),
                }
            ),
        ),
    ]
