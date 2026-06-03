from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

from ....artifacts import preview_report_rows, section
from .shared import artifact_link


def build_sensitivity_sections(
    *,
    storage_risk: Any,
    thinning_sensitivity: Any | None,
    consensus_sensitivity: Any | None,
    benchmark: Any | None,
    benchmark_tree_count: int | None,
    benchmark_taxon_count: int | None,
    maturity: Any | None,
    summary: Any,
    scaled_report_note: dict[str, object],
    artifact_paths: dict[str, Path],
    out_path: Path,
    preview_limit: int,
    thinning_rows: list[dict[str, Any]],
    thinning_truncated: int,
    consensus_rows: list[dict[str, Any]],
    consensus_truncated: int,
    benchmark_rows: list[dict[str, Any]],
    benchmark_truncated: int,
) -> list[tuple[str, object]]:
    """Build sensitivity and readiness sections for the tree uncertainty report."""
    return [
        section(
            "storage-risk",
            {
                **asdict(storage_risk),
                "artifact_path": artifact_link(
                    artifact_paths=artifact_paths,
                    artifact_key="storage_risk",
                    out_path=out_path,
                ),
            },
        ),
        section(
            "thinning-sensitivity",
            {
                **(
                    {
                        "path": str(thinning_sensitivity.path),
                        "original_tree_count": thinning_sensitivity.original_tree_count,
                        "original_rooted_topology_count": thinning_sensitivity.original_rooted_topology_count,
                        "original_dominant_topology_frequency": thinning_sensitivity.original_dominant_topology_frequency,
                        "warning_count": len(thinning_sensitivity.warnings),
                        "warnings": thinning_sensitivity.warnings,
                        "row_count": len(thinning_sensitivity.rows),
                        "truncated_row_count": thinning_truncated,
                        "preview_row_count": min(len(thinning_rows), preview_limit),
                        "preview_rows": preview_report_rows(
                            thinning_rows, limit=preview_limit
                        ),
                    }
                    if thinning_sensitivity is not None
                    else scaled_report_note
                ),
                "artifact_path": artifact_link(
                    artifact_paths=artifact_paths,
                    artifact_key="thinning_sensitivity",
                    out_path=out_path,
                ),
            },
        ),
        section(
            "consensus-threshold-sensitivity",
            {
                **(
                    {
                        "path": str(consensus_sensitivity.path),
                        "tree_count": consensus_sensitivity.tree_count,
                        "warning_count": len(consensus_sensitivity.warnings),
                        "warnings": consensus_sensitivity.warnings,
                        "row_count": len(consensus_sensitivity.rows),
                        "truncated_row_count": consensus_truncated,
                        "preview_row_count": min(len(consensus_rows), preview_limit),
                        "preview_rows": preview_report_rows(
                            consensus_rows, limit=preview_limit
                        ),
                    }
                    if consensus_sensitivity is not None
                    else scaled_report_note
                ),
                "artifact_path": artifact_link(
                    artifact_paths=artifact_paths,
                    artifact_key="consensus_threshold_sensitivity",
                    out_path=out_path,
                ),
            },
        ),
        section(
            "tree-set-benchmark",
            {
                **(
                    {
                        "tree_counts": benchmark.tree_counts,
                        "taxon_counts": benchmark.taxon_counts,
                        "sampled_tree_count": benchmark_tree_count,
                        "sampled_taxon_count": benchmark_taxon_count,
                        "benchmark_capped": benchmark_tree_count != summary.tree_count
                        or benchmark_taxon_count != max(len(summary.shared_taxa), 2),
                        "row_count": len(benchmark.rows),
                        "truncated_row_count": benchmark_truncated,
                        "preview_row_count": min(len(benchmark_rows), preview_limit),
                        "preview_rows": preview_report_rows(
                            benchmark_rows, limit=preview_limit
                        ),
                    }
                    if benchmark is not None
                    else scaled_report_note
                ),
                "artifact_path": artifact_link(
                    artifact_paths=artifact_paths,
                    artifact_key="tree_set_benchmark",
                    out_path=out_path,
                ),
            },
        ),
        section(
            "maturity-gate",
            (
                {
                    **asdict(maturity),
                    "artifact_path": artifact_link(
                        artifact_paths=artifact_paths,
                        artifact_key="maturity_gate",
                        out_path=out_path,
                    ),
                }
                if maturity is not None
                else {
                    **scaled_report_note,
                    "artifact_path": artifact_link(
                        artifact_paths=artifact_paths,
                        artifact_key="maturity_gate",
                        out_path=out_path,
                    ),
                }
            ),
        ),
    ]
