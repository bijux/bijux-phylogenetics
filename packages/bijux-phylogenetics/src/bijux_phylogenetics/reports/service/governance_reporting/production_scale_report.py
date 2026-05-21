from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.validation import build_production_scale_readiness_report

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..models import ProductionScaleReadinessReportBuildResult


def render_production_scale_readiness_report(
    *,
    out_path: Path,
    replicates: int = 1,
    tree_tip_counts: list[int] | None = None,
    alignment_size_classes: list[tuple[str, int, int]] | None = None,
    tree_set_size_classes: list[tuple[str, int, int]] | None = None,
    stress_tiers: list[str] | None = None,
) -> ProductionScaleReadinessReportBuildResult:
    """Render one reviewer-facing production-scale readiness report from governed benchmark evidence."""
    production_scale_readiness = build_production_scale_readiness_report(
        replicates=replicates,
        tree_tip_counts=tree_tip_counts,
        alignment_size_classes=alignment_size_classes,
        tree_set_size_classes=tree_set_size_classes,
        stress_tiers=stress_tiers,
    )
    title = "Bijux Production-Scale Readiness Report"
    highest_ready_scale_counts = {
        scale: sum(
            1
            for entry in production_scale_readiness.entries
            if entry.highest_ready_scale == scale
        )
        for scale in sorted(
            {
                "below-small",
                *(
                    threshold.scale
                    for threshold in production_scale_readiness.scale_definitions
                ),
            }
        )
    }
    scale_coverage = [
        {
            "scale": threshold.scale,
            "description": threshold.description,
            "minimum_taxa": threshold.minimum_taxa,
            "minimum_sites": threshold.minimum_sites,
            "minimum_tree_count": threshold.minimum_tree_count,
            "minimum_posterior_size": threshold.minimum_posterior_size,
            "ready_workflow_count": sum(
                1
                for entry in production_scale_readiness.entries
                for decision in entry.scale_decisions
                if decision.scale == threshold.scale and decision.ready
            ),
            "ready_workflows": sorted(
                entry.workflow
                for entry in production_scale_readiness.entries
                for decision in entry.scale_decisions
                if decision.scale == threshold.scale and decision.ready
            ),
        }
        for threshold in production_scale_readiness.scale_definitions
    ]
    reviewer_summary = [
        f"workflow count: {len(production_scale_readiness.entries)}",
        "highest ready scale distribution: "
        + ", ".join(
            f"{scale}={count}"
            for scale, count in highest_ready_scale_counts.items()
            if count > 0
        ),
        f"stress tiers: {', '.join(production_scale_readiness.stress_tiers)}",
    ]
    sections = [
        section("reviewer-summary", reviewer_summary),
        section(
            "scale-definitions",
            [asdict(item) for item in production_scale_readiness.scale_definitions],
        ),
        section("scale-coverage", scale_coverage),
        section(
            "production-scale-readiness",
            [asdict(item) for item in production_scale_readiness.entries],
        ),
        section("known-limitations", production_scale_readiness.limitations),
    ]
    machine_manifest = {
        "report_kind": "production-scale-readiness",
        "title": title,
        "input_paths": [],
        "input_checksums": {},
        "sections": [name for name, _ in sections],
        "metrics": {
            "goal_id": production_scale_readiness.goal_id,
            "workflow_count": len(production_scale_readiness.entries),
            "replicates": production_scale_readiness.replicates,
            "stress_tier_count": len(production_scale_readiness.stress_tiers),
            "scale_definition_count": len(production_scale_readiness.scale_definitions),
            "below_small_workflow_count": highest_ready_scale_counts.get(
                "below-small", 0
            ),
            **{
                f"{threshold.scale}_ready_workflow_count": sum(
                    1
                    for entry in production_scale_readiness.entries
                    for decision in entry.scale_decisions
                    if decision.scale == threshold.scale and decision.ready
                )
                for threshold in production_scale_readiness.scale_definitions
            },
        },
        "reviewer_summary": reviewer_summary,
        "limitations": production_scale_readiness.limitations,
    }
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path),
        machine_manifest,
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ProductionScaleReadinessReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="production-scale-readiness",
        title=title,
        production_scale_readiness=production_scale_readiness,
        machine_manifest=machine_manifest,
    )
