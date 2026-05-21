from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.distance import (
    assess_imported_distance_method_assumptions,
    build_distance_method_report,
    build_tree_from_imported_distance_matrix,
    inspect_imported_distance_matrix_quality,
)
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.render.html import write_html_report

from .artifacts import section
from .models import DistanceReportBuildResult
from .summary import distance_method_limitations


def render_distance_report(
    *,
    out_path: Path,
    alignment_path: Path | None = None,
    matrix_path: Path | None = None,
) -> DistanceReportBuildResult:
    """Build a deterministic HTML report for computed or imported distance analysis."""
    if (alignment_path is None) == (matrix_path is None):
        raise ValueError(
            "render_distance_report requires exactly one of alignment_path or matrix_path"
        )

    method_limitations = distance_method_limitations()
    if alignment_path is not None:
        report = build_distance_method_report(alignment_path)
        title = "Bijux Distance Analysis Report"
        sections = [
            section("computed-distance-matrix", asdict(report.matrix)),
            section("distance-quality", asdict(report.quality)),
            section("distance-method-assumptions", asdict(report.assumptions)),
            section(
                "distance-reference-validation", asdict(report.reference_validation)
            ),
            section("neighbor-joining-tree", {"newick": report.built_tree_newick}),
            section("upgma-tree", {"newick": report.alternative_tree_newick}),
            section("distance-tree-comparison", asdict(report.topology_comparison)),
            section("distance-bootstrap-summary", asdict(report.bootstrap_summary)),
            section("distance-model-comparison", asdict(report.model_comparison)),
            section(
                "distance-gap-policy-sensitivity", asdict(report.gap_policy_sensitivity)
            ),
            section("distance-maturity-gate", asdict(report.maturity_gate)),
            section("distance-method-limitations", method_limitations),
        ]
        machine_manifest = {
            "report_kind": "distance-analysis",
            "source_kind": "alignment",
            "source_path": str(alignment_path),
            "method_limitations": method_limitations,
            "sections": [name for name, _ in sections],
        }
        write_html_report(
            title=title,
            sections=sections,
            out_path=out_path,
            embedded_json=machine_manifest,
        )
        return DistanceReportBuildResult(
            output_path=out_path,
            report_kind="distance-analysis",
            title=title,
            source_path=alignment_path,
            source_kind="alignment",
            method_limitations=method_limitations,
            machine_manifest=machine_manifest,
        )

    quality = inspect_imported_distance_matrix_quality(matrix_path)
    assumptions = assess_imported_distance_method_assumptions(matrix_path)
    title = "Bijux Imported Distance Report"
    sections = [
        section("imported-distance-matrix-quality", asdict(quality)),
        section("distance-method-assumptions", asdict(assumptions)),
        section("distance-method-limitations", method_limitations),
    ]
    validation = quality.validation
    if (
        validation.complete
        and validation.symmetric
        and validation.zero_diagonal
        and validation.nonnegative
    ):
        nj_tree, _ = build_tree_from_imported_distance_matrix(
            matrix_path, method="neighbor-joining"
        )
        upgma_tree, _ = build_tree_from_imported_distance_matrix(
            matrix_path, method="upgma"
        )
        sections.extend(
            [
                section("neighbor-joining-tree", {"newick": dumps_newick(nj_tree)}),
                section("upgma-tree", {"newick": dumps_newick(upgma_tree)}),
            ]
        )
    machine_manifest = {
        "report_kind": "distance-analysis",
        "source_kind": "imported-distance-matrix",
        "source_path": str(matrix_path),
        "method_limitations": method_limitations,
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return DistanceReportBuildResult(
        output_path=out_path,
        report_kind="distance-analysis",
        title=title,
        source_path=matrix_path,
        source_kind="imported-distance-matrix",
        method_limitations=method_limitations,
        machine_manifest=machine_manifest,
    )
