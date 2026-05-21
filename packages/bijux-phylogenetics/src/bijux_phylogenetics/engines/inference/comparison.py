from __future__ import annotations

from dataclasses import asdict
from datetime import UTC, datetime
import json
from pathlib import Path
import shutil

from bijux_phylogenetics.compare.influence import analyze_taxon_influence
from bijux_phylogenetics.compare.presentation import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.compare.topology import write_tree_comparison_table
from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet
from bijux_phylogenetics.render.html import write_html_report

from ..common import build_file_checksums, write_engine_manifest
from ..validation import compare_inferred_trees_across_engines
from ..workflows.fasttree import run_fast_tree_inference
from ..workflows.iqtree import (
    run_bootstrap_support_estimation,
    run_model_selection,
)
from .tree_inference_comparison import (
    InferenceComparisonConclusionRow,
    InferenceComparisonConclusionSummary,
    InferenceComparisonConflictRow,
    InferenceComparisonSharedCladeRow,
    InferenceComparisonWeightedConflictRow,
    InferenceComparisonWorkflowReport,
    build_inference_comparison_conclusion_rows,
    build_inference_comparison_conflict_rows,
    build_inference_comparison_shared_clade_rows,
    build_inference_comparison_weighted_conflict_rows,
    summarize_inference_comparison_conclusions,
)

__all__ = [
    "InferenceComparisonConflictRow",
    "InferenceComparisonConclusionRow",
    "InferenceComparisonConclusionSummary",
    "InferenceComparisonSharedCladeRow",
    "InferenceComparisonWeightedConflictRow",
    "build_inference_comparison_conclusion_rows",
    "InferenceComparisonWorkflowReport",
    "build_inference_comparison_conflict_rows",
    "build_inference_comparison_shared_clade_rows",
    "build_inference_comparison_weighted_conflict_rows",
    "rewrite_inference_comparison_report_html",
    "run_tree_inference_comparison",
    "write_inference_comparison_clade_table",
    "write_inference_comparison_conclusion_table",
    "write_inference_comparison_summary_table",
    "write_inference_comparison_weighted_conflict_table",
]



def _artifact_prefix(engine_artifact_dir: Path, step: str) -> Path:
    return engine_artifact_dir / step / step


def _copy_output(source: Path, target: Path) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    return target


def _write_inference_comparison_taxon_influence_table(
    path: Path,
    report: TaxonInfluenceReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "influence_rank",
                "taxon",
                "retained_taxa",
                "baseline_support_disagreements",
                "leave_one_out_support_disagreements",
                "support_disagreement_delta",
                "baseline_conflicting_clades",
                "leave_one_out_conflicting_clades",
                "conflicting_clade_delta",
                "baseline_high_support_conflicts",
                "leave_one_out_high_support_conflicts",
                "high_support_conflict_delta",
                "topology_changed",
                "support_changed",
                "influence_score",
            ]
        )
    ]
    for row in report.rows:
        lines.append(
            "\t".join(
                [
                    str(row.influence_rank),
                    row.taxon,
                    "|".join(row.retained_taxa),
                    str(row.baseline_support_disagreements),
                    str(row.leave_one_out_support_disagreements),
                    str(row.support_disagreement_delta),
                    str(row.baseline_conflicting_clades),
                    str(row.leave_one_out_conflicting_clades),
                    str(row.conflicting_clade_delta),
                    str(row.baseline_high_support_conflicts),
                    str(row.leave_one_out_high_support_conflicts),
                    str(row.high_support_conflict_delta),
                    "true" if row.topology_changed else "false",
                    "true" if row.support_changed else "false",
                    _render_float(row.influence_score),
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_inference_comparison_weighted_conflict_table(
    path: Path,
    rows: list[InferenceComparisonWeightedConflictRow],
) -> Path:
    """Write one ranked support-weighted conflict ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "split_id",
                "comparison_status",
                "conflict_kind",
                "severity_class",
                "fasttree_support_fraction",
                "iqtree_support_fraction",
                "support_fraction_delta",
                "strongest_support_fraction",
                "support_weight",
                "serious_conflict",
                "detail",
            ]
        )
    ]
    for row in rows:
        lines.append(
            "\t".join(
                [
                    row.split_id,
                    row.comparison_status,
                    row.conflict_kind,
                    row.severity_class,
                    _render_float(row.fasttree_support_fraction),
                    _render_float(row.iqtree_support_fraction),
                    _render_float(row.support_fraction_delta),
                    _render_float(row.strongest_support_fraction),
                    _render_float(row.support_weight),
                    "true" if row.serious_conflict else "false",
                    row.detail,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_inference_comparison_conclusion_table(
    path: Path,
    rows: list[InferenceComparisonConclusionRow],
) -> Path:
    """Write one reviewer-facing clade stability ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "split_id",
                "conclusion_class",
                "evidence_class",
                "comparison_status",
                "fasttree_present",
                "iqtree_present",
                "fasttree_support_fraction",
                "iqtree_support_fraction",
                "support_fraction_delta",
                "serious_conflict",
                "detail",
            ]
        )
    ]
    for row in rows:
        lines.append(
            "\t".join(
                [
                    row.split_id,
                    row.conclusion_class,
                    row.evidence_class,
                    row.comparison_status,
                    "true" if row.fasttree_present else "false",
                    "true" if row.iqtree_present else "false",
                    _render_float(row.fasttree_support_fraction),
                    _render_float(row.iqtree_support_fraction),
                    _render_float(row.support_fraction_delta),
                    "true" if row.serious_conflict else "false",
                    row.detail,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_inference_comparison_summary_table(
    path: Path,
    summary: InferenceComparisonConclusionSummary,
) -> Path:
    """Write one compact summary row for the compared inference conclusions."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "shared_taxa_count",
                "robinson_foulds_distance",
                "normalized_robinson_foulds",
                "branch_score_distance",
                "stable_clade_count",
                "unstable_clade_count",
                "engine_specific_clade_count",
                "support_weighted_conflict_count",
                "low_support_disagreement_count",
                "moderate_support_disagreement_count",
                "high_support_conflict_count",
                "high_support_disagreement_count",
                "serious_conflict_count",
                "top_conflict_driver_taxa",
            ]
        ),
        "\t".join(
            [
                str(summary.shared_taxa_count),
                str(summary.robinson_foulds_distance),
                _render_float(summary.normalized_robinson_foulds),
                _render_float(summary.branch_score_distance),
                str(summary.stable_clade_count),
                str(summary.unstable_clade_count),
                str(summary.engine_specific_clade_count),
                str(summary.support_weighted_conflict_count),
                str(summary.low_support_disagreement_count),
                str(summary.moderate_support_disagreement_count),
                str(summary.high_support_conflict_count),
                str(summary.high_support_disagreement_count),
                str(summary.serious_conflict_count),
                "|".join(summary.top_conflict_driver_taxa),
            ]
        ),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def rewrite_inference_comparison_report_html(
    *,
    base_report: ComparisonReportBuildResult,
    summary: InferenceComparisonConclusionSummary,
    conclusion_rows: list[InferenceComparisonConclusionRow],
    weighted_conflict_rows: list[InferenceComparisonWeightedConflictRow],
    taxon_influence_report: TaxonInfluenceReport | None,
) -> Path:
    """Rewrite the comparison HTML with reviewer-facing stability synthesis."""
    sections = [
        (
            "comparison-summary",
            json.dumps(asdict(summary), indent=2, sort_keys=True),
        ),
        (
            "clade-conclusions",
            json.dumps(
                [asdict(row) for row in conclusion_rows],
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "support-weighted-conflicts",
            json.dumps(
                [asdict(row) for row in weighted_conflict_rows],
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "taxon-influence",
            json.dumps(
                None
                if taxon_influence_report is None
                else asdict(taxon_influence_report),
                indent=2,
                sort_keys=True,
                default=str,
            ),
        ),
        (
            "topology-metrics",
            json.dumps(
                asdict(base_report.topology), indent=2, sort_keys=True, default=str
            ),
        ),
        (
            "clade-comparison",
            json.dumps(
                asdict(base_report.clades), indent=2, sort_keys=True, default=str
            ),
        ),
        (
            "support-comparison",
            json.dumps(
                asdict(base_report.support), indent=2, sort_keys=True, default=str
            ),
        ),
        (
            "branch-length-comparison",
            json.dumps(
                asdict(base_report.branch_lengths),
                indent=2,
                sort_keys=True,
                default=str,
            ),
        ),
    ]
    return write_html_report(
        title="Bijux Tree Inference Comparison Report",
        sections=sections,
        out_path=base_report.output_path,
        embedded_json={
            "summary": asdict(summary),
            "conclusions": [asdict(row) for row in conclusion_rows],
            "weighted_conflicts": [asdict(row) for row in weighted_conflict_rows],
        },
    )


def write_inference_comparison_clade_table(
    path: Path,
    *,
    shared_rows: list[InferenceComparisonSharedCladeRow] | None = None,
    conflict_rows: list[InferenceComparisonConflictRow] | None = None,
) -> Path:
    """Write one shared-clade or conflicting-clade comparison table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str]
    if shared_rows is not None:
        lines = [
            "\t".join(
                [
                    "split_id",
                    "fasttree_support",
                    "fasttree_support_fraction",
                    "fasttree_support_label_kind",
                    "iqtree_support",
                    "iqtree_support_fraction",
                    "iqtree_support_label_kind",
                    "support_fraction_delta",
                    "support_disagreement",
                ]
            )
        ]
        for row in shared_rows:
            lines.append(
                "\t".join(
                    [
                        row.split_id,
                        _render_float(row.fasttree_support),
                        _render_float(row.fasttree_support_fraction),
                        row.fasttree_support_label_kind,
                        _render_float(row.iqtree_support),
                        _render_float(row.iqtree_support_fraction),
                        row.iqtree_support_label_kind,
                        _render_float(row.support_fraction_delta),
                        "true" if row.support_disagreement else "false",
                    ]
                )
            )
    elif conflict_rows is not None:
        lines = [
            "\t".join(
                [
                    "split_id",
                    "conflict_kind",
                    "fasttree_present",
                    "iqtree_present",
                    "fasttree_support",
                    "fasttree_support_fraction",
                    "iqtree_support",
                    "iqtree_support_fraction",
                    "detail",
                ]
            )
        ]
        for row in conflict_rows:
            lines.append(
                "\t".join(
                    [
                        row.split_id,
                        row.conflict_kind,
                        "true" if row.fasttree_present else "false",
                        "true" if row.iqtree_present else "false",
                        _render_float(row.fasttree_support),
                        _render_float(row.fasttree_support_fraction),
                        _render_float(row.iqtree_support),
                        _render_float(row.iqtree_support_fraction),
                        row.detail,
                    ]
                )
            )
    else:
        raise ValueError("one row set must be provided when writing comparison tables")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _render_float(value: float | None) -> str:
    return "" if value is None else format(value, ".12g")


def run_tree_inference_comparison(
    input_path: Path,
    *,
    out_dir: Path,
    prefix: str = "engine-comparison",
    sequence_type: AlignmentAlphabet | None = None,
    iqtree_executable: str | Path = "iqtree2",
    fasttree_executable: str | Path = "FastTree",
    iqtree_seed: int = 1,
    iqtree_threads: int = 1,
    bootstrap_replicates: int = 1000,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> InferenceComparisonWorkflowReport:
    """Run IQ-TREE and FastTree on one alignment and compare the inferred trees."""
    started_at = datetime.now(UTC)
    workflow_prefix = prefix
    out_dir.mkdir(parents=True, exist_ok=True)
    engine_artifact_dir = out_dir / "engine-artifacts" / workflow_prefix
    engine_artifact_dir.mkdir(parents=True, exist_ok=True)

    final_outputs = {
        "fasttree_tree": out_dir / f"{workflow_prefix}.fasttree.nwk",
        "iqtree_support_tree": out_dir / f"{workflow_prefix}.iqtree-support.nwk",
        "comparison_report": out_dir / f"{workflow_prefix}.comparison.html",
        "stability_summary": out_dir / f"{workflow_prefix}.stability-summary.tsv",
        "conclusion_table": out_dir / f"{workflow_prefix}.conclusions.tsv",
        "support_weighted_conflicts": out_dir
        / f"{workflow_prefix}.support-weighted-conflicts.tsv",
        "taxon_influence": out_dir / f"{workflow_prefix}.taxon-influence.tsv",
        "comparison_table": out_dir / f"{workflow_prefix}.comparison.tsv",
        "shared_clades": out_dir / f"{workflow_prefix}.shared-clades.tsv",
        "conflicting_clades": out_dir / f"{workflow_prefix}.conflicting-clades.tsv",
        "manifest": out_dir / f"{workflow_prefix}.manifest.json",
    }

    model_selection_workflow = run_model_selection(
        input_path,
        out_dir=_artifact_prefix(engine_artifact_dir, "model-selection").parent,
        executable=iqtree_executable,
        prefix="model-selection",
        sequence_type=sequence_type,
        resume=resume,
        seed=iqtree_seed,
        threads=iqtree_threads,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )
    if model_selection_workflow.selected_model is None:
        raise ValueError("model-selection workflow did not expose a selected model")
    iqtree_support_workflow = run_bootstrap_support_estimation(
        input_path,
        out_dir=_artifact_prefix(engine_artifact_dir, "bootstrap-support").parent,
        model=model_selection_workflow.selected_model,
        replicates=bootstrap_replicates,
        prefix="bootstrap-support",
        executable=iqtree_executable,
        sequence_type=sequence_type,
        resume=resume,
        seed=iqtree_seed,
        threads=iqtree_threads,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )
    fasttree_workflow = run_fast_tree_inference(
        input_path,
        _artifact_prefix(engine_artifact_dir, "fasttree").with_suffix(".nwk"),
        executable=fasttree_executable,
        sequence_type=sequence_type,
        resume=resume,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )

    _copy_output(fasttree_workflow.output_paths["tree"], final_outputs["fasttree_tree"])
    _copy_output(
        iqtree_support_workflow.output_paths["support_tree"],
        final_outputs["iqtree_support_tree"],
    )
    comparison_report = build_tree_comparison_report(
        final_outputs["fasttree_tree"],
        final_outputs["iqtree_support_tree"],
        out_path=final_outputs["comparison_report"],
    )
    write_tree_comparison_table(
        final_outputs["comparison_table"],
        final_outputs["fasttree_tree"],
        final_outputs["iqtree_support_tree"],
    )
    engine_comparison = compare_inferred_trees_across_engines(
        fasttree_workflow.manifest_path,
        iqtree_support_workflow.manifest_path,
    )
    shared_clade_rows = build_inference_comparison_shared_clade_rows(engine_comparison)
    conflicting_clade_rows = build_inference_comparison_conflict_rows(
        engine_comparison,
        comparison_report,
    )
    weighted_conflict_rows = build_inference_comparison_weighted_conflict_rows(
        engine_comparison.support
    )
    conclusion_rows = build_inference_comparison_conclusion_rows(
        engine_comparison.support
    )
    taxon_influence_report: TaxonInfluenceReport | None = None
    try:
        taxon_influence_report = analyze_taxon_influence(
            final_outputs["fasttree_tree"],
            final_outputs["iqtree_support_tree"],
        )
    except ValueError:
        taxon_influence_report = None
    conclusion_summary = summarize_inference_comparison_conclusions(
        engine_comparison.topology,
        engine_comparison.branch_lengths,
        weighted_conflict_rows=weighted_conflict_rows,
        conclusion_rows=conclusion_rows,
        taxon_influence_report=taxon_influence_report,
    )
    write_inference_comparison_clade_table(
        final_outputs["shared_clades"],
        shared_rows=shared_clade_rows,
    )
    write_inference_comparison_clade_table(
        final_outputs["conflicting_clades"],
        conflict_rows=conflicting_clade_rows,
    )
    write_inference_comparison_summary_table(
        final_outputs["stability_summary"],
        conclusion_summary,
    )
    write_inference_comparison_conclusion_table(
        final_outputs["conclusion_table"],
        conclusion_rows,
    )
    write_inference_comparison_weighted_conflict_table(
        final_outputs["support_weighted_conflicts"],
        weighted_conflict_rows,
    )
    if taxon_influence_report is not None:
        _write_inference_comparison_taxon_influence_table(
            final_outputs["taxon_influence"],
            taxon_influence_report,
        )
    else:
        final_outputs.pop("taxon_influence")
    rewrite_inference_comparison_report_html(
        base_report=comparison_report,
        summary=conclusion_summary,
        conclusion_rows=conclusion_rows,
        weighted_conflict_rows=weighted_conflict_rows,
        taxon_influence_report=taxon_influence_report,
    )

    warnings = list(
        dict.fromkeys(
            model_selection_workflow.run.warning_lines
            + iqtree_support_workflow.run.warning_lines
            + fasttree_workflow.run.warning_lines
            + engine_comparison.warnings
            + (
                []
                if fasttree_workflow.fasttree_support_summary is None
                else fasttree_workflow.fasttree_support_summary.warnings
            )
        )
    )
    notes = [
        "comparison workflow runs IQ-TREE model selection, IQ-TREE bootstrap-supported inference, and FastTree approximate inference on the same alignment",
        "comparison table exports one flat split ledger across topology, support, and branch-length comparisons",
        "shared-clade and conflicting-clade ledgers separate topological agreement from topological or support disagreement",
        "support-weighted conflict and conclusion ledgers separate stable clades, unstable shared clades, and engine-specific clades",
        "FastTree SH-like local support and IQ-TREE ultrafast bootstrap support are normalized to fractions only for side-by-side review and not as proof that the methods are interchangeable",
        *(
            []
            if not conclusion_summary.top_conflict_driver_taxa
            else [
                "leave-one-out taxon influence identified the strongest conflict-driving taxa"
            ]
        ),
        f"iqtree random seed: {iqtree_seed}",
        f"iqtree threads: {iqtree_threads}",
        f"ultrafast bootstrap replicates: {bootstrap_replicates}",
    ]
    report = InferenceComparisonWorkflowReport(
        workflow="tree-inference-comparison",
        input_path=input_path,
        out_dir=out_dir,
        prefix=workflow_prefix,
        sequence_type=sequence_type,
        selected_model=model_selection_workflow.selected_model,
        iqtree_seed=iqtree_seed,
        iqtree_threads=iqtree_threads,
        bootstrap_replicates=bootstrap_replicates,
        timeout_seconds=timeout_seconds,
        started_at_utc=started_at.replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z"),
        ended_at_utc="",
        runtime_seconds=0.0,
        engine_artifact_dir=engine_artifact_dir,
        manifest_path=final_outputs["manifest"],
        output_paths=final_outputs,
        step_manifests={
            "model_selection": model_selection_workflow.manifest_path,
            "iqtree_support": iqtree_support_workflow.manifest_path,
            "fasttree": fasttree_workflow.manifest_path,
        },
        config={
            "sequence_type": sequence_type,
            "iqtree_seed": iqtree_seed,
            "iqtree_threads": iqtree_threads,
            "bootstrap_replicates": bootstrap_replicates,
            "timeout_seconds": timeout_seconds,
            "resume": resume,
            "incomplete_run_policy": incomplete_run_policy,
        },
        commands={
            "model_selection": model_selection_workflow.run.command,
            "iqtree_support": iqtree_support_workflow.run.command,
            "fasttree": fasttree_workflow.run.command,
        },
        engine_versions={
            "iqtree_model_selection": model_selection_workflow.run.version.text,
            "iqtree_support": iqtree_support_workflow.run.version.text,
            "fasttree": fasttree_workflow.run.version.text,
        },
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        model_selection_workflow=model_selection_workflow,
        iqtree_support_workflow=iqtree_support_workflow,
        fasttree_workflow=fasttree_workflow,
        engine_comparison=engine_comparison,
        comparison_report=comparison_report,
        shared_clade_rows=shared_clade_rows,
        conflicting_clade_rows=conflicting_clade_rows,
        weighted_conflict_rows=weighted_conflict_rows,
        conclusion_rows=conclusion_rows,
        conclusion_summary=conclusion_summary,
        taxon_influence_report=taxon_influence_report,
        warnings=warnings,
        notes=notes,
    )
    ended_at = datetime.now(UTC)
    report.ended_at_utc = (
        ended_at.replace(microsecond=0).isoformat().replace("+00:00", "Z")
    )
    report.runtime_seconds = max(
        0.0,
        round((ended_at - started_at).total_seconds(), 6),
    )
    report.output_checksums = build_file_checksums(
        [
            final_outputs["fasttree_tree"],
            final_outputs["iqtree_support_tree"],
            final_outputs["comparison_report"],
            final_outputs["stability_summary"],
            final_outputs["conclusion_table"],
            final_outputs["support_weighted_conflicts"],
            final_outputs["comparison_table"],
            final_outputs["shared_clades"],
            final_outputs["conflicting_clades"],
            *(
                []
                if taxon_influence_report is None
                else [final_outputs["taxon_influence"]]
            ),
        ]
    )
    write_engine_manifest(report.manifest_path, report)
    return report
