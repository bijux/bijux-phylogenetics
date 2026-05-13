from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil

from bijux_phylogenetics.compare.reports import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.compare.topology import write_tree_comparison_table
from bijux_phylogenetics.core.alignment import AlignmentAlphabet
from bijux_phylogenetics.io.iqtree_support import support_fraction

from .common import build_file_checksums, write_engine_manifest
from .validation import InferenceTreeComparisonReport, compare_inferred_trees_across_engines
from .workflows import (
    EngineWorkflowReport,
    run_bootstrap_support_estimation,
    run_fast_tree_inference,
    run_model_selection,
)

__all__ = [
    "InferenceComparisonConflictRow",
    "InferenceComparisonSharedCladeRow",
    "InferenceComparisonWorkflowReport",
    "build_inference_comparison_conflict_rows",
    "build_inference_comparison_shared_clade_rows",
    "run_tree_inference_comparison",
    "write_inference_comparison_clade_table",
]

_SUPPORT_DISAGREEMENT_THRESHOLD = 0.15


@dataclass(frozen=True, slots=True)
class InferenceComparisonSharedCladeRow:
    """One shared clade across FastTree and IQ-TREE support workflows."""

    split_id: str
    fasttree_support: float | None
    fasttree_support_fraction: float | None
    fasttree_support_label_kind: str
    iqtree_support: float | None
    iqtree_support_fraction: float | None
    iqtree_support_label_kind: str
    support_fraction_delta: float | None
    support_disagreement: bool


@dataclass(frozen=True, slots=True)
class InferenceComparisonConflictRow:
    """One clade-level conflict record across the two inference engines."""

    split_id: str
    conflict_kind: str
    fasttree_present: bool
    iqtree_present: bool
    fasttree_support: float | None
    fasttree_support_fraction: float | None
    iqtree_support: float | None
    iqtree_support_fraction: float | None
    detail: str


@dataclass(slots=True)
class InferenceComparisonWorkflowReport:
    """End-to-end result for one engine-comparison workflow on one alignment."""

    input_path: Path
    out_dir: Path
    prefix: str
    sequence_type: AlignmentAlphabet | None
    selected_model: str
    engine_artifact_dir: Path
    manifest_path: Path
    output_paths: dict[str, Path]
    step_manifests: dict[str, Path]
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    model_selection_workflow: EngineWorkflowReport
    iqtree_support_workflow: EngineWorkflowReport
    fasttree_workflow: EngineWorkflowReport
    engine_comparison: InferenceTreeComparisonReport
    comparison_report: ComparisonReportBuildResult
    shared_clade_rows: list[InferenceComparisonSharedCladeRow]
    conflicting_clade_rows: list[InferenceComparisonConflictRow]
    warnings: list[str]
    notes: list[str]


def _artifact_prefix(engine_artifact_dir: Path, step: str) -> Path:
    return engine_artifact_dir / step / step


def _copy_output(source: Path, target: Path) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    return target


def _normalize_fasttree_support(value: float | None) -> float | None:
    return value


def _normalize_iqtree_support(value: float | None) -> float | None:
    return support_fraction(value)


def build_inference_comparison_shared_clade_rows(
    comparison: InferenceTreeComparisonReport,
    *,
    support_disagreement_threshold: float = _SUPPORT_DISAGREEMENT_THRESHOLD,
) -> list[InferenceComparisonSharedCladeRow]:
    """Convert one compared engine pair into shared-clade rows."""
    rows: list[InferenceComparisonSharedCladeRow] = []
    for pair in comparison.support.shared_clades:
        fasttree_support_fraction = _normalize_fasttree_support(pair.left_support)
        iqtree_support_fraction = _normalize_iqtree_support(pair.right_support)
        support_fraction_delta = (
            None
            if fasttree_support_fraction is None or iqtree_support_fraction is None
            else abs(fasttree_support_fraction - iqtree_support_fraction)
        )
        rows.append(
            InferenceComparisonSharedCladeRow(
                split_id=pair.split_id,
                fasttree_support=pair.left_support,
                fasttree_support_fraction=fasttree_support_fraction,
                fasttree_support_label_kind="sh-like-local-support",
                iqtree_support=pair.right_support,
                iqtree_support_fraction=iqtree_support_fraction,
                iqtree_support_label_kind="ufboot-support",
                support_fraction_delta=support_fraction_delta,
                support_disagreement=(
                    support_fraction_delta is not None
                    and support_fraction_delta >= support_disagreement_threshold
                ),
            )
        )
    return rows


def build_inference_comparison_conflict_rows(
    comparison: InferenceTreeComparisonReport,
    comparison_report: ComparisonReportBuildResult,
    *,
    support_disagreement_threshold: float = _SUPPORT_DISAGREEMENT_THRESHOLD,
) -> list[InferenceComparisonConflictRow]:
    """Build one combined topology-plus-support conflict ledger."""
    del comparison_report
    rows: list[InferenceComparisonConflictRow] = []
    shared_rows = build_inference_comparison_shared_clade_rows(
        comparison,
        support_disagreement_threshold=support_disagreement_threshold,
    )
    for support_conflict in comparison.support.conflicting_clades:
        if support_conflict.comparison_status == "left_only":
            conflict_kind = "fasttree_only"
            fasttree_present = True
            iqtree_present = False
        elif support_conflict.comparison_status == "right_only":
            conflict_kind = "iqtree_only"
            fasttree_present = False
            iqtree_present = True
        else:
            continue
        rows.append(
            InferenceComparisonConflictRow(
                split_id=support_conflict.split_id,
                conflict_kind=conflict_kind,
                fasttree_present=fasttree_present,
                iqtree_present=iqtree_present,
                fasttree_support=support_conflict.left_support,
                fasttree_support_fraction=support_conflict.left_support_fraction,
                iqtree_support=support_conflict.right_support,
                iqtree_support_fraction=support_conflict.right_support_fraction,
                detail=support_conflict.detail,
            )
        )
    for row in shared_rows:
        if not row.support_disagreement:
            continue
        rows.append(
            InferenceComparisonConflictRow(
                split_id=row.split_id,
                conflict_kind="support_disagreement",
                fasttree_present=True,
                iqtree_present=True,
                fasttree_support=row.fasttree_support,
                fasttree_support_fraction=row.fasttree_support_fraction,
                iqtree_support=row.iqtree_support,
                iqtree_support_fraction=row.iqtree_support_fraction,
                detail=(
                    "normalized support fractions differ by at least "
                    f"{support_disagreement_threshold:.2f}"
                ),
            )
        )
    return rows


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
    workflow_prefix = prefix
    out_dir.mkdir(parents=True, exist_ok=True)
    engine_artifact_dir = out_dir / "engine-artifacts" / workflow_prefix
    engine_artifact_dir.mkdir(parents=True, exist_ok=True)

    final_outputs = {
        "fasttree_tree": out_dir / f"{workflow_prefix}.fasttree.nwk",
        "iqtree_support_tree": out_dir / f"{workflow_prefix}.iqtree-support.nwk",
        "comparison_report": out_dir / f"{workflow_prefix}.comparison.html",
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
    write_inference_comparison_clade_table(
        final_outputs["shared_clades"],
        shared_rows=shared_clade_rows,
    )
    write_inference_comparison_clade_table(
        final_outputs["conflicting_clades"],
        conflict_rows=conflicting_clade_rows,
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
        "FastTree SH-like local support and IQ-TREE ultrafast bootstrap support are normalized to fractions only for side-by-side review and not as proof that the methods are interchangeable",
        f"iqtree random seed: {iqtree_seed}",
        f"iqtree threads: {iqtree_threads}",
        f"ultrafast bootstrap replicates: {bootstrap_replicates}",
    ]
    report = InferenceComparisonWorkflowReport(
        input_path=input_path,
        out_dir=out_dir,
        prefix=workflow_prefix,
        sequence_type=sequence_type,
        selected_model=model_selection_workflow.selected_model,
        engine_artifact_dir=engine_artifact_dir,
        manifest_path=final_outputs["manifest"],
        output_paths=final_outputs,
        step_manifests={
            "model_selection": model_selection_workflow.manifest_path,
            "iqtree_support": iqtree_support_workflow.manifest_path,
            "fasttree": fasttree_workflow.manifest_path,
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
        warnings=warnings,
        notes=notes,
    )
    report.output_checksums = build_file_checksums(
        [
            final_outputs["fasttree_tree"],
            final_outputs["iqtree_support_tree"],
            final_outputs["comparison_report"],
            final_outputs["comparison_table"],
            final_outputs["shared_clades"],
            final_outputs["conflicting_clades"],
        ]
    )
    write_engine_manifest(report.manifest_path, report)
    return report
