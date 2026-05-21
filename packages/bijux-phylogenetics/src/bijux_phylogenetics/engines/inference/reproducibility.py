from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from bijux_phylogenetics.compare.topology import (
    compare_support_values,
    compare_tree_paths,
)
from bijux_phylogenetics.io.iqtree_support import support_fraction
from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet

from ..common import build_file_checksums, write_engine_manifest
from ..workflows.iqtree import (
    run_bootstrap_support_estimation,
    run_model_selection,
)
from ..workflows.models import EngineWorkflowReport

__all__ = [
    "InferenceReproducibilityComparisonRow",
    "InferenceReproducibilityRunRow",
    "InferenceReproducibilitySupportDeltaRow",
    "InferenceReproducibilityWorkflowReport",
    "run_inference_reproducibility_check",
    "write_inference_reproducibility_table",
]

_LOG_LIKELIHOOD_TOLERANCE = 1e-6
_SUPPORT_FRACTION_TOLERANCE = 1e-9


@dataclass(frozen=True, slots=True)
class InferenceReproducibilityRunRow:
    """One repeated supported-inference run tracked by the reproducibility workflow."""

    run_index: int
    run_label: str
    run_role: str
    manifest_path: Path
    support_tree_path: Path
    log_likelihood: float | None
    support_value_count: int
    minimum_support: float | None
    maximum_support: float | None


@dataclass(frozen=True, slots=True)
class InferenceReproducibilitySupportDeltaRow:
    """One clade-level support delta between the baseline and one rerun."""

    rerun_index: int
    rerun_label: str
    split_id: str
    baseline_support: float | None
    baseline_support_fraction: float | None
    rerun_support: float | None
    rerun_support_fraction: float | None
    support_fraction_delta: float | None
    exceeds_tolerance: bool


@dataclass(frozen=True, slots=True)
class InferenceReproducibilityComparisonRow:
    """One baseline-versus-rerun reproducibility comparison record."""

    baseline_run_index: int
    baseline_run_label: str
    rerun_index: int
    rerun_label: str
    topology_equal: bool
    same_unrooted_topology: bool
    same_topology_different_branch_lengths: bool
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    shared_taxa_count: int
    shared_clade_count: int
    support_difference_count: int
    maximum_support_fraction_delta: float | None
    log_likelihood_delta: float | None
    classification: str
    detail: str


@dataclass(slots=True)
class InferenceReproducibilityWorkflowReport:
    """End-to-end reproducibility result for repeated supported IQ-TREE inference."""

    workflow: str
    input_path: Path
    out_dir: Path
    prefix: str
    sequence_type: AlignmentAlphabet | None
    selected_model: str
    repeat_count: int
    bootstrap_replicates: int
    iqtree_seed: int
    iqtree_threads: int
    started_at_utc: str
    ended_at_utc: str
    runtime_seconds: float
    engine_artifact_dir: Path
    manifest_path: Path
    output_paths: dict[str, Path]
    step_manifests: dict[str, Path]
    config: dict[str, object]
    commands: dict[str, list[str]]
    engine_versions: dict[str, str]
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    model_selection_workflow: EngineWorkflowReport
    run_workflows: list[EngineWorkflowReport]
    run_rows: list[InferenceReproducibilityRunRow]
    comparison_rows: list[InferenceReproducibilityComparisonRow]
    support_delta_rows: list[InferenceReproducibilitySupportDeltaRow]
    overall_status: str
    warnings: list[str]
    notes: list[str]


def _artifact_prefix(engine_artifact_dir: Path, run_label: str) -> Path:
    return engine_artifact_dir / run_label / run_label


def _run_label(run_index: int) -> str:
    return "baseline" if run_index == 1 else f"rerun-{run_index}"


def _run_role(run_index: int) -> str:
    return "baseline" if run_index == 1 else "rerun"


def _render_float(value: float | None) -> str:
    return "" if value is None else format(value, ".12g")


def _support_delta(
    left_support: float | None,
    right_support: float | None,
) -> tuple[float | None, bool]:
    left_fraction = support_fraction(left_support)
    right_fraction = support_fraction(right_support)
    if left_fraction is None and right_fraction is None:
        return None, False
    if left_fraction is None or right_fraction is None:
        return None, True
    delta = abs(left_fraction - right_fraction)
    return delta, delta > _SUPPORT_FRACTION_TOLERANCE


def _build_run_row(
    run_index: int,
    workflow: EngineWorkflowReport,
) -> InferenceReproducibilityRunRow:
    if workflow.iqtree_summary is None:
        raise ValueError("bootstrap-support workflow did not expose an IQ-TREE summary")
    return InferenceReproducibilityRunRow(
        run_index=run_index,
        run_label=_run_label(run_index),
        run_role=_run_role(run_index),
        manifest_path=workflow.manifest_path,
        support_tree_path=workflow.output_paths["support_tree"],
        log_likelihood=workflow.log_likelihood,
        support_value_count=workflow.iqtree_summary.support_value_count,
        minimum_support=workflow.iqtree_summary.minimum_support,
        maximum_support=workflow.iqtree_summary.maximum_support,
    )


def _build_support_delta_rows(
    baseline_row: InferenceReproducibilityRunRow,
    rerun_row: InferenceReproducibilityRunRow,
) -> list[InferenceReproducibilitySupportDeltaRow]:
    comparison = compare_support_values(
        baseline_row.support_tree_path,
        rerun_row.support_tree_path,
    )
    rows: list[InferenceReproducibilitySupportDeltaRow] = []
    for pair in comparison.shared_clades:
        delta, exceeds_tolerance = _support_delta(pair.left_support, pair.right_support)
        rows.append(
            InferenceReproducibilitySupportDeltaRow(
                rerun_index=rerun_row.run_index,
                rerun_label=rerun_row.run_label,
                split_id=pair.split_id,
                baseline_support=pair.left_support,
                baseline_support_fraction=support_fraction(pair.left_support),
                rerun_support=pair.right_support,
                rerun_support_fraction=support_fraction(pair.right_support),
                support_fraction_delta=delta,
                exceeds_tolerance=exceeds_tolerance,
            )
        )
    return rows


def _build_comparison_row(
    baseline_row: InferenceReproducibilityRunRow,
    rerun_row: InferenceReproducibilityRunRow,
    support_delta_rows: list[InferenceReproducibilitySupportDeltaRow],
) -> InferenceReproducibilityComparisonRow:
    topology = compare_tree_paths(
        baseline_row.support_tree_path,
        rerun_row.support_tree_path,
    )
    support_difference_count = sum(
        1 for row in support_delta_rows if row.exceeds_tolerance
    )
    comparable_deltas = [
        row.support_fraction_delta
        for row in support_delta_rows
        if row.support_fraction_delta is not None
    ]
    maximum_support_fraction_delta = (
        None if not comparable_deltas else max(comparable_deltas)
    )
    log_likelihood_delta = (
        None
        if baseline_row.log_likelihood is None or rerun_row.log_likelihood is None
        else abs(baseline_row.log_likelihood - rerun_row.log_likelihood)
    )
    likelihood_matches = (
        log_likelihood_delta is None
        or log_likelihood_delta <= _LOG_LIKELIHOOD_TOLERANCE
    )
    support_matches = support_difference_count == 0
    if topology.topology_equal and likelihood_matches and support_matches:
        classification = "deterministic"
        detail = "topology, likelihood, and support values match the baseline"
    elif topology.same_unrooted_topology and likelihood_matches and support_matches:
        classification = "equivalent"
        detail = "unrooted topology and support values match but rooting or branch lengths differ"
    else:
        detail_parts: list[str] = []
        if not topology.same_unrooted_topology:
            detail_parts.append("unrooted topology differs")
        elif not topology.topology_equal:
            detail_parts.append("rooting or rooted clade interpretation differs")
        if not likelihood_matches:
            detail_parts.append("log-likelihood differs")
        if support_difference_count > 0:
            detail_parts.append(
                f"{support_difference_count} shared clades differ in support beyond tolerance"
            )
        classification = "unstable"
        detail = "; ".join(detail_parts)
    return InferenceReproducibilityComparisonRow(
        baseline_run_index=baseline_row.run_index,
        baseline_run_label=baseline_row.run_label,
        rerun_index=rerun_row.run_index,
        rerun_label=rerun_row.run_label,
        topology_equal=topology.topology_equal,
        same_unrooted_topology=topology.same_unrooted_topology,
        same_topology_different_branch_lengths=topology.same_topology_different_branch_lengths,
        robinson_foulds_distance=topology.robinson_foulds_distance,
        normalized_robinson_foulds=topology.normalized_robinson_foulds,
        shared_taxa_count=len(topology.shared_taxa),
        shared_clade_count=len(support_delta_rows),
        support_difference_count=support_difference_count,
        maximum_support_fraction_delta=maximum_support_fraction_delta,
        log_likelihood_delta=log_likelihood_delta,
        classification=classification,
        detail=detail,
    )


def _overall_status(
    comparison_rows: list[InferenceReproducibilityComparisonRow],
) -> str:
    if all(row.classification == "deterministic" for row in comparison_rows):
        return "deterministic"
    if all(
        row.classification in {"deterministic", "equivalent"} for row in comparison_rows
    ):
        return "equivalent"
    return "unstable"


def write_inference_reproducibility_table(
    path: Path,
    *,
    run_rows: list[InferenceReproducibilityRunRow] | None = None,
    comparison_rows: list[InferenceReproducibilityComparisonRow] | None = None,
    support_delta_rows: list[InferenceReproducibilitySupportDeltaRow] | None = None,
) -> Path:
    """Write one reproducibility run, comparison, or support-delta TSV ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if run_rows is not None:
        lines = [
            "\t".join(
                [
                    "run_index",
                    "run_label",
                    "run_role",
                    "manifest_path",
                    "support_tree_path",
                    "log_likelihood",
                    "support_value_count",
                    "minimum_support",
                    "maximum_support",
                ]
            )
        ]
        for row in run_rows:
            lines.append(
                "\t".join(
                    [
                        str(row.run_index),
                        row.run_label,
                        row.run_role,
                        str(row.manifest_path),
                        str(row.support_tree_path),
                        _render_float(row.log_likelihood),
                        str(row.support_value_count),
                        _render_float(row.minimum_support),
                        _render_float(row.maximum_support),
                    ]
                )
            )
    elif comparison_rows is not None:
        lines = [
            "\t".join(
                [
                    "baseline_run_index",
                    "baseline_run_label",
                    "rerun_index",
                    "rerun_label",
                    "topology_equal",
                    "same_unrooted_topology",
                    "same_topology_different_branch_lengths",
                    "robinson_foulds_distance",
                    "normalized_robinson_foulds",
                    "shared_taxa_count",
                    "shared_clade_count",
                    "support_difference_count",
                    "maximum_support_fraction_delta",
                    "log_likelihood_delta",
                    "classification",
                    "detail",
                ]
            )
        ]
        for row in comparison_rows:
            lines.append(
                "\t".join(
                    [
                        str(row.baseline_run_index),
                        row.baseline_run_label,
                        str(row.rerun_index),
                        row.rerun_label,
                        "true" if row.topology_equal else "false",
                        "true" if row.same_unrooted_topology else "false",
                        "true"
                        if row.same_topology_different_branch_lengths
                        else "false",
                        str(row.robinson_foulds_distance),
                        _render_float(row.normalized_robinson_foulds),
                        str(row.shared_taxa_count),
                        str(row.shared_clade_count),
                        str(row.support_difference_count),
                        _render_float(row.maximum_support_fraction_delta),
                        _render_float(row.log_likelihood_delta),
                        row.classification,
                        row.detail,
                    ]
                )
            )
    elif support_delta_rows is not None:
        lines = [
            "\t".join(
                [
                    "rerun_index",
                    "rerun_label",
                    "split_id",
                    "baseline_support",
                    "baseline_support_fraction",
                    "rerun_support",
                    "rerun_support_fraction",
                    "support_fraction_delta",
                    "exceeds_tolerance",
                ]
            )
        ]
        for row in support_delta_rows:
            lines.append(
                "\t".join(
                    [
                        str(row.rerun_index),
                        row.rerun_label,
                        row.split_id,
                        _render_float(row.baseline_support),
                        _render_float(row.baseline_support_fraction),
                        _render_float(row.rerun_support),
                        _render_float(row.rerun_support_fraction),
                        _render_float(row.support_fraction_delta),
                        "true" if row.exceeds_tolerance else "false",
                    ]
                )
            )
    else:
        raise ValueError(
            "one row set must be provided when writing reproducibility tables"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run_inference_reproducibility_check(
    input_path: Path,
    *,
    out_dir: Path,
    prefix: str = "inference-reproducibility",
    sequence_type: AlignmentAlphabet | None = None,
    executable: str | Path = "iqtree2",
    repeats: int = 3,
    bootstrap_replicates: int = 1000,
    seed: int = 1,
    threads: int = 1,
) -> InferenceReproducibilityWorkflowReport:
    """Rerun supported IQ-TREE inference and classify deterministic versus unstable output."""
    started_at = datetime.now(UTC)
    if repeats < 2:
        raise ValueError(f"repeats must be at least 2, got {repeats}")
    workflow_prefix = prefix
    out_dir.mkdir(parents=True, exist_ok=True)
    engine_artifact_dir = out_dir / "engine-artifacts" / workflow_prefix
    engine_artifact_dir.mkdir(parents=True, exist_ok=True)

    final_outputs = {
        "runs_table": out_dir / f"{workflow_prefix}.runs.tsv",
        "comparison_table": out_dir / f"{workflow_prefix}.comparisons.tsv",
        "support_delta_table": out_dir / f"{workflow_prefix}.support-deltas.tsv",
        "manifest": out_dir / f"{workflow_prefix}.manifest.json",
    }

    model_selection_workflow = run_model_selection(
        input_path,
        out_dir=_artifact_prefix(engine_artifact_dir, "model-selection").parent,
        executable=executable,
        prefix="model-selection",
        sequence_type=sequence_type,
        seed=seed,
        threads=threads,
    )
    if model_selection_workflow.selected_model is None:
        raise ValueError("model-selection workflow did not expose a selected model")

    run_workflows: list[EngineWorkflowReport] = []
    run_rows: list[InferenceReproducibilityRunRow] = []
    step_manifests = {"model_selection": model_selection_workflow.manifest_path}
    for run_index in range(1, repeats + 1):
        run_label = _run_label(run_index)
        workflow = run_bootstrap_support_estimation(
            input_path,
            out_dir=_artifact_prefix(engine_artifact_dir, run_label).parent,
            model=model_selection_workflow.selected_model,
            replicates=bootstrap_replicates,
            prefix=run_label,
            executable=executable,
            sequence_type=sequence_type,
            seed=seed,
            threads=threads,
        )
        run_workflows.append(workflow)
        run_rows.append(_build_run_row(run_index, workflow))
        step_manifests[run_label] = workflow.manifest_path

    baseline_row = run_rows[0]
    comparison_rows: list[InferenceReproducibilityComparisonRow] = []
    support_delta_rows: list[InferenceReproducibilitySupportDeltaRow] = []
    for rerun_row in run_rows[1:]:
        rerun_support_delta_rows = _build_support_delta_rows(baseline_row, rerun_row)
        support_delta_rows.extend(rerun_support_delta_rows)
        comparison_rows.append(
            _build_comparison_row(
                baseline_row,
                rerun_row,
                rerun_support_delta_rows,
            )
        )

    write_inference_reproducibility_table(
        final_outputs["runs_table"],
        run_rows=run_rows,
    )
    write_inference_reproducibility_table(
        final_outputs["comparison_table"],
        comparison_rows=comparison_rows,
    )
    write_inference_reproducibility_table(
        final_outputs["support_delta_table"],
        support_delta_rows=support_delta_rows,
    )

    overall_status = _overall_status(comparison_rows)
    warnings = list(
        dict.fromkeys(
            [
                *model_selection_workflow.run.warning_lines,
                *(
                    warning
                    for workflow in run_workflows
                    for warning in workflow.run.warning_lines
                ),
                *(
                    []
                    if overall_status != "unstable"
                    else [
                        "one or more reruns changed topology, likelihood, or support beyond the governed tolerances"
                    ]
                ),
            ]
        )
    )
    notes = [
        "the reproducibility workflow runs model selection once to choose a fixed model, then reruns identical supported IQ-TREE inference settings against that chosen model",
        "comparison rows classify each rerun as deterministic, equivalent, or unstable relative to the baseline run",
        "support deltas are compared as normalized fractions so IQ-TREE percentage labels can be reviewed under a stable 0..1 tolerance contract",
        f"iqtree random seed: {seed}",
        f"iqtree threads: {threads}",
        f"ultrafast bootstrap replicates: {bootstrap_replicates}",
        f"repeat count: {repeats}",
    ]
    report = InferenceReproducibilityWorkflowReport(
        workflow="inference-reproducibility",
        input_path=input_path,
        out_dir=out_dir,
        prefix=workflow_prefix,
        sequence_type=sequence_type,
        selected_model=model_selection_workflow.selected_model,
        repeat_count=repeats,
        bootstrap_replicates=bootstrap_replicates,
        iqtree_seed=seed,
        iqtree_threads=threads,
        started_at_utc=started_at.replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z"),
        ended_at_utc="",
        runtime_seconds=0.0,
        engine_artifact_dir=engine_artifact_dir,
        manifest_path=final_outputs["manifest"],
        output_paths=final_outputs,
        step_manifests=step_manifests,
        config={
            "sequence_type": sequence_type,
            "repeats": repeats,
            "bootstrap_replicates": bootstrap_replicates,
            "seed": seed,
            "threads": threads,
        },
        commands={
            "model_selection": model_selection_workflow.run.command,
            **{
                row.run_label: workflow.run.command
                for row, workflow in zip(run_rows, run_workflows, strict=True)
            },
        },
        engine_versions={
            "iqtree_model_selection": model_selection_workflow.run.version.text,
            **{
                row.run_label: workflow.run.version.text
                for row, workflow in zip(run_rows, run_workflows, strict=True)
            },
        },
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        model_selection_workflow=model_selection_workflow,
        run_workflows=run_workflows,
        run_rows=run_rows,
        comparison_rows=comparison_rows,
        support_delta_rows=support_delta_rows,
        overall_status=overall_status,
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
            final_outputs["runs_table"],
            final_outputs["comparison_table"],
            final_outputs["support_delta_table"],
        ]
    )
    write_engine_manifest(report.manifest_path, report)
    return report
