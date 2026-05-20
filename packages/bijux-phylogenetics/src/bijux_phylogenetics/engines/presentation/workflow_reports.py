from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from ..artifacts.bootstrap import build_low_support_bootstrap_rows
from ..common import load_engine_manifest
from ..artifacts.fasttree import build_fasttree_low_support_rows
from ..artifacts.sh_alrt import build_conflicting_sh_alrt_support_rows
from ..validation import (
    classify_inference_workflow_failure,
    compare_inferred_trees_across_engines,
    compare_ml_trees_across_models,
    detect_weakly_supported_backbone,
    summarize_bootstrap_support_distribution,
    summarize_fasttree_support_distribution,
    summarize_sh_alrt_support_distribution,
    validate_bootstrap_tree_set,
    validate_inference_engine_outputs,
    validate_ml_tree_contains_expected_taxa,
    validate_model_selection_against_engine_outputs,
)


@dataclass(slots=True)
class InferenceWorkflowReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    manifest_path: Path
    engine_name: str
    workflow: str
    warning_count: int
    supplement_sections: list[str]
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class ModelSelectionLimitationsReport:
    manifest_path: Path
    selected_model: str | None
    selected_criterion: str | None
    candidate_model_count: int
    best_model_aic: str | None
    best_model_aicc: str | None
    best_model_bic: str | None
    validation_issues: list[str]
    limitations: list[str]
    interpretation_limits: list[str]


@dataclass(slots=True)
class ModelSelectionLimitationsReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    manifest_path: Path
    selected_model: str | None
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class InferenceSensitivityReport:
    baseline_manifest_path: Path
    baseline_tree_path: Path
    baseline_engine_name: str
    baseline_selected_model: str | None
    alignment_filtering_sensitivity: object | None
    model_sensitivity: object | None
    engine_sensitivity: object | None
    bootstrap_support: object | None
    weak_backbone: object | None
    notes: list[str]


@dataclass(slots=True)
class InferenceSensitivityReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    baseline_manifest_path: Path
    machine_manifest: dict[str, object]


def render_inference_workflow_report(
    *, manifest_path: Path, out_path: Path
) -> InferenceWorkflowReportBuildResult:
    """Render a deterministic HTML report for one engine workflow manifest."""
    manifest = load_engine_manifest(manifest_path)
    output_paths = {
        key: Path(value) for key, value in dict(manifest["output_paths"]).items()
    }
    input_paths = [Path(path) for path in manifest["input_paths"]]
    failure = classify_inference_workflow_failure(
        workflow=str(manifest["workflow"]),
        input_paths=input_paths,
        output_paths=output_paths,
        run_exit_code=int(manifest["run"].get("exit_code", 0)),
    )
    consistency = validate_inference_engine_outputs(manifest_path)
    title = f"Bijux Inference Workflow Report: {manifest['workflow']}"
    sections = [
        (
            "workflow-summary",
            json.dumps(
                {
                    "engine_name": manifest["engine_name"],
                    "workflow": manifest["workflow"],
                    "resumed": manifest.get("resumed", False),
                    "selected_model": manifest.get("selected_model"),
                    "notes": manifest.get("notes", []),
                },
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "workflow-failure-taxonomy",
            json.dumps(asdict(failure), default=str, indent=2, sort_keys=True),
        ),
        (
            "workflow-consistency",
            json.dumps(asdict(consistency), default=str, indent=2, sort_keys=True),
        ),
        (
            "inputs",
            json.dumps(
                {
                    "input_paths": manifest["input_paths"],
                    "input_checksums": manifest.get("input_checksums", {}),
                },
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "outputs",
            json.dumps(
                {
                    "output_paths": manifest["output_paths"],
                    "output_checksums": manifest.get("output_checksums", {}),
                },
                indent=2,
                sort_keys=True,
            ),
        ),
        ("engine-run", json.dumps(manifest["run"], indent=2, sort_keys=True)),
    ]
    supplement_sections = [
        "workflow-summary",
        "workflow-failure-taxonomy",
        "workflow-consistency",
        "inputs",
        "outputs",
        "engine-run",
    ]
    if manifest["workflow"] == "model-selection":
        sections.append(
            (
                "model-selection-limitations",
                json.dumps(
                    asdict(build_model_selection_limitations_report(manifest_path)),
                    default=str,
                    indent=2,
                    sort_keys=True,
                ),
            )
        )
        supplement_sections.append("model-selection-limitations")
    if (
        manifest["workflow"] == "alignment-trimming"
        and manifest.get("trimming_summary") is not None
    ):
        sections.append(
            (
                "alignment-trimming-summary",
                json.dumps(
                    manifest["trimming_summary"],
                    indent=2,
                    sort_keys=True,
                ),
            )
        )
        supplement_sections.append("alignment-trimming-summary")
    if manifest["workflow"] == "maximum-likelihood-tree":
        sections.append(
            (
                "ml-tree-validation",
                json.dumps(
                    asdict(validate_ml_tree_contains_expected_taxa(manifest_path)),
                    default=str,
                    indent=2,
                    sort_keys=True,
                ),
            )
        )
        supplement_sections.append("ml-tree-validation")
    if manifest["workflow"] == "bootstrap-support":
        bootstrap_path = output_paths.get("bootstrap_trees")
        if bootstrap_path is not None:
            sections.append(
                (
                    "bootstrap-tree-set-validation",
                    json.dumps(
                        asdict(validate_bootstrap_tree_set(bootstrap_path)),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            supplement_sections.append("bootstrap-tree-set-validation")
        tree_path = output_paths.get("tree") or output_paths.get("support_tree")
        if tree_path is not None:
            bootstrap_support_summary = summarize_bootstrap_support_distribution(
                tree_path
            )
            sections.append(
                (
                    "bootstrap-support-summary",
                    json.dumps(
                        asdict(bootstrap_support_summary),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "bootstrap-support-histogram",
                    json.dumps(
                        bootstrap_support_summary.support_histogram,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "low-support-branches",
                    json.dumps(
                        [
                            asdict(row)
                            for row in build_low_support_bootstrap_rows(
                                bootstrap_support_summary
                            )
                        ],
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "weak-backbone",
                    json.dumps(
                        asdict(detect_weakly_supported_backbone(tree_path)),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            supplement_sections.extend(
                [
                    "bootstrap-support-summary",
                    "bootstrap-support-histogram",
                    "low-support-branches",
                    "weak-backbone",
                ]
            )
    if manifest["workflow"] == "fast-approximate-tree":
        tree_path = output_paths.get("tree") or output_paths.get("support_tree")
        if tree_path is not None:
            fasttree_support_summary = summarize_fasttree_support_distribution(
                tree_path
            )
            sections.append(
                (
                    "fasttree-approximation-limits",
                    json.dumps(
                        {
                            "approximate_method": fasttree_support_summary.approximate_method,
                            "support_label_kind": fasttree_support_summary.support_label_kind,
                            "support_scale": fasttree_support_summary.support_scale,
                            "limitations": [
                                "FastTree uses an approximately maximum-likelihood search and should be reviewed as exploratory or large-alignment evidence rather than as a direct substitute for a fully optimized ML workflow",
                                "FastTree local support values are SH-like proportions on a 0-to-1 scale and are not bootstrap percentages",
                            ],
                        },
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "fasttree-support-summary",
                    json.dumps(
                        asdict(fasttree_support_summary),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "fasttree-support-histogram",
                    json.dumps(
                        fasttree_support_summary.support_histogram,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "fasttree-low-support-branches",
                    json.dumps(
                        [
                            asdict(row)
                            for row in build_fasttree_low_support_rows(
                                fasttree_support_summary
                            )
                        ],
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            supplement_sections.extend(
                [
                    "fasttree-approximation-limits",
                    "fasttree-support-summary",
                    "fasttree-support-histogram",
                    "fasttree-low-support-branches",
                ]
            )
    if manifest["workflow"] == "sh-alrt-support":
        bootstrap_path = output_paths.get("bootstrap_trees")
        if bootstrap_path is not None:
            sections.append(
                (
                    "bootstrap-tree-set-validation",
                    json.dumps(
                        asdict(validate_bootstrap_tree_set(bootstrap_path)),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            supplement_sections.append("bootstrap-tree-set-validation")
        tree_path = output_paths.get("tree") or output_paths.get("support_tree")
        if tree_path is not None:
            sh_alrt_support_summary = summarize_sh_alrt_support_distribution(tree_path)
            sections.append(
                (
                    "sh-alrt-support-summary",
                    json.dumps(
                        asdict(sh_alrt_support_summary),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "conflicting-support-branches",
                    json.dumps(
                        [
                            asdict(row)
                            for row in build_conflicting_sh_alrt_support_rows(
                                sh_alrt_support_summary
                            )
                        ],
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "weak-backbone",
                    json.dumps(
                        asdict(detect_weakly_supported_backbone(tree_path)),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            supplement_sections.extend(
                [
                    "sh-alrt-support-summary",
                    "conflicting-support-branches",
                    "weak-backbone",
                ]
            )
    machine_manifest = {
        "report_kind": "inference-workflow",
        "title": title,
        "manifest_path": str(manifest_path),
        "engine_name": manifest["engine_name"],
        "workflow": manifest["workflow"],
        "warning_count": len(manifest["run"].get("warning_lines", [])),
        "sections": [name for name, _ in sections],
        "supplement_sections": supplement_sections,
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return InferenceWorkflowReportBuildResult(
        output_path=out_path,
        report_kind="inference-workflow",
        title=title,
        manifest_path=manifest_path,
        engine_name=str(manifest["engine_name"]),
        workflow=str(manifest["workflow"]),
        warning_count=len(manifest["run"].get("warning_lines", [])),
        supplement_sections=supplement_sections,
        machine_manifest=machine_manifest,
    )


def build_model_selection_limitations_report(
    manifest_path: Path,
) -> ModelSelectionLimitationsReport:
    """Describe why one selected substitution model should not be treated as biological truth."""
    manifest = load_engine_manifest(manifest_path)
    validation = validate_model_selection_against_engine_outputs(manifest_path)
    selected_model = manifest.get("selected_model")
    limitations = [
        "the selected model is the best fit among the engine's tested candidates, not proof of the true evolutionary process",
        "closely scoring models can yield similar likelihood support while implying different branch lengths or uncertainty",
        "model fit reflects the provided alignment, trimming policy, and missing-data pattern as much as the underlying biology",
        "an apparently optimal model can still be inappropriate if the alignment is compositionally biased, saturated, or poorly aligned",
    ]
    interpretation_limits = [
        "do not describe the selected model as the biologically correct mutation process without additional validation",
        "treat downstream topology or branch-length differences across plausible models as sensitivity, not noise",
        "report whether alternative near-best models were considered when scientific conclusions depend on branch support or branch lengths",
    ]
    return ModelSelectionLimitationsReport(
        manifest_path=manifest_path,
        selected_model=None if selected_model is None else str(selected_model),
        selected_criterion=validation.report_selected_criterion,
        candidate_model_count=validation.candidate_model_count,
        best_model_aic=validation.best_model_aic,
        best_model_aicc=validation.best_model_aicc,
        best_model_bic=validation.best_model_bic,
        validation_issues=validation.issues,
        limitations=limitations,
        interpretation_limits=interpretation_limits,
    )


def render_model_selection_limitations_report(
    *,
    manifest_path: Path,
    out_path: Path,
) -> ModelSelectionLimitationsReportBuildResult:
    """Render a reviewer-facing model-selection limitations report."""
    report = build_model_selection_limitations_report(manifest_path)
    title = "Bijux Model Selection Limitations Report"
    sections = [
        (
            "model-selection-summary",
            json.dumps(asdict(report), default=str, indent=2, sort_keys=True),
        ),
    ]
    machine_manifest = {
        "report_kind": "model-selection-limitations",
        "title": title,
        "manifest_path": str(manifest_path),
        "selected_model": report.selected_model,
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ModelSelectionLimitationsReportBuildResult(
        output_path=out_path,
        report_kind="model-selection-limitations",
        title=title,
        manifest_path=manifest_path,
        selected_model=report.selected_model,
        machine_manifest=machine_manifest,
    )


def build_inference_sensitivity_report(
    baseline_manifest_path: Path,
    *,
    filtered_manifest_path: Path | None = None,
    compare_model_manifest_path: Path | None = None,
    compare_engine_manifest_path: Path | None = None,
    bootstrap_manifest_path: Path | None = None,
) -> InferenceSensitivityReport:
    """Summarize how filtering, model choice, engine choice, and bootstrap support affect one inference workflow."""
    baseline_manifest = load_engine_manifest(baseline_manifest_path)
    baseline_output_paths = {
        key: Path(value)
        for key, value in dict(baseline_manifest["output_paths"]).items()
    }
    baseline_tree_path = baseline_output_paths["tree"]
    notes: list[str] = []
    alignment_filtering_sensitivity = None
    if filtered_manifest_path is not None:
        alignment_filtering_sensitivity = compare_ml_trees_across_models(
            baseline_manifest_path, filtered_manifest_path
        )
        notes.append(
            "alignment filtering sensitivity compares baseline and filtered-tree workflows"
        )
    model_sensitivity = None
    if compare_model_manifest_path is not None:
        model_sensitivity = compare_ml_trees_across_models(
            baseline_manifest_path, compare_model_manifest_path
        )
        notes.append("model sensitivity compares alternative selected-model workflows")
    engine_sensitivity = None
    if compare_engine_manifest_path is not None:
        engine_sensitivity = compare_inferred_trees_across_engines(
            baseline_manifest_path, compare_engine_manifest_path
        )
        notes.append(
            "engine sensitivity compares baseline inference against an alternative engine workflow"
        )
    bootstrap_support = None
    weak_backbone = None
    if bootstrap_manifest_path is not None:
        bootstrap_manifest = load_engine_manifest(bootstrap_manifest_path)
        bootstrap_output_paths = {
            key: Path(value)
            for key, value in dict(bootstrap_manifest["output_paths"]).items()
        }
        tree_path = bootstrap_output_paths.get("tree", baseline_tree_path)
        bootstrap_support = summarize_bootstrap_support_distribution(tree_path)
        weak_backbone = detect_weakly_supported_backbone(tree_path)
        notes.append(
            "bootstrap sensitivity summarizes support dispersion and backbone weakness on the supported tree"
        )
    return InferenceSensitivityReport(
        baseline_manifest_path=baseline_manifest_path,
        baseline_tree_path=baseline_tree_path,
        baseline_engine_name=str(baseline_manifest["engine_name"]),
        baseline_selected_model=None
        if baseline_manifest.get("selected_model") is None
        else str(baseline_manifest["selected_model"]),
        alignment_filtering_sensitivity=alignment_filtering_sensitivity,
        model_sensitivity=model_sensitivity,
        engine_sensitivity=engine_sensitivity,
        bootstrap_support=bootstrap_support,
        weak_backbone=weak_backbone,
        notes=notes,
    )


def render_inference_sensitivity_report(
    *,
    baseline_manifest_path: Path,
    out_path: Path,
    filtered_manifest_path: Path | None = None,
    compare_model_manifest_path: Path | None = None,
    compare_engine_manifest_path: Path | None = None,
    bootstrap_manifest_path: Path | None = None,
) -> InferenceSensitivityReportBuildResult:
    """Render a deterministic HTML sensitivity report for inference workflows."""
    report = build_inference_sensitivity_report(
        baseline_manifest_path,
        filtered_manifest_path=filtered_manifest_path,
        compare_model_manifest_path=compare_model_manifest_path,
        compare_engine_manifest_path=compare_engine_manifest_path,
        bootstrap_manifest_path=bootstrap_manifest_path,
    )
    title = "Bijux Inference Sensitivity Report"
    sections = [
        (
            "inference-sensitivity",
            json.dumps(asdict(report), default=str, indent=2, sort_keys=True),
        )
    ]
    machine_manifest = {
        "report_kind": "inference-sensitivity",
        "title": title,
        "baseline_manifest_path": str(baseline_manifest_path),
        "sections": [name for name, _ in sections],
        "has_alignment_filtering_sensitivity": report.alignment_filtering_sensitivity
        is not None,
        "has_model_sensitivity": report.model_sensitivity is not None,
        "has_engine_sensitivity": report.engine_sensitivity is not None,
        "has_bootstrap_support": report.bootstrap_support is not None,
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return InferenceSensitivityReportBuildResult(
        output_path=out_path,
        report_kind="inference-sensitivity",
        title=title,
        baseline_manifest_path=baseline_manifest_path,
        machine_manifest=machine_manifest,
    )
