from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from ..common import load_engine_manifest
from ..validation import (
    compare_inferred_trees_across_engines,
    compare_ml_trees_across_models,
    detect_weakly_supported_backbone,
    summarize_bootstrap_support_distribution,
)


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
