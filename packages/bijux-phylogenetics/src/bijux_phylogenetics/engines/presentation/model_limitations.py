from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from ..common import load_engine_manifest
from ..validation import validate_model_selection_against_engine_outputs


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
