from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path


@dataclass(frozen=True, slots=True)
class ReviewerAuditChecklistItem:
    """Reviewer-facing pass/risk/block item for one analysis package surface."""

    section: str
    status: str
    summary: str
    evidence: list[str]
    artifact_paths: list[str]


@dataclass(slots=True)
class ReviewerAuditChecklist:
    """Reviewer-readable checklist over one machine-produced analysis package."""

    report_kind: str
    items: list[ReviewerAuditChecklistItem]


@dataclass(slots=True)
class ReviewerAuditChecklistWriteResult:
    """Result of writing one reviewer audit checklist artifact."""

    output_path: Path
    checklist: ReviewerAuditChecklist


def _mapping(payload: dict[str, object], key: str) -> dict[str, object]:
    value = payload.get(key)
    if isinstance(value, dict):
        return value
    return {}


def _text_list(payload: object) -> list[str]:
    if not isinstance(payload, list):
        return []
    return [str(item) for item in payload if str(item).strip()]


def _artifact_paths(mapping: dict[str, object], *keys: str) -> list[str]:
    paths: list[str] = []
    for key in keys:
        value = mapping.get(key)
        if isinstance(value, str) and value.strip():
            paths.append(value)
    return paths


def _status(*, blocked: bool = False, risk: bool = False) -> str:
    if blocked:
        return "blocked"
    if risk:
        return "risk"
    return "pass"


def _tree_checklist(manifest: dict[str, object]) -> ReviewerAuditChecklist:
    outputs = _mapping(manifest, "outputs")
    metrics = _mapping(manifest, "metrics")
    validation = _mapping(manifest, "validation")
    inspection = _mapping(manifest, "inspection")
    support_audit = _mapping(manifest, "support_audit")
    limitations = _text_list(manifest.get("limitations"))
    validation_warnings = _text_list(validation.get("warnings"))
    inspection_warnings = _text_list(inspection.get("warnings"))
    validity_decision = str(validation.get("validity_decision", "unknown"))
    rendered_support_count = int(metrics.get("rendered_support_count", 0))
    supported_branch_count = int(metrics.get("supported_branch_count", 0))
    support_validated = bool(support_audit.get("validated"))
    items = [
        ReviewerAuditChecklistItem(
            section="inputs",
            status="pass",
            summary="input tree path and checksum were recorded",
            evidence=[
                f"input path: {manifest.get('input_path', '')}",
                f"input checksum: {manifest.get('input_checksum', '')}",
            ],
            artifact_paths=[],
        ),
        ReviewerAuditChecklistItem(
            section="methods",
            status=_status(
                blocked=not bool(outputs.get("methods_summary_path")),
                risk=not bool(str(manifest.get("methods_summary_text", "")).strip()),
            ),
            summary="tree validation methods text is available for reproducibility review",
            evidence=[
                f"methods summary path: {outputs.get('methods_summary_path', '')}",
            ],
            artifact_paths=_artifact_paths(outputs, "methods_summary_path"),
        ),
        ReviewerAuditChecklistItem(
            section="validity",
            status=_status(
                blocked=validity_decision != "valid",
                risk=bool(validation_warnings or inspection_warnings),
            ),
            summary=f"tree validity decision is {validity_decision}",
            evidence=(validation_warnings + inspection_warnings)[:8],
            artifact_paths=_artifact_paths(
                outputs,
                "report_path",
                "clade_table_path",
                "branch_stats_path",
            ),
        ),
        ReviewerAuditChecklistItem(
            section="support_surface",
            status=_status(
                blocked=not support_validated,
                risk=rendered_support_count < supported_branch_count,
            ),
            summary="support evidence was checked against the rendered review figure",
            evidence=[
                f"support audit validated: {str(support_validated).lower()}",
                f"supported branches: {supported_branch_count}",
                f"rendered support labels: {rendered_support_count}",
            ]
            + _text_list(support_audit.get("warnings")),
            artifact_paths=_artifact_paths(
                outputs,
                "figure_path",
                "support_table_path",
            ),
        ),
        ReviewerAuditChecklistItem(
            section="interpretation_limits",
            status=_status(risk=bool(limitations)),
            summary="interpretation limits were recorded explicitly for reviewer inspection",
            evidence=limitations[:8],
            artifact_paths=_artifact_paths(outputs, "report_path"),
        ),
    ]
    return ReviewerAuditChecklist(report_kind="tree_package", items=items)


def _comparative_checklist(manifest: dict[str, object]) -> ReviewerAuditChecklist:
    outputs = _mapping(manifest, "outputs")
    metrics = _mapping(manifest, "metrics")
    summary = _mapping(manifest, "summary")
    limitations = _text_list(manifest.get("limitations"))
    delta = float(summary.get("better_model_aicc_delta", 0.0))
    items = [
        ReviewerAuditChecklistItem(
            section="inputs",
            status="pass",
            summary="tree and trait input paths plus checksums were recorded",
            evidence=[
                f"input count: {len(_text_list(manifest.get('input_paths')))}",
                f"checksum count: {len(_mapping(manifest, 'input_checksums'))}",
            ],
            artifact_paths=[],
        ),
        ReviewerAuditChecklistItem(
            section="methods",
            status=_status(blocked=not bool(outputs.get("methods_summary_path"))),
            summary="comparative methods text is available for reproducibility review",
            evidence=[
                f"methods summary path: {outputs.get('methods_summary_path', '')}",
                f"methods summary warnings: {metrics.get('methods_summary_warning_count', 0)}",
            ],
            artifact_paths=_artifact_paths(outputs, "methods_summary_path"),
        ),
        ReviewerAuditChecklistItem(
            section="model_selection",
            status=_status(risk=delta < 2.0),
            summary="comparative model support was quantified with an explicit AICc gap",
            evidence=[
                f"selected model: {summary.get('selected_model', '')}",
                f"runner-up AICc delta: {delta:.6f}",
            ],
            artifact_paths=_artifact_paths(
                outputs,
                "model_comparison_table_path",
                "summary_table_path",
            ),
        ),
        ReviewerAuditChecklistItem(
            section="diagnostics",
            status=_status(
                risk=int(metrics.get("limitation_count", 0)) > 0
                or int(metrics.get("methods_summary_warning_count", 0)) > 0
            ),
            summary="comparative residual, signal, and exclusion diagnostics were retained",
            evidence=[
                f"limitation count: {metrics.get('limitation_count', 0)}",
                f"coefficient count: {metrics.get('coefficient_count', 0)}",
                f"contrast count: {metrics.get('contrast_count', 0)}",
            ]
            + limitations[:5],
            artifact_paths=_artifact_paths(
                outputs,
                "residual_table_path",
                "signal_table_path",
                "audit_table_path",
                "contrast_table_path",
            ),
        ),
        ReviewerAuditChecklistItem(
            section="interpretation_limits",
            status=_status(risk=bool(limitations)),
            summary="comparative limitations remain explicit instead of being implied by coefficients alone",
            evidence=limitations[:8],
            artifact_paths=_artifact_paths(
                outputs,
                "interpretation_table_path",
                "report_path",
            ),
        ),
    ]
    return ReviewerAuditChecklist(report_kind="comparative_package", items=items)


def _ancestral_checklist(manifest: dict[str, object]) -> ReviewerAuditChecklist:
    outputs = _mapping(manifest, "outputs")
    metrics = _mapping(manifest, "metrics")
    machine_report_manifest = _mapping(manifest, "machine_report_manifest")
    limitations = _text_list(machine_report_manifest.get("limitations"))
    warning_count = int(metrics.get("warning_count", 0))
    items = [
        ReviewerAuditChecklistItem(
            section="inputs",
            status="pass",
            summary="tree and trait input paths plus checksums were recorded",
            evidence=[
                f"input count: {len(_text_list(manifest.get('input_paths')))}",
                f"checksum count: {len(_mapping(manifest, 'input_checksums'))}",
            ],
            artifact_paths=[],
        ),
        ReviewerAuditChecklistItem(
            section="methods",
            status=_status(blocked=not bool(outputs.get("methods_summary_path"))),
            summary="ancestral reconstruction methods text is available for reproducibility review",
            evidence=[
                f"methods summary path: {outputs.get('methods_summary_path', '')}",
                f"methods summary warnings: {metrics.get('methods_summary_warning_count', 0)}",
            ],
            artifact_paths=_artifact_paths(outputs, "methods_summary_path"),
        ),
        ReviewerAuditChecklistItem(
            section="uncertainty_surface",
            status=_status(
                blocked=not bool(outputs.get("uncertainty_table_path"))
                or not bool(outputs.get("figure_path"))
            ),
            summary="node uncertainty remained visible in both figure and ledger form",
            evidence=[
                f"uncertainty table path: {outputs.get('uncertainty_table_path', '')}",
                f"figure path: {outputs.get('figure_path', '')}",
            ],
            artifact_paths=_artifact_paths(
                outputs,
                "figure_path",
                "uncertainty_table_path",
                "node_table_path",
            ),
        ),
        ReviewerAuditChecklistItem(
            section="transition_review",
            status=_status(
                blocked=not bool(outputs.get("transition_count_table_path"))
                or not bool(outputs.get("transition_branch_table_path"))
            ),
            summary="branchwise transition or change-review ledgers were retained",
            evidence=[
                f"transition count rows: {metrics.get('transition_count_row_count', 0)}",
                f"transition branch rows: {metrics.get('transition_branch_row_count', 0)}",
            ],
            artifact_paths=_artifact_paths(
                outputs,
                "transition_count_table_path",
                "transition_branch_table_path",
                "exclusion_table_path",
            ),
        ),
        ReviewerAuditChecklistItem(
            section="interpretation_limits",
            status=_status(risk=warning_count > 0 or bool(limitations)),
            summary="ancestral interpretation limits and exclusions were recorded explicitly",
            evidence=[
                f"warning count: {warning_count}",
            ]
            + limitations[:7],
            artifact_paths=_artifact_paths(outputs, "report_path"),
        ),
    ]
    return ReviewerAuditChecklist(report_kind="ancestral_report_package", items=items)


def _alignment_checklist(manifest: dict[str, object]) -> ReviewerAuditChecklist:
    output_paths = _text_list(manifest.get("output_paths"))
    metrics = _mapping(manifest, "metrics")
    audit = _mapping(manifest, "audit")
    items = [
        ReviewerAuditChecklistItem(
            section="inputs",
            status="pass",
            summary="alignment input path and checksum were recorded",
            evidence=[
                f"input path: {manifest.get('input_path', '')}",
                f"input checksum: {manifest.get('input_checksum', '')}",
            ],
            artifact_paths=[],
        ),
        ReviewerAuditChecklistItem(
            section="reproducibility",
            status=_status(
                blocked=not bool(manifest.get("reproducibility_manifest_path"))
            ),
            summary="the figure package includes a reproducibility manifest with settings and surfaces",
            evidence=[
                f"reproducibility manifest path: {manifest.get('reproducibility_manifest_path', '')}",
                f"output artifact count: {len(output_paths)}",
            ],
            artifact_paths=[str(manifest.get("reproducibility_manifest_path", ""))],
        ),
        ReviewerAuditChecklistItem(
            section="publication_readiness",
            status=_status(
                blocked=not bool(metrics.get("publication_ready")),
                risk=bool(audit.get("suspicious_alignment")),
            ),
            summary="publication readiness stayed tied to explicit quality and suspicious-alignment checks",
            evidence=[
                f"publication ready: {str(bool(metrics.get('publication_ready'))).lower()}",
                f"quality score: {metrics.get('quality_score', 0)}",
                f"suspicious alignment: {str(bool(audit.get('suspicious_alignment'))).lower()}",
            ]
            + _text_list(audit.get("limitations"))[:5],
            artifact_paths=[
                path for path in output_paths if path.endswith((".html", ".tsv"))
            ][:5],
        ),
        ReviewerAuditChecklistItem(
            section="visible_surfaces",
            status=_status(
                blocked=not all(
                    bool(audit.get(key))
                    for key in (
                        "heatmap_visible",
                        "site_summary_visible",
                        "sequence_panel_visible",
                    )
                )
            ),
            summary="all three reviewer figure surfaces remained visible",
            evidence=[
                f"heatmap visible: {str(bool(audit.get('heatmap_visible'))).lower()}",
                f"site summary visible: {str(bool(audit.get('site_summary_visible'))).lower()}",
                f"sequence panel visible: {str(bool(audit.get('sequence_panel_visible'))).lower()}",
            ],
            artifact_paths=[path for path in output_paths if path.endswith(".svg")][:3],
        ),
        ReviewerAuditChecklistItem(
            section="interpretation_limits",
            status=_status(risk=bool(_text_list(audit.get("limitations")))),
            summary="alignment limitations were surfaced explicitly for reviewer interpretation",
            evidence=_text_list(audit.get("limitations"))[:8],
            artifact_paths=[path for path in output_paths if path.endswith(".html")][
                :1
            ],
        ),
    ]
    return ReviewerAuditChecklist(
        report_kind="alignment_quality_figure_package", items=items
    )


def build_reviewer_audit_checklist(
    manifest: dict[str, object],
) -> ReviewerAuditChecklist:
    """Build one reviewer audit checklist from a supported package manifest."""
    report_kind = str(manifest.get("report_kind", ""))
    if report_kind == "tree_package":
        return _tree_checklist(manifest)
    if report_kind == "comparative_package":
        return _comparative_checklist(manifest)
    if report_kind == "ancestral_report_package":
        return _ancestral_checklist(manifest)
    if report_kind == "alignment_quality_figure_package":
        return _alignment_checklist(manifest)
    raise ValueError("reviewer audit checklist requires a supported package manifest")


def write_reviewer_audit_checklist(
    output_path: Path,
    manifest: dict[str, object],
) -> ReviewerAuditChecklistWriteResult:
    """Write one reviewer-facing audit checklist from a machine manifest."""
    checklist = build_reviewer_audit_checklist(manifest)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["section", "status", "summary", "evidence", "artifact_paths"]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for item in checklist.items:
            writer.writerow(
                {
                    "section": item.section,
                    "status": item.status,
                    "summary": item.summary,
                    "evidence": "|".join(item.evidence),
                    "artifact_paths": "|".join(item.artifact_paths),
                }
            )
    return ReviewerAuditChecklistWriteResult(
        output_path=output_path,
        checklist=checklist,
    )


def write_reviewer_audit_checklist_from_manifest(
    output_path: Path,
    manifest_path: Path,
) -> ReviewerAuditChecklistWriteResult:
    """Write one reviewer audit checklist from a stored package manifest."""
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    return write_reviewer_audit_checklist(output_path, manifest)
