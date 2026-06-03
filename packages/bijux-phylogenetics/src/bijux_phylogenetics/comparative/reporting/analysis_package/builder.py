from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path

from bijux_phylogenetics.comparative.reporting import (
    build_comparative_method_report,
    write_comparative_methods_summary_text,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    comparative_report_method_tier,
)
from bijux_phylogenetics.reports.review import write_reviewer_audit_checklist

from .contracts import ComparativeReportPackageResult
from .presentation import write_comparative_report_html
from .summaries import (
    summarize_comparative_analysis,
    summarize_comparative_audit,
    summarize_comparative_coefficients,
    summarize_comparative_interpretation,
    summarize_comparative_residuals,
    summarize_comparative_signal,
)
from .tables import (
    write_comparative_audit_table,
    write_comparative_coefficient_table,
    write_comparative_contrast_table,
    write_comparative_interpretation_table,
    write_comparative_model_comparison_table,
    write_comparative_residual_table,
    write_comparative_signal_table,
    write_comparative_summary_table,
)


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_comparative_report_package(
    tree_path: Path,
    traits_path: Path,
    *,
    out_dir: Path,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeReportPackageResult:
    out_dir.mkdir(parents=True, exist_ok=True)
    report = build_comparative_method_report(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    summary_row = summarize_comparative_analysis(report)
    coefficient_rows = summarize_comparative_coefficients(report)
    residual_rows = summarize_comparative_residuals(report)
    signal_row = summarize_comparative_signal(report)
    interpretation_rows = summarize_comparative_interpretation(report)
    audit_rows = summarize_comparative_audit(report)
    method_tier = comparative_report_method_tier()

    report_path = out_dir / "comparative-report.html"
    methods_summary_path = out_dir / "comparative-methods-summary.md"
    reviewer_audit_checklist_path = out_dir / "reviewer-audit-checklist.tsv"
    summary_table_path = out_dir / "comparative-summary.tsv"
    coefficient_table_path = out_dir / "coefficient-table.tsv"
    residual_table_path = out_dir / "residual-summary.tsv"
    signal_table_path = out_dir / "signal-summary.tsv"
    model_comparison_table_path = out_dir / "model-comparison.tsv"
    interpretation_table_path = out_dir / "interpretation-table.tsv"
    audit_table_path = out_dir / "audit-table.tsv"
    contrast_table_path = out_dir / "contrast-table.tsv"
    manifest_path = out_dir / "comparative-report.manifest.json"

    methods_summary = write_comparative_methods_summary_text(
        methods_summary_path, report
    )
    write_comparative_summary_table(summary_table_path, summary_row)
    write_comparative_coefficient_table(coefficient_table_path, coefficient_rows)
    write_comparative_residual_table(residual_table_path, residual_rows)
    write_comparative_signal_table(signal_table_path, signal_row)
    write_comparative_model_comparison_table(model_comparison_table_path, report)
    write_comparative_interpretation_table(
        interpretation_table_path, interpretation_rows
    )
    write_comparative_audit_table(audit_table_path, audit_rows)
    write_comparative_contrast_table(contrast_table_path, report)

    machine_manifest = {
        "report_kind": "comparative_package",
        "input_paths": [str(tree_path), str(traits_path)],
        "input_checksums": {
            str(tree_path): _checksum(tree_path),
            str(traits_path): _checksum(traits_path),
        },
        "outputs": {
            "report_path": str(report_path),
            "methods_summary_path": str(methods_summary_path),
            "reviewer_audit_checklist_path": str(reviewer_audit_checklist_path),
            "summary_table_path": str(summary_table_path),
            "coefficient_table_path": str(coefficient_table_path),
            "residual_table_path": str(residual_table_path),
            "signal_table_path": str(signal_table_path),
            "model_comparison_table_path": str(model_comparison_table_path),
            "interpretation_table_path": str(interpretation_table_path),
            "audit_table_path": str(audit_table_path),
            "contrast_table_path": str(contrast_table_path),
        },
        "metrics": {
            "analysis_taxa": summary_row.analysis_taxa,
            "selected_model": summary_row.selected_model,
            "coefficient_count": len(coefficient_rows),
            "contrast_count": signal_row.independent_contrast_count,
            "limitation_count": len(report.snapshot.limitations),
            "methods_summary_warning_count": methods_summary.warning_count,
        },
        "summary": asdict(summary_row),
        "limitations": report.snapshot.limitations,
    }
    reviewer_audit_checklist = write_reviewer_audit_checklist(
        reviewer_audit_checklist_path,
        machine_manifest,
    ).checklist
    machine_manifest["reviewer_audit_checklist"] = asdict(reviewer_audit_checklist)
    manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_comparative_report_html(
        path=report_path,
        report=report,
        methods_summary_text=methods_summary.text,
        summary_row=summary_row,
        coefficient_rows=coefficient_rows,
        residual_rows=residual_rows,
        signal_row=signal_row,
        interpretation_rows=interpretation_rows,
        method_tier=method_tier,
        reviewer_audit_checklist=reviewer_audit_checklist,
        manifest=machine_manifest,
    )
    return ComparativeReportPackageResult(
        output_dir=out_dir,
        report_path=report_path,
        methods_summary_path=methods_summary_path,
        reviewer_audit_checklist_path=reviewer_audit_checklist_path,
        summary_table_path=summary_table_path,
        coefficient_table_path=coefficient_table_path,
        residual_table_path=residual_table_path,
        signal_table_path=signal_table_path,
        model_comparison_table_path=model_comparison_table_path,
        interpretation_table_path=interpretation_table_path,
        audit_table_path=audit_table_path,
        contrast_table_path=contrast_table_path,
        manifest_path=manifest_path,
        report=report,
        methods_summary=methods_summary,
        summary_row=summary_row,
        coefficient_rows=coefficient_rows,
        residual_rows=residual_rows,
        signal_row=signal_row,
        interpretation_rows=interpretation_rows,
        audit_rows=audit_rows,
        method_tier=method_tier,
        reviewer_audit_checklist=reviewer_audit_checklist,
        machine_manifest=machine_manifest,
    )
