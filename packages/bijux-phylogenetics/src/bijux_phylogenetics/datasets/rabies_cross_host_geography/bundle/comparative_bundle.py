from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.pgls import write_pgls_model_matrix_table
from bijux_phylogenetics.comparative.pgls.categorical_contrasts import (
    write_pgls_categorical_contrast_table,
)
from bijux_phylogenetics.comparative.pgls.lambda_fit import (
    write_pgls_lambda_profile_table,
)
from bijux_phylogenetics.comparative.reporting.analysis_package import (
    ComparativeAnalysisSummaryRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeResidualTableRow,
    ComparativeSignalTableRow,
    summarize_comparative_analysis,
    summarize_comparative_audit,
    summarize_comparative_coefficients,
    summarize_comparative_interpretation,
    summarize_comparative_residuals,
    summarize_comparative_signal,
    write_comparative_audit_table,
    write_comparative_coefficient_table,
    write_comparative_contrast_table,
    write_comparative_interpretation_table,
    write_comparative_model_comparison_table,
    write_comparative_residual_table,
    write_comparative_signal_table,
    write_comparative_summary_table,
)
from bijux_phylogenetics.datasets.rabies_cross_host_geography.models import (
    RabiesCrossHostGeographyPanelWorkflowReport,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .comparative_review import _write_comparative_manifest, _write_comparative_report
from .input_artifacts import _copy_output


@dataclass(frozen=True)
class ComparativeBundleArtifacts:
    output_root: Path
    traits_path: Path
    tree_path: Path
    report_path: Path
    summary_row: ComparativeAnalysisSummaryRow
    coefficient_rows: list[ComparativeCoefficientTableRow]
    residual_rows: list[ComparativeResidualTableRow]
    signal_row: ComparativeSignalTableRow
    interpretation_rows: list[ComparativeInterpretationRow]
    audit_rows: list[dict[str, str]]
    summary_path: Path
    coefficients_path: Path
    residuals_path: Path
    signal_path: Path
    model_comparison_path: Path
    interpretation_path: Path
    audit_path: Path
    contrasts_path: Path
    model_matrix_path: Path
    categorical_contrasts_path: Path
    lambda_profile_path: Path
    manifest_path: Path


def _write_comparative_bundle_artifacts(
    output_root: Path,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
) -> ComparativeBundleArtifacts:
    comparative_traits_path = write_taxon_rows(
        output_root / "comparative-traits.tsv",
        columns=list(report.comparative_traits_rows[0].keys()),
        rows=report.comparative_traits_rows,
    )
    comparative_tree_path = _copy_output(
        report.comparative_tree_path,
        output_root / "comparative-tree.nwk",
    )
    comparative_output_root = output_root / "comparative"
    comparative_output_root.mkdir(parents=True, exist_ok=True)
    comparative_report = report.comparative_report
    comparative_summary_row = summarize_comparative_analysis(comparative_report)
    comparative_coefficient_rows = summarize_comparative_coefficients(
        comparative_report
    )
    comparative_residual_rows = summarize_comparative_residuals(comparative_report)
    comparative_signal_row = summarize_comparative_signal(comparative_report)
    comparative_interpretation_rows = summarize_comparative_interpretation(
        comparative_report
    )
    comparative_audit_rows = summarize_comparative_audit(comparative_report)
    comparative_report_path = _write_comparative_report(
        comparative_output_root / "comparative-report.html",
        summary_row=comparative_summary_row,
        coefficient_rows=comparative_coefficient_rows,
        residual_rows=comparative_residual_rows,
        signal_row=comparative_signal_row,
        interpretation_rows=comparative_interpretation_rows,
        branch_repairs=report.comparative_branch_repairs,
    )
    comparative_summary_path = write_comparative_summary_table(
        comparative_output_root / "comparative-summary.tsv",
        comparative_summary_row,
    )
    comparative_coefficients_path = write_comparative_coefficient_table(
        comparative_output_root / "coefficient-table.tsv",
        comparative_coefficient_rows,
    )
    comparative_residuals_path = write_comparative_residual_table(
        comparative_output_root / "residual-summary.tsv",
        comparative_residual_rows,
    )
    comparative_signal_path = write_comparative_signal_table(
        comparative_output_root / "signal-summary.tsv",
        comparative_signal_row,
    )
    comparative_model_comparison_path = write_comparative_model_comparison_table(
        comparative_output_root / "model-comparison.tsv",
        comparative_report,
    )
    comparative_interpretation_path = write_comparative_interpretation_table(
        comparative_output_root / "interpretation-table.tsv",
        comparative_interpretation_rows,
    )
    comparative_audit_path = write_comparative_audit_table(
        comparative_output_root / "audit-table.tsv",
        comparative_audit_rows,
    )
    comparative_contrasts_path = write_comparative_contrast_table(
        comparative_output_root / "contrast-table.tsv",
        comparative_report,
    )
    comparative_model_matrix_path = comparative_output_root / "model-matrix.tsv"
    write_pgls_model_matrix_table(
        comparative_model_matrix_path,
        comparative_report.snapshot.pgls_inputs.model_matrix,
    )
    comparative_categorical_contrasts_path = write_pgls_categorical_contrast_table(
        comparative_output_root / "categorical-contrasts.tsv",
        report.comparative_categorical_contrasts,
    )
    comparative_lambda_profile_path = write_pgls_lambda_profile_table(
        comparative_output_root / "lambda-profile.tsv",
        comparative_report.snapshot.pgls_model.lambda_fit,
    )
    comparative_manifest_path = _write_comparative_manifest(
        comparative_output_root / "comparative.manifest.json",
        comparative_summary_row=comparative_summary_row,
        branch_repairs=report.comparative_branch_repairs,
        output_paths={
            "comparative_report": comparative_report_path,
            "summary_table": comparative_summary_path,
            "coefficient_table": comparative_coefficients_path,
            "residual_table": comparative_residuals_path,
            "signal_table": comparative_signal_path,
            "model_comparison_table": comparative_model_comparison_path,
            "interpretation_table": comparative_interpretation_path,
            "audit_table": comparative_audit_path,
            "contrast_table": comparative_contrasts_path,
            "model_matrix_table": comparative_model_matrix_path,
            "categorical_contrast_table": comparative_categorical_contrasts_path,
            "lambda_profile_table": comparative_lambda_profile_path,
        },
    )
    return ComparativeBundleArtifacts(
        output_root=comparative_output_root,
        traits_path=comparative_traits_path,
        tree_path=comparative_tree_path,
        report_path=comparative_report_path,
        summary_row=comparative_summary_row,
        coefficient_rows=comparative_coefficient_rows,
        residual_rows=comparative_residual_rows,
        signal_row=comparative_signal_row,
        interpretation_rows=comparative_interpretation_rows,
        audit_rows=comparative_audit_rows,
        summary_path=comparative_summary_path,
        coefficients_path=comparative_coefficients_path,
        residuals_path=comparative_residuals_path,
        signal_path=comparative_signal_path,
        model_comparison_path=comparative_model_comparison_path,
        interpretation_path=comparative_interpretation_path,
        audit_path=comparative_audit_path,
        contrasts_path=comparative_contrasts_path,
        model_matrix_path=comparative_model_matrix_path,
        categorical_contrasts_path=comparative_categorical_contrasts_path,
        lambda_profile_path=comparative_lambda_profile_path,
        manifest_path=comparative_manifest_path,
    )
