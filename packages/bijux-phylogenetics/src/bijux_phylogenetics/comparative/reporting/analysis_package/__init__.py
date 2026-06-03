from __future__ import annotations

from .builder import build_comparative_report_package
from .contracts import (
    ComparativeAnalysisSummaryRow,
    ComparativeAuditTableRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeReportPackageResult,
    ComparativeResidualTableRow,
    ComparativeSignalTableRow,
)
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
