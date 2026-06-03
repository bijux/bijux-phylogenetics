from __future__ import annotations

from .builder import build_generated_geiger_parity_report
from .contracts import (
    GeigerBenchmarkSummaryRow,
    GeigerBoundaryWarningSummaryRow,
    GeigerExcludedModelRow,
    GeigerGoalCoverageRow,
    GeigerOptimizerMismatchCategoryRow,
    GeigerSimulationRecoveryRow,
    GeigerToleranceRuleRow,
    GeneratedGeigerParityReport,
)
from .presentation import (
    write_generated_geiger_parity_report_json,
    write_generated_geiger_parity_report_markdown,
)

__all__ = [
    "build_generated_geiger_parity_report",
    "GeigerBenchmarkSummaryRow",
    "GeigerBoundaryWarningSummaryRow",
    "GeigerExcludedModelRow",
    "GeigerGoalCoverageRow",
    "GeigerOptimizerMismatchCategoryRow",
    "GeigerSimulationRecoveryRow",
    "GeigerToleranceRuleRow",
    "GeneratedGeigerParityReport",
    "write_generated_geiger_parity_report_json",
    "write_generated_geiger_parity_report_markdown",
]
