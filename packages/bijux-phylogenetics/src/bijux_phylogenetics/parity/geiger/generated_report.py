from __future__ import annotations

from bijux_phylogenetics.parity.geiger.generated_parity_report import (
    GeigerBenchmarkSummaryRow,
    GeigerBoundaryWarningSummaryRow,
    GeigerExcludedModelRow,
    GeigerGoalCoverageRow,
    GeigerOptimizerMismatchCategoryRow,
    GeigerSimulationRecoveryRow,
    GeigerToleranceRuleRow,
    GeneratedGeigerParityReport,
    build_generated_geiger_parity_report,
    write_generated_geiger_parity_report_json,
    write_generated_geiger_parity_report_markdown,
)
from bijux_phylogenetics.parity.geiger.generated_parity_report.governed_artifacts import (
    load_large_tree_benchmark_summary as _load_large_tree_benchmark_summary,
)
from bijux_phylogenetics.parity.geiger.generated_parity_report.governed_artifacts import (
    load_real_dataset_benchmark_summary as _load_real_dataset_benchmark_summary,
)
from bijux_phylogenetics.parity.geiger.generated_parity_report.governed_artifacts import (
    load_recovery_summary as _load_recovery_summary,
)
from bijux_phylogenetics.parity.geiger.generated_parity_report.governed_artifacts import (
    load_sim_char_summary as _load_sim_char_summary,
)
from bijux_phylogenetics.parity.geiger.generated_parity_report.governed_artifacts import (
    repository_root as _repository_root,
)

# Preserve the compatibility seam that governed tests patch on the public facade.
_COMPATIBILITY_TEST_SEAMS = (
    _load_recovery_summary,
    _load_sim_char_summary,
    _load_large_tree_benchmark_summary,
    _load_real_dataset_benchmark_summary,
    _repository_root,
)

__all__ = [
    "GeneratedGeigerParityReport",
    "GeigerGoalCoverageRow",
    "GeigerExcludedModelRow",
    "GeigerToleranceRuleRow",
    "GeigerOptimizerMismatchCategoryRow",
    "GeigerBoundaryWarningSummaryRow",
    "GeigerSimulationRecoveryRow",
    "GeigerBenchmarkSummaryRow",
    "build_generated_geiger_parity_report",
    "write_generated_geiger_parity_report_json",
    "write_generated_geiger_parity_report_markdown",
]
