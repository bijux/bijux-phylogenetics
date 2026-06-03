from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.comparative.continuous_mode_recovery.reference_payloads import (
    GEIGER_FITCONTINUOUS_RECOVERY_REFERENCE_PAYLOADS,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .models import ContinuousModeRecoveryReport


def geiger_fitcontinuous_recovery_reference_payload(case_id: str) -> dict[str, object]:
    """Expose the governed stored geiger recovery payload for one recovery case."""
    return GEIGER_FITCONTINUOUS_RECOVERY_REFERENCE_PAYLOADS[case_id]


def write_geiger_fitcontinuous_recovery_reference_payload_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write the stored governed geiger payload summaries used in the benchmark."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "fit_summary_json",
            "comparison_summary_json",
        ],
        rows=[
            {
                "case_id": case.scenario.case_id,
                "fit_summary_json": json.dumps(
                    GEIGER_FITCONTINUOUS_RECOVERY_REFERENCE_PAYLOADS[
                        case.scenario.case_id
                    ]["fit_summary"],
                    sort_keys=True,
                ),
                "comparison_summary_json": json.dumps(
                    GEIGER_FITCONTINUOUS_RECOVERY_REFERENCE_PAYLOADS[
                        case.scenario.case_id
                    ]["comparison_summary"],
                    sort_keys=True,
                ),
            }
            for case in report.case_reports
        ],
    )
