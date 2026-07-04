from __future__ import annotations

from pathlib import Path
import shutil

import pytest

from bijux_phylogenetics.parity import run_phytools_parity_cases
from tests.support.live_r_packages import r_package_available


@pytest.mark.slow
def test_run_phytools_parity_cases_passes_sim_history_cases_against_live_phytools_when_available(
    tmp_path: Path,
) -> None:
    rscript = shutil.which("Rscript")
    if rscript is None or not r_package_available(rscript, "phytools"):
        pytest.skip("live phytools package is not available")

    report = run_phytools_parity_cases(
        case_ids=[
            "sim-history-binary-no-change-example-tree",
            "sim-history-binary-high-rate-example-tree",
            "sim-history-multistate-no-change-six-taxa",
            "sim-history-multistate-high-rate-six-taxa",
            "sim-history-binary-root-prior-example-tree",
        ],
        rscript_executable=rscript,
        failure_root=tmp_path / "phytools-live-sim-history-failures",
    )

    assert report.all_passed is True
    assert report.case_count == 5
    assert report.failed_case_count == 0
    assert report.skipped_case_count == 0
