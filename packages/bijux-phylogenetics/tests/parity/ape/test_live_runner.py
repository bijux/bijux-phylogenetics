from __future__ import annotations

from pathlib import Path
import shutil

import pytest

from bijux_phylogenetics.parity import run_ape_parity_cases
from tests.support.live_r_packages import r_package_available


def test_run_ape_parity_cases_records_live_environment_status(tmp_path: Path) -> None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not available")
    if not r_package_available(rscript, "jsonlite"):
        pytest.skip("jsonlite is required for live ape parity validation")

    report = run_ape_parity_cases(
        case_ids=["read-tree-balanced-rooted-ultrametric"],
        rscript_executable=rscript,
        failure_root=tmp_path / "ape-parity-failures",
    )

    observation = report.observations[0]
    if r_package_available(rscript, "ape"):
        assert observation.status == "passed"
        assert observation.ape_version
    else:
        assert observation.status == "skipped"
        assert observation.mismatch_reason == "ape_package_unavailable"
        assert observation.reproducible_artifact_root is not None
