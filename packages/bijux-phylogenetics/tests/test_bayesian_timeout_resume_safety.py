from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import json
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.beast import (
    prepare_beast_time_tree_analysis,
    run_beast_posterior_inference,
)
from bijux_phylogenetics.bayesian.mrbayes import (
    prepare_mrbayes_analysis,
    run_mrbayes_posterior_inference,
)
from bijux_phylogenetics.runtime.errors import (
    EngineUnavailableError,
    EngineWorkflowError,
)

from .support.fake_bayesian_engines import (
    fake_beast,
    fake_beast_killed,
    fake_beast_malformed_outputs,
    fake_beast_timeout,
    fake_mrbayes,
    fake_mrbayes_killed,
    fake_mrbayes_malformed_outputs,
    fake_mrbayes_timeout,
)

pytestmark = pytest.mark.engine_contract

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


@dataclass(frozen=True, slots=True)
class BayesianPosteriorSafetyCase:
    name: str
    prepare_analysis: Callable[[Path], Path]
    run_inference: Callable[..., object]
    fake_valid: Callable[[Path], Path]
    fake_timeout: Callable[[Path], Path]
    fake_killed: Callable[[Path], Path]
    fake_malformed: Callable[[Path], Path]
    manifest_path: Callable[[Path], Path]


def _prepare_beast_analysis(tmp_path: Path) -> Path:
    xml_path = tmp_path / "analysis.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        xml_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )
    return xml_path


def _prepare_mrbayes_analysis(tmp_path: Path) -> Path:
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(
        fixture("alignments/example_alignment.fasta"),
        nexus_path,
    )
    return nexus_path


SAFETY_CASES = (
    BayesianPosteriorSafetyCase(
        name="beast",
        prepare_analysis=_prepare_beast_analysis,
        run_inference=run_beast_posterior_inference,
        fake_valid=fake_beast,
        fake_timeout=fake_beast_timeout,
        fake_killed=fake_beast_killed,
        fake_malformed=fake_beast_malformed_outputs,
        manifest_path=lambda analysis_path: analysis_path.with_suffix(".manifest.json"),
    ),
    BayesianPosteriorSafetyCase(
        name="mrbayes",
        prepare_analysis=_prepare_mrbayes_analysis,
        run_inference=run_mrbayes_posterior_inference,
        fake_valid=fake_mrbayes,
        fake_timeout=fake_mrbayes_timeout,
        fake_killed=fake_mrbayes_killed,
        fake_malformed=fake_mrbayes_malformed_outputs,
        manifest_path=lambda analysis_path: analysis_path.with_suffix("").with_suffix(
            ".manifest.json"
        ),
    ),
)


@pytest.mark.parametrize("case", SAFETY_CASES, ids=lambda case: case.name)
def test_bayesian_posterior_timeout_marks_incomplete(
    tmp_path: Path,
    case: BayesianPosteriorSafetyCase,
) -> None:
    analysis_path = case.prepare_analysis(tmp_path)
    executable = case.fake_timeout(tmp_path / f"{case.name}-timeout")

    with pytest.raises(EngineWorkflowError, match="timed out"):
        case.run_inference(
            analysis_path,
            executable=executable,
            timeout_seconds=0.5,
        )

    marker_path = case.manifest_path(analysis_path).with_suffix(".incomplete.json")
    assert marker_path.exists()
    payload = json.loads(marker_path.read_text(encoding="utf-8"))
    assert payload["timed_out"] is True


@pytest.mark.parametrize("case", SAFETY_CASES, ids=lambda case: case.name)
def test_bayesian_posterior_killed_process_marks_incomplete(
    tmp_path: Path,
    case: BayesianPosteriorSafetyCase,
) -> None:
    analysis_path = case.prepare_analysis(tmp_path)
    executable = case.fake_killed(tmp_path / f"{case.name}-killed")

    with pytest.raises(EngineWorkflowError, match="failed with exit code"):
        case.run_inference(
            analysis_path,
            executable=executable,
        )

    marker_path = case.manifest_path(analysis_path).with_suffix(".incomplete.json")
    payload = json.loads(marker_path.read_text(encoding="utf-8"))
    assert payload["exit_code"] == -15


@pytest.mark.parametrize("case", SAFETY_CASES, ids=lambda case: case.name)
def test_bayesian_posterior_missing_executable_leaves_no_incomplete_marker(
    tmp_path: Path,
    case: BayesianPosteriorSafetyCase,
) -> None:
    analysis_path = case.prepare_analysis(tmp_path)

    with pytest.raises(EngineUnavailableError, match="was not found"):
        case.run_inference(
            analysis_path,
            executable=tmp_path / f"missing-{case.name}",
        )

    marker_path = case.manifest_path(analysis_path).with_suffix(".incomplete.json")
    assert marker_path.exists() is False


@pytest.mark.parametrize("case", SAFETY_CASES, ids=lambda case: case.name)
def test_bayesian_posterior_incomplete_outputs_are_rejected_until_cleaned(
    tmp_path: Path,
    case: BayesianPosteriorSafetyCase,
) -> None:
    analysis_path = case.prepare_analysis(tmp_path)
    malformed = case.fake_malformed(tmp_path / f"{case.name}-malformed")
    valid = case.fake_valid(tmp_path / f"{case.name}-valid")

    with pytest.raises(EngineWorkflowError):
        case.run_inference(
            analysis_path,
            executable=malformed,
        )

    marker_path = case.manifest_path(analysis_path).with_suffix(".incomplete.json")
    assert marker_path.exists()

    with pytest.raises(EngineWorkflowError, match="incomplete outputs") as rejected:
        case.run_inference(
            analysis_path,
            executable=valid,
            resume=True,
            incomplete_run_policy="reject",
        )

    assert rejected.value.code == "engine_incomplete_outputs_present"

    cleaned = case.run_inference(
        analysis_path,
        executable=valid,
        resume=True,
        incomplete_run_policy="clean",
    )
    assert cleaned.resumed is False
    assert marker_path.exists() is False


@pytest.mark.parametrize("case", SAFETY_CASES, ids=lambda case: case.name)
def test_bayesian_posterior_resume_reuses_only_verified_complete_outputs(
    tmp_path: Path,
    case: BayesianPosteriorSafetyCase,
) -> None:
    analysis_path = case.prepare_analysis(tmp_path)
    executable = case.fake_valid(tmp_path / f"{case.name}-valid")

    first = case.run_inference(
        analysis_path,
        executable=executable,
    )
    resumed = case.run_inference(
        analysis_path,
        executable=executable,
        resume=True,
    )

    assert first.resumed is False
    assert resumed.resumed is True
