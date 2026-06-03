from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.engines.validation.preflight import (
    inspect_external_engine_preflight,
    inspect_external_engine_surface,
    require_external_engine_surface,
    require_preflight_workflow,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError


@pytest.fixture
def fake_engine_metadata(monkeypatch: pytest.MonkeyPatch):
    installed = {
        "mafft": ("MAFFT v7.520", "/opt/tools/mafft"),
        "trimal": ("trimAl v1.5.0", "/opt/tools/trimal"),
        "iqtree2": ("IQ-TREE multicore version 2.4.0", "/opt/tools/iqtree2"),
        "FastTree": ("FastTree Version 2.2.0", "/opt/tools/FastTree"),
        "mb": ("MrBayes v3.2.7a", "/opt/tools/mb"),
        "beast": ("BEAST v2.7.7", "/opt/tools/beast"),
    }

    def fake_resolve(executable: str | Path) -> str:
        key = str(executable)
        if key not in installed:
            raise Exception(f"unexpected executable lookup: {key}")
        return installed[key][1]

    def fake_read(engine_name: str, executable: str | Path, *, version_args):
        key = (
            Path(executable).name
            if Path(executable).name in installed
            else str(executable)
        )
        if key == "mb":
            assert version_args == ("-v",)
        version_text = installed[key][0]
        from bijux_phylogenetics.engines.common import EngineVersionInfo

        return EngineVersionInfo(
            engine_name=engine_name,
            executable=str(executable),
            command=[str(executable), *version_args],
            text=version_text,
        )

    monkeypatch.setattr(
        "bijux_phylogenetics.engines.validation.preflight.resolve_engine_executable",
        fake_resolve,
    )
    monkeypatch.setattr(
        "bijux_phylogenetics.engines.validation.preflight.read_engine_version",
        fake_read,
    )
    return installed


def test_engine_preflight_marks_supported_versions_ready(fake_engine_metadata) -> None:
    report = inspect_external_engine_preflight()

    assert report.overall_status == "ready"
    assert all(engine.support_status == "tested" for engine in report.engines)
    fasta_to_tree = next(
        workflow
        for workflow in report.workflows
        if workflow.workflow_id == "fasta-to-tree"
    )
    assert fasta_to_tree.runnable is True
    assert fasta_to_tree.readiness_status == "ready"


def test_engine_preflight_marks_newer_versions_as_caution(fake_engine_metadata) -> None:
    fake_engine_metadata["iqtree2"] = (
        "IQ-TREE multicore version 4.0.0",
        "/opt/tools/iqtree2",
    )

    report = inspect_external_engine_preflight(selected_workflow="fasta-to-tree")

    iqtree = next(engine for engine in report.engines if engine.engine_id == "iqtree")
    workflow = next(
        workflow
        for workflow in report.workflows
        if workflow.workflow_id == "fasta-to-tree"
    )
    assert iqtree.support_status == "untested"
    assert iqtree.compatible is True
    assert workflow.readiness_status == "caution"
    assert workflow.runnable is True
    assert workflow.caution_engines == ["IQ-TREE"]


def test_engine_preflight_blocks_missing_selected_workflow_dependencies(
    fake_engine_metadata,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_resolve(executable: str | Path) -> str:
        if str(executable) == "beast":
            from bijux_phylogenetics.runtime.errors import EngineUnavailableError

            raise EngineUnavailableError(
                "engine executable is not available on PATH: beast"
            )
        return fake_engine_metadata[str(executable)][1]

    monkeypatch.setattr(
        "bijux_phylogenetics.engines.validation.preflight.resolve_engine_executable",
        fake_resolve,
    )

    report = inspect_external_engine_preflight(selected_workflow="beast-posterior")
    workflow = next(
        workflow
        for workflow in report.workflows
        if workflow.workflow_id == "beast-posterior"
    )

    assert workflow.runnable is False
    assert workflow.readiness_status == "blocked"
    assert workflow.blocking_engines == ["BEAST"]
    with pytest.raises(
        EngineWorkflowError, match="workflow 'beast-posterior' is blocked"
    ) as error:
        require_preflight_workflow(report, workflow_id="beast-posterior")
    assert error.value.code == "engine_preflight_workflow_blocked"


def test_engine_preflight_marks_older_versions_unsupported(
    fake_engine_metadata,
) -> None:
    fake_engine_metadata["mb"] = ("MrBayes v3.1.2", "/opt/tools/mb")

    report = inspect_external_engine_preflight(selected_workflow="mrbayes-posterior")

    mrbayes = next(engine for engine in report.engines if engine.engine_id == "mrbayes")
    workflow = next(
        workflow
        for workflow in report.workflows
        if workflow.workflow_id == "mrbayes-posterior"
    )
    assert mrbayes.support_status == "unsupported"
    assert mrbayes.compatible is False
    assert workflow.runnable is False
    assert workflow.blocking_engines == ["MrBayes"]


def test_engine_preflight_inspects_direct_surface_requirements(
    fake_engine_metadata,
) -> None:
    statuses, workflow = inspect_external_engine_surface(
        workflow_id="maximum-likelihood-tree",
        summary="IQ-TREE maximum-likelihood tree-inference workflow.",
        required_engines=("iqtree",),
        executables={"iqtree": "iqtree2"},
    )

    assert [status.engine_id for status in statuses] == ["iqtree"]
    assert statuses[0].support_status == "tested"
    assert workflow.workflow_id == "maximum-likelihood-tree"
    assert workflow.required_engines == ["IQ-TREE"]
    assert workflow.readiness_status == "ready"


def test_engine_preflight_blocks_direct_surface_when_engine_is_unsupported(
    fake_engine_metadata,
) -> None:
    fake_engine_metadata["iqtree2"] = (
        "IQ-TREE multicore version 1.6.0",
        "/opt/tools/iqtree2",
    )

    with pytest.raises(
        EngineWorkflowError,
        match="workflow 'maximum-likelihood-tree' is blocked",
    ) as error:
        require_external_engine_surface(
            workflow_id="maximum-likelihood-tree",
            summary="IQ-TREE maximum-likelihood tree-inference workflow.",
            required_engines=("iqtree",),
            executables={"iqtree": "iqtree2"},
        )

    assert error.value.code == "engine_preflight_workflow_blocked"
    assert error.value.details["blocking_engines"] == ["IQ-TREE"]
