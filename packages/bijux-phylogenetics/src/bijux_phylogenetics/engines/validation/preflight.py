from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re

from bijux_phylogenetics.runtime.errors import (
    EngineUnavailableError,
    EngineWorkflowError,
)

from ..common import EngineVersionInfo, read_engine_version, resolve_engine_executable

VersionTuple = tuple[int, ...]


@dataclass(frozen=True, slots=True)
class EngineSupportPolicy:
    engine_id: str
    display_name: str
    default_executable: str
    version_args: tuple[str, ...]
    minimum_supported: VersionTuple
    tested_before: VersionTuple


@dataclass(frozen=True, slots=True)
class WorkflowSupportPolicy:
    workflow_id: str
    summary: str
    required_engines: tuple[str, ...]


@dataclass(slots=True)
class ExternalEnginePreflightStatus:
    engine_id: str
    engine_name: str
    requested_executable: str
    executable_path: Path | None
    available: bool
    version_text: str | None
    parsed_version: str | None
    support_status: str
    compatible: bool
    issues: list[str] = field(default_factory=list)


@dataclass(slots=True)
class WorkflowPreflightStatus:
    workflow_id: str
    summary: str
    required_engines: list[str]
    readiness_status: str
    runnable: bool
    blocking_engines: list[str]
    caution_engines: list[str]
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ExternalEnginePreflightReport:
    engines: list[ExternalEnginePreflightStatus]
    workflows: list[WorkflowPreflightStatus]
    selected_workflow: str | None
    overall_status: str


_VERSION_PATTERN = re.compile(r"(\d+)\.(\d+)(?:\.(\d+))?")

_ENGINE_POLICIES: tuple[EngineSupportPolicy, ...] = (
    EngineSupportPolicy(
        engine_id="mafft",
        display_name="MAFFT",
        default_executable="mafft",
        version_args=("--version",),
        minimum_supported=(7, 0, 0),
        tested_before=(8, 0, 0),
    ),
    EngineSupportPolicy(
        engine_id="trimal",
        display_name="trimAl",
        default_executable="trimal",
        version_args=("--version",),
        minimum_supported=(1, 0, 0),
        tested_before=(2, 0, 0),
    ),
    EngineSupportPolicy(
        engine_id="iqtree",
        display_name="IQ-TREE",
        default_executable="iqtree2",
        version_args=("--version",),
        minimum_supported=(2, 0, 0),
        tested_before=(4, 0, 0),
    ),
    EngineSupportPolicy(
        engine_id="fasttree",
        display_name="FastTree",
        default_executable="FastTree",
        version_args=("-help",),
        minimum_supported=(2, 0, 0),
        tested_before=(3, 0, 0),
    ),
    EngineSupportPolicy(
        engine_id="mrbayes",
        display_name="MrBayes",
        default_executable="mb",
        version_args=("-v",),
        minimum_supported=(3, 2, 0),
        tested_before=(4, 0, 0),
    ),
    EngineSupportPolicy(
        engine_id="beast",
        display_name="BEAST",
        default_executable="beast",
        version_args=("-version",),
        minimum_supported=(2, 0, 0),
        tested_before=(3, 0, 0),
    ),
)

_WORKFLOW_POLICIES: tuple[WorkflowSupportPolicy, ...] = (
    WorkflowSupportPolicy(
        workflow_id="fasta-to-tree",
        summary="Canonical unaligned FASTA to trimmed alignment plus IQ-TREE inference workflow.",
        required_engines=("mafft", "trimal", "iqtree"),
    ),
    WorkflowSupportPolicy(
        workflow_id="model-selection",
        summary="IQ-TREE model-selection workflow.",
        required_engines=("iqtree",),
    ),
    WorkflowSupportPolicy(
        workflow_id="maximum-likelihood",
        summary="IQ-TREE maximum-likelihood inference workflow.",
        required_engines=("iqtree",),
    ),
    WorkflowSupportPolicy(
        workflow_id="bootstrap-support",
        summary="IQ-TREE bootstrap support workflow.",
        required_engines=("iqtree",),
    ),
    WorkflowSupportPolicy(
        workflow_id="fasttree-inference",
        summary="FastTree approximate inference workflow.",
        required_engines=("fasttree",),
    ),
    WorkflowSupportPolicy(
        workflow_id="tree-inference-comparison",
        summary="Paired IQ-TREE and FastTree comparison workflow.",
        required_engines=("iqtree", "fasttree"),
    ),
    WorkflowSupportPolicy(
        workflow_id="mrbayes-posterior",
        summary="MrBayes posterior inference workflow.",
        required_engines=("mrbayes",),
    ),
    WorkflowSupportPolicy(
        workflow_id="beast-posterior",
        summary="BEAST posterior inference workflow.",
        required_engines=("beast",),
    ),
)


def list_external_engine_workflows() -> tuple[str, ...]:
    """Return the stable workflow ids supported by the engine preflight surface."""
    return tuple(policy.workflow_id for policy in _WORKFLOW_POLICIES)


def _policy_by_engine_id(engine_id: str) -> EngineSupportPolicy:
    for policy in _ENGINE_POLICIES:
        if policy.engine_id == engine_id:
            return policy
    raise KeyError(engine_id)


def _parse_version_tuple(text: str) -> VersionTuple | None:
    match = _VERSION_PATTERN.search(text)
    if match is None:
        return None
    groups = [int(value) for value in match.groups(default="0")]
    return tuple(groups)


def _version_text(version: VersionTuple | None) -> str | None:
    if version is None:
        return None
    return ".".join(str(part) for part in version)


def _classify_engine_support(
    policy: EngineSupportPolicy,
    *,
    resolved_path: Path,
    version: EngineVersionInfo,
    requested_executable: str,
) -> ExternalEnginePreflightStatus:
    parsed_version = _parse_version_tuple(version.text)
    issues: list[str] = []
    if parsed_version is None:
        issues.append(
            f"{policy.display_name} version could not be parsed from '{version.text.splitlines()[0]}'"
        )
        return ExternalEnginePreflightStatus(
            engine_id=policy.engine_id,
            engine_name=policy.display_name,
            requested_executable=requested_executable,
            executable_path=resolved_path,
            available=True,
            version_text=version.text,
            parsed_version=None,
            support_status="untested",
            compatible=True,
            issues=issues,
        )
    if parsed_version < policy.minimum_supported:
        issues.append(
            f"{policy.display_name} {'.'.join(str(part) for part in parsed_version)} is older than the supported floor {'.'.join(str(part) for part in policy.minimum_supported)}"
        )
        return ExternalEnginePreflightStatus(
            engine_id=policy.engine_id,
            engine_name=policy.display_name,
            requested_executable=requested_executable,
            executable_path=resolved_path,
            available=True,
            version_text=version.text,
            parsed_version=_version_text(parsed_version),
            support_status="unsupported",
            compatible=False,
            issues=issues,
        )
    if parsed_version >= policy.tested_before:
        issues.append(
            f"{policy.display_name} {'.'.join(str(part) for part in parsed_version)} is newer than the governed tested range before {'.'.join(str(part) for part in policy.tested_before)}"
        )
        return ExternalEnginePreflightStatus(
            engine_id=policy.engine_id,
            engine_name=policy.display_name,
            requested_executable=requested_executable,
            executable_path=resolved_path,
            available=True,
            version_text=version.text,
            parsed_version=_version_text(parsed_version),
            support_status="untested",
            compatible=True,
            issues=issues,
        )
    return ExternalEnginePreflightStatus(
        engine_id=policy.engine_id,
        engine_name=policy.display_name,
        requested_executable=requested_executable,
        executable_path=resolved_path,
        available=True,
        version_text=version.text,
        parsed_version=_version_text(parsed_version),
        support_status="tested",
        compatible=True,
        issues=[],
    )


def _inspect_engine_status(
    policy: EngineSupportPolicy,
    *,
    requested_executable: str,
) -> ExternalEnginePreflightStatus:
    try:
        resolved = Path(resolve_engine_executable(requested_executable))
        version = read_engine_version(
            policy.display_name,
            resolved,
            version_args=policy.version_args,
        )
    except EngineUnavailableError as error:
        return ExternalEnginePreflightStatus(
            engine_id=policy.engine_id,
            engine_name=policy.display_name,
            requested_executable=requested_executable,
            executable_path=None,
            available=False,
            version_text=None,
            parsed_version=None,
            support_status="missing",
            compatible=False,
            issues=[error.message],
        )
    return _classify_engine_support(
        policy,
        resolved_path=resolved,
        version=version,
        requested_executable=requested_executable,
    )


def _build_workflow_preflight_status(
    policy: WorkflowSupportPolicy,
    *,
    by_engine_id: dict[str, ExternalEnginePreflightStatus],
) -> WorkflowPreflightStatus:
    blocking_engines = [
        engine_id
        for engine_id in policy.required_engines
        if by_engine_id[engine_id].support_status in {"missing", "unsupported"}
    ]
    caution_engines = [
        engine_id
        for engine_id in policy.required_engines
        if by_engine_id[engine_id].support_status == "untested"
    ]
    if blocking_engines:
        readiness_status = "blocked"
        notes = [
            f"{_policy_by_engine_id(engine_id).display_name}: {by_engine_id[engine_id].issues[0]}"
            for engine_id in blocking_engines
        ]
    elif caution_engines:
        readiness_status = "caution"
        notes = [
            f"{_policy_by_engine_id(engine_id).display_name}: {by_engine_id[engine_id].issues[0]}"
            for engine_id in caution_engines
        ]
    else:
        readiness_status = "ready"
        notes = []
    return WorkflowPreflightStatus(
        workflow_id=policy.workflow_id,
        summary=policy.summary,
        required_engines=[
            _policy_by_engine_id(engine_id).display_name
            for engine_id in policy.required_engines
        ],
        readiness_status=readiness_status,
        runnable=readiness_status != "blocked",
        blocking_engines=[
            _policy_by_engine_id(engine_id).display_name
            for engine_id in blocking_engines
        ],
        caution_engines=[
            _policy_by_engine_id(engine_id).display_name
            for engine_id in caution_engines
        ],
        notes=notes,
    )


def inspect_external_engine_surface(
    *,
    workflow_id: str,
    summary: str,
    required_engines: tuple[str, ...],
    executables: dict[str, str | Path | None] | None = None,
) -> tuple[list[ExternalEnginePreflightStatus], WorkflowPreflightStatus]:
    """Inspect one concrete external-engine workflow surface without scanning every engine."""
    configured = {} if executables is None else dict(executables)
    statuses: list[ExternalEnginePreflightStatus] = []
    by_engine_id: dict[str, ExternalEnginePreflightStatus] = {}
    for engine_id in required_engines:
        policy = _policy_by_engine_id(engine_id)
        requested = str(configured.get(engine_id) or policy.default_executable)
        status = _inspect_engine_status(policy, requested_executable=requested)
        statuses.append(status)
        by_engine_id[engine_id] = status
    workflow_status = _build_workflow_preflight_status(
        WorkflowSupportPolicy(
            workflow_id=workflow_id,
            summary=summary,
            required_engines=required_engines,
        ),
        by_engine_id=by_engine_id,
    )
    return statuses, workflow_status


def inspect_external_engine_preflight(
    *,
    executables: dict[str, str | Path | None] | None = None,
    selected_workflow: str | None = None,
) -> ExternalEnginePreflightReport:
    """Inspect external engine availability and workflow readiness."""
    configured = {} if executables is None else dict(executables)
    statuses: list[ExternalEnginePreflightStatus] = []
    by_engine_id: dict[str, ExternalEnginePreflightStatus] = {}
    for policy in _ENGINE_POLICIES:
        requested = str(configured.get(policy.engine_id) or policy.default_executable)
        status = _inspect_engine_status(policy, requested_executable=requested)
        statuses.append(status)
        by_engine_id[policy.engine_id] = status

    workflows: list[WorkflowPreflightStatus] = []
    selected_seen = selected_workflow is None
    for policy in _WORKFLOW_POLICIES:
        if selected_workflow == policy.workflow_id:
            selected_seen = True
        workflows.append(
            _build_workflow_preflight_status(policy, by_engine_id=by_engine_id)
        )
    if not selected_seen:
        available = ", ".join(list_external_engine_workflows())
        raise ValueError(
            f"selected_workflow must be one of: {available}; got {selected_workflow}"
        )
    statuses_by_readiness = {workflow.readiness_status for workflow in workflows}
    overall_status = (
        "blocked"
        if "blocked" in statuses_by_readiness
        else "caution"
        if "caution" in statuses_by_readiness
        else "ready"
    )
    return ExternalEnginePreflightReport(
        engines=statuses,
        workflows=workflows,
        selected_workflow=selected_workflow,
        overall_status=overall_status,
    )


def require_preflight_workflow(
    report: ExternalEnginePreflightReport,
    *,
    workflow_id: str,
) -> WorkflowPreflightStatus:
    """Return one selected workflow or raise when the environment blocks it."""
    for workflow in report.workflows:
        if workflow.workflow_id != workflow_id:
            continue
        if workflow.runnable:
            return workflow
        raise EngineWorkflowError(
            f"workflow '{workflow.workflow_id}' is blocked by external engine availability or compatibility",
            code="engine_preflight_workflow_blocked",
            details={
                "workflow_id": workflow.workflow_id,
                "blocking_engines": workflow.blocking_engines,
                "required_engines": workflow.required_engines,
                "notes": workflow.notes,
            },
        )
    available = ", ".join(policy.workflow_id for policy in _WORKFLOW_POLICIES)
    raise ValueError(f"unknown workflow '{workflow_id}', expected one of: {available}")


def require_external_engine_surface(
    *,
    workflow_id: str,
    summary: str,
    required_engines: tuple[str, ...],
    executables: dict[str, str | Path | None] | None = None,
    preserve_missing_error: bool = False,
) -> WorkflowPreflightStatus:
    """Require one concrete external-engine workflow surface to be runnable."""
    statuses, workflow_status = inspect_external_engine_surface(
        workflow_id=workflow_id,
        summary=summary,
        required_engines=required_engines,
        executables=executables,
    )
    if workflow_status.runnable:
        return workflow_status
    if preserve_missing_error and len(statuses) == 1:
        status = statuses[0]
        if status.support_status == "missing":
            raise EngineUnavailableError(status.issues[0])
    raise EngineWorkflowError(
        f"workflow '{workflow_status.workflow_id}' is blocked by external engine availability or compatibility",
        code="engine_preflight_workflow_blocked",
        details={
            "workflow_id": workflow_status.workflow_id,
            "blocking_engines": workflow_status.blocking_engines,
            "required_engines": workflow_status.required_engines,
            "notes": workflow_status.notes,
        },
    )
