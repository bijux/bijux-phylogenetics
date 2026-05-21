from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, cast

import yaml

from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from ..bundles import (
    WorkflowResultBundleExtraInput,
    WorkflowResultBundleReport,
    WorkflowResultBundleValidationReport,
    export_workflow_result_bundle,
    validate_workflow_result_bundle,
)
from ..validation import (
    ExternalEnginePreflightReport,
    WorkflowPreflightStatus,
    inspect_external_engine_preflight,
    require_preflight_workflow,
)
from ..workflows.alignment import (
    list_mafft_alignment_modes,
    list_trimal_trimming_modes,
)
from .fasta_to_tree import FastaToTreeWorkflowReport, run_fasta_to_tree_workflow

__all__ = [
    "PhyloWorkflowConfig",
    "WorkflowConfigRunReport",
    "load_phylo_workflow_config",
    "run_phylo_workflow_config",
]


@dataclass(frozen=True, slots=True)
class PhyloWorkflowConfig:
    """One validated one-command workflow configuration."""

    workflow: str
    config_path: Path
    input_fasta_path: Path
    metadata_path: Path | None
    metadata_taxon_column: str | None
    traits_path: Path | None
    traits_taxon_column: str | None
    out_dir: Path
    bundle_root: Path
    prefix: str | None
    mafft_executable: str
    trimal_executable: str
    iqtree_executable: str
    sequence_type: AlignmentAlphabet | None
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float
    bootstrap_replicates: int
    normalize_identifiers: bool
    remove_invalid_records: bool
    iqtree_seed: int
    iqtree_threads: int
    timeout_seconds: float | None
    resume: bool
    incomplete_run_policy: str
    resolved_payload: dict[str, object]

    def extra_bundle_inputs(self) -> list[WorkflowResultBundleExtraInput]:
        """Return auxiliary files that belong in the exported result bundle."""
        extras = [
            WorkflowResultBundleExtraInput(
                label=f"workflow-config-source{self.config_path.suffix or '.yaml'}",
                source_path=self.config_path,
            )
        ]
        if self.metadata_path is not None:
            extras.append(
                WorkflowResultBundleExtraInput(
                    label=f"metadata-{self.metadata_path.name}",
                    source_path=self.metadata_path,
                )
            )
        if self.traits_path is not None:
            extras.append(
                WorkflowResultBundleExtraInput(
                    label=f"traits-{self.traits_path.name}",
                    source_path=self.traits_path,
                )
            )
        return extras

    def engine_executables(self) -> dict[str, str]:
        """Return the configured external-engine executables keyed by engine id."""
        return {
            "mafft": self.mafft_executable,
            "trimal": self.trimal_executable,
            "iqtree": self.iqtree_executable,
        }


@dataclass(slots=True)
class WorkflowConfigRunReport:
    """One end-to-end config-driven workflow execution result."""

    workflow: str
    config_path: Path
    workflow_config: PhyloWorkflowConfig
    preflight: ExternalEnginePreflightReport
    selected_workflow_status: WorkflowPreflightStatus
    fasta_to_tree_report: FastaToTreeWorkflowReport
    bundle_report: WorkflowResultBundleReport
    bundle_validation: WorkflowResultBundleValidationReport
    output_paths: dict[str, Path]
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def load_phylo_workflow_config(config_path: Path) -> PhyloWorkflowConfig:
    """Load and validate one serious workflow config file."""
    config_path = config_path.resolve()
    payload = _load_yaml_payload(config_path)
    workflow = _required_string(
        payload,
        "workflow",
        config_path=config_path,
        context="config",
    )
    if workflow != "fasta-to-tree":
        raise PhylogeneticsError(
            "workflow config currently supports only the fasta-to-tree workflow",
            code="workflow_config_unsupported_workflow",
            details={
                "workflow": workflow,
                "supported_workflows": ["fasta-to-tree"],
                "config_path": str(config_path),
            },
        )

    inputs = _required_mapping(
        payload,
        "inputs",
        config_path=config_path,
        context="config",
    )
    engines = _optional_mapping(
        payload,
        "engines",
        config_path=config_path,
        context="config",
    )
    alignment = _optional_mapping(
        payload,
        "alignment",
        config_path=config_path,
        context="config",
    )
    trimming = _optional_mapping(
        payload,
        "trimming",
        config_path=config_path,
        context="config",
    )
    inference = _optional_mapping(
        payload,
        "inference",
        config_path=config_path,
        context="config",
    )
    output = _required_mapping(
        payload,
        "output",
        config_path=config_path,
        context="config",
    )
    resources = _optional_mapping(
        payload,
        "resources",
        config_path=config_path,
        context="config",
    )

    input_fasta_path = _resolve_existing_path(
        _required_string(
            inputs,
            "fasta",
            config_path=config_path,
            context="inputs",
        ),
        config_path=config_path,
        field_name="inputs.fasta",
    )
    metadata_path = _resolve_optional_existing_path(
        inputs.get("metadata"),
        config_path=config_path,
        field_name="inputs.metadata",
    )
    traits_path = _resolve_optional_existing_path(
        inputs.get("traits"),
        config_path=config_path,
        field_name="inputs.traits",
    )
    metadata_taxon_column = _optional_string(
        inputs.get("metadata_taxon_column"),
        config_path=config_path,
        field_name="inputs.metadata_taxon_column",
    )
    traits_taxon_column = _optional_string(
        inputs.get("traits_taxon_column"),
        config_path=config_path,
        field_name="inputs.traits_taxon_column",
    )

    out_dir = _resolve_output_path(
        _required_string(
            output,
            "out_dir",
            config_path=config_path,
            context="output",
        ),
        config_path=config_path,
    )
    prefix = _optional_string(
        output.get("prefix"),
        config_path=config_path,
        field_name="output.prefix",
    )
    bundle_root_text = _optional_string(
        output.get("bundle_root"),
        config_path=config_path,
        field_name="output.bundle_root",
    )
    bundle_root = (
        out_dir / f"{input_fasta_path.stem if prefix is None else prefix}.result-bundle"
        if bundle_root_text is None
        else _resolve_output_path(bundle_root_text, config_path=config_path)
    )

    sequence_type = _optional_sequence_type(
        alignment.get("sequence_type"),
        config_path=config_path,
        field_name="alignment.sequence_type",
    )
    alignment_mode = _string_with_default(
        alignment.get("mode"),
        default="auto",
        config_path=config_path,
        field_name="alignment.mode",
    )
    trimming_mode = _string_with_default(
        trimming.get("mode"),
        default="gap-threshold",
        config_path=config_path,
        field_name="trimming.mode",
    )
    trim_gap_threshold = _float_with_default(
        trimming.get("gap_threshold"),
        default=0.1,
        config_path=config_path,
        field_name="trimming.gap_threshold",
    )
    bootstrap_replicates = _int_with_default(
        inference.get("bootstrap_replicates"),
        default=1000,
        config_path=config_path,
        field_name="inference.bootstrap_replicates",
    )
    iqtree_seed = _int_with_default(
        inference.get("seed"),
        default=1,
        config_path=config_path,
        field_name="inference.seed",
    )
    iqtree_threads = _int_with_default(
        inference.get("threads"),
        default=1,
        config_path=config_path,
        field_name="inference.threads",
    )
    normalize_identifiers = _bool_with_default(
        alignment.get("normalize_identifiers"),
        default=False,
        config_path=config_path,
        field_name="alignment.normalize_identifiers",
    )
    remove_invalid_records = _bool_with_default(
        alignment.get("remove_invalid_records"),
        default=False,
        config_path=config_path,
        field_name="alignment.remove_invalid_records",
    )
    timeout_seconds = _optional_float(
        resources.get("timeout_seconds"),
        config_path=config_path,
        field_name="resources.timeout_seconds",
    )
    resume = _bool_with_default(
        resources.get("resume"),
        default=False,
        config_path=config_path,
        field_name="resources.resume",
    )
    incomplete_run_policy = _string_with_default(
        resources.get("incomplete_run_policy"),
        default="reject",
        config_path=config_path,
        field_name="resources.incomplete_run_policy",
    )

    mafft_executable = _string_with_default(
        engines.get("mafft_executable"),
        default="mafft",
        config_path=config_path,
        field_name="engines.mafft_executable",
    )
    trimal_executable = _string_with_default(
        engines.get("trimal_executable"),
        default="trimal",
        config_path=config_path,
        field_name="engines.trimal_executable",
    )
    iqtree_executable = _string_with_default(
        engines.get("iqtree_executable"),
        default="iqtree2",
        config_path=config_path,
        field_name="engines.iqtree_executable",
    )

    available_alignment_modes = set(list_mafft_alignment_modes())
    if alignment_mode not in available_alignment_modes:
        raise _invalid_config_error(
            "alignment.mode must be one supported MAFFT alignment mode",
            config_path=config_path,
            details={
                "field": "alignment.mode",
                "value": alignment_mode,
                "available_modes": sorted(available_alignment_modes),
            },
        )
    available_trimming_modes = set(list_trimal_trimming_modes())
    if trimming_mode not in available_trimming_modes:
        raise _invalid_config_error(
            "trimming.mode must be one supported trimAl trimming mode",
            config_path=config_path,
            details={
                "field": "trimming.mode",
                "value": trimming_mode,
                "available_modes": sorted(available_trimming_modes),
            },
        )
    if incomplete_run_policy not in {"clean", "reject"}:
        raise _invalid_config_error(
            "resources.incomplete_run_policy must be one of: clean, reject",
            config_path=config_path,
            details={"field": "resources.incomplete_run_policy"},
        )
    if not 0.0 <= trim_gap_threshold <= 1.0:
        raise _invalid_config_error(
            "trimming.gap_threshold must be between 0.0 and 1.0",
            config_path=config_path,
            details={
                "field": "trimming.gap_threshold",
                "value": trim_gap_threshold,
            },
        )
    if bootstrap_replicates < 1000:
        raise _invalid_config_error(
            "inference.bootstrap_replicates must be at least 1000 for IQ-TREE ultrafast bootstrap support",
            config_path=config_path,
            details={
                "field": "inference.bootstrap_replicates",
                "value": bootstrap_replicates,
            },
        )
    if iqtree_threads < 1:
        raise _invalid_config_error(
            "inference.threads must be at least 1",
            config_path=config_path,
            details={"field": "inference.threads", "value": iqtree_threads},
        )
    if timeout_seconds is not None and timeout_seconds <= 0.0:
        raise _invalid_config_error(
            "resources.timeout_seconds must be greater than zero when configured",
            config_path=config_path,
            details={
                "field": "resources.timeout_seconds",
                "value": timeout_seconds,
            },
        )

    resolved_payload = {
        "workflow": workflow,
        "source_config_path": str(config_path),
        "inputs": {
            "fasta": str(input_fasta_path),
            "metadata": None if metadata_path is None else str(metadata_path),
            "metadata_taxon_column": metadata_taxon_column,
            "traits": None if traits_path is None else str(traits_path),
            "traits_taxon_column": traits_taxon_column,
        },
        "engines": {
            "mafft_executable": mafft_executable,
            "trimal_executable": trimal_executable,
            "iqtree_executable": iqtree_executable,
        },
        "alignment": {
            "sequence_type": sequence_type,
            "mode": alignment_mode,
            "normalize_identifiers": normalize_identifiers,
            "remove_invalid_records": remove_invalid_records,
        },
        "trimming": {
            "mode": trimming_mode,
            "gap_threshold": trim_gap_threshold,
        },
        "inference": {
            "bootstrap_replicates": bootstrap_replicates,
            "seed": iqtree_seed,
            "threads": iqtree_threads,
        },
        "output": {
            "out_dir": str(out_dir),
            "prefix": prefix,
            "bundle_root": str(bundle_root),
        },
        "resources": {
            "timeout_seconds": timeout_seconds,
            "resume": resume,
            "incomplete_run_policy": incomplete_run_policy,
        },
    }

    return PhyloWorkflowConfig(
        workflow=workflow,
        config_path=config_path,
        input_fasta_path=input_fasta_path,
        metadata_path=metadata_path,
        metadata_taxon_column=metadata_taxon_column,
        traits_path=traits_path,
        traits_taxon_column=traits_taxon_column,
        out_dir=out_dir,
        bundle_root=bundle_root,
        prefix=prefix,
        mafft_executable=mafft_executable,
        trimal_executable=trimal_executable,
        iqtree_executable=iqtree_executable,
        sequence_type=sequence_type,
        alignment_mode=alignment_mode,
        trimming_mode=trimming_mode,
        trim_gap_threshold=trim_gap_threshold,
        bootstrap_replicates=bootstrap_replicates,
        normalize_identifiers=normalize_identifiers,
        remove_invalid_records=remove_invalid_records,
        iqtree_seed=iqtree_seed,
        iqtree_threads=iqtree_threads,
        timeout_seconds=timeout_seconds,
        resume=resume,
        incomplete_run_policy=incomplete_run_policy,
        resolved_payload=resolved_payload,
    )


def run_phylo_workflow_config(config_path: Path) -> WorkflowConfigRunReport:
    """Run one serious workflow from one validated config file."""
    workflow_config = load_phylo_workflow_config(config_path)
    preflight = inspect_external_engine_preflight(
        executables=workflow_config.engine_executables(),
        selected_workflow=workflow_config.workflow,
    )
    selected_workflow_status = require_preflight_workflow(
        preflight,
        workflow_id=workflow_config.workflow,
    )
    workflow_report = run_fasta_to_tree_workflow(
        workflow_config.input_fasta_path,
        out_dir=workflow_config.out_dir,
        prefix=workflow_config.prefix,
        sequence_type=workflow_config.sequence_type,
        mafft_executable=workflow_config.mafft_executable,
        alignment_mode=workflow_config.alignment_mode,
        trimal_executable=workflow_config.trimal_executable,
        trimming_mode=workflow_config.trimming_mode,
        iqtree_executable=workflow_config.iqtree_executable,
        iqtree_seed=workflow_config.iqtree_seed,
        iqtree_threads=workflow_config.iqtree_threads,
        trim_gap_threshold=workflow_config.trim_gap_threshold,
        bootstrap_replicates=workflow_config.bootstrap_replicates,
        normalize_identifiers=workflow_config.normalize_identifiers,
        remove_invalid_records=workflow_config.remove_invalid_records,
        resume=workflow_config.resume,
        timeout_seconds=workflow_config.timeout_seconds,
        incomplete_run_policy=workflow_config.incomplete_run_policy,
    )
    bundle_notes = [
        "workflow config source was copied into the result bundle alongside the resolved workflow config",
    ]
    if (
        workflow_config.metadata_path is not None
        or workflow_config.traits_path is not None
    ):
        bundle_notes.append(
            "metadata and traits files were bundled for reviewer context and downstream comparative work; the current fasta-to-tree execution path does not consume them during tree building"
        )
    bundle_report = export_workflow_result_bundle(
        workflow_report.manifest_path,
        bundle_root=workflow_config.bundle_root,
        config_payload=workflow_config.resolved_payload,
        extra_inputs=workflow_config.extra_bundle_inputs(),
        extra_notes=bundle_notes,
    )
    bundle_validation = validate_workflow_result_bundle(workflow_config.bundle_root)
    if not bundle_validation.valid:
        raise PhylogeneticsError(
            "workflow config run produced a result bundle that failed validation",
            code="workflow_config_bundle_validation_failed",
            details={
                "bundle_root": str(workflow_config.bundle_root),
                "issues": [
                    {
                        "kind": issue.kind,
                        "label": issue.label,
                        "detail": issue.detail,
                        "relative_path": (
                            None
                            if issue.relative_path is None
                            else issue.relative_path.as_posix()
                        ),
                    }
                    for issue in bundle_validation.issues
                ],
            },
        )

    output_paths = dict(workflow_report.output_paths)
    output_paths["bundle_root"] = workflow_config.bundle_root
    output_paths["bundle_manifest"] = bundle_report.bundle_manifest_path
    output_paths["bundle_report"] = bundle_report.report_path
    notes = [
        "config-driven phylo run executed through the governed fasta-to-tree workflow",
        f"selected workflow readiness: {selected_workflow_status.readiness_status}",
        f"bundle validation passed for {workflow_config.bundle_root}",
    ]
    warnings = list(dict.fromkeys(workflow_report.warnings))
    return WorkflowConfigRunReport(
        workflow=workflow_config.workflow,
        config_path=config_path.resolve(),
        workflow_config=workflow_config,
        preflight=preflight,
        selected_workflow_status=selected_workflow_status,
        fasta_to_tree_report=workflow_report,
        bundle_report=bundle_report,
        bundle_validation=bundle_validation,
        output_paths=output_paths,
        warnings=warnings,
        notes=notes,
    )


def _load_yaml_payload(config_path: Path) -> dict[str, Any]:
    if not config_path.exists():
        raise _invalid_config_error(
            f"workflow config file does not exist: {config_path}",
            config_path=config_path,
            code="workflow_config_missing_file",
        )
    try:
        payload = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as error:
        details: dict[str, object] = {}
        mark = getattr(error, "problem_mark", None)
        if mark is not None:
            details["line"] = mark.line + 1
            details["column"] = mark.column + 1
        raise _invalid_config_error(
            f"workflow config could not be parsed: {error}",
            config_path=config_path,
            code="workflow_config_parse_error",
            details=details,
        ) from error
    if not isinstance(payload, dict):
        raise _invalid_config_error(
            "workflow config must contain one top-level mapping",
            config_path=config_path,
        )
    return cast(dict[str, Any], payload)


def _required_mapping(
    payload: dict[str, Any],
    key: str,
    *,
    config_path: Path,
    context: str,
) -> dict[str, Any]:
    value = payload.get(key)
    if value is None:
        raise _invalid_config_error(
            f"{context} is missing required mapping '{key}'",
            config_path=config_path,
            details={"field": f"{context}.{key}"},
        )
    if not isinstance(value, dict):
        raise _invalid_config_error(
            f"{context}.{key} must be a mapping",
            config_path=config_path,
            details={"field": f"{context}.{key}"},
        )
    return cast(dict[str, Any], value)


def _optional_mapping(
    payload: dict[str, Any],
    key: str,
    *,
    config_path: Path,
    context: str,
) -> dict[str, Any]:
    value = payload.get(key)
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise _invalid_config_error(
            f"{context}.{key} must be a mapping",
            config_path=config_path,
            details={"field": f"{context}.{key}"},
        )
    return cast(dict[str, Any], value)


def _required_string(
    payload: dict[str, Any],
    key: str,
    *,
    config_path: Path,
    context: str,
) -> str:
    value = payload.get(key)
    if isinstance(value, str) and value.strip():
        return value
    raise _invalid_config_error(
        f"{context}.{key} must be a non-empty string",
        config_path=config_path,
        details={"field": f"{context}.{key}"},
    )


def _optional_string(
    value: object,
    *,
    config_path: Path,
    field_name: str,
) -> str | None:
    if value is None:
        return None
    if isinstance(value, str) and value.strip():
        return value
    raise _invalid_config_error(
        f"{field_name} must be a non-empty string",
        config_path=config_path,
        details={"field": field_name},
    )


def _string_with_default(
    value: object,
    *,
    default: str,
    config_path: Path,
    field_name: str,
) -> str:
    resolved = _optional_string(
        value,
        config_path=config_path,
        field_name=field_name,
    )
    return default if resolved is None else resolved


def _bool_with_default(
    value: object,
    *,
    default: bool,
    config_path: Path,
    field_name: str,
) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    raise _invalid_config_error(
        f"{field_name} must be a boolean",
        config_path=config_path,
        details={"field": field_name},
    )


def _int_with_default(
    value: object,
    *,
    default: int,
    config_path: Path,
    field_name: str,
) -> int:
    if value is None:
        return default
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    raise _invalid_config_error(
        f"{field_name} must be an integer",
        config_path=config_path,
        details={"field": field_name},
    )


def _float_with_default(
    value: object,
    *,
    default: float,
    config_path: Path,
    field_name: str,
) -> float:
    if value is None:
        return default
    return _coerce_float(
        value,
        config_path=config_path,
        field_name=field_name,
    )


def _optional_float(
    value: object,
    *,
    config_path: Path,
    field_name: str,
) -> float | None:
    if value is None:
        return None
    return _coerce_float(
        value,
        config_path=config_path,
        field_name=field_name,
    )


def _coerce_float(
    value: object,
    *,
    config_path: Path,
    field_name: str,
) -> float:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)
    raise _invalid_config_error(
        f"{field_name} must be a number",
        config_path=config_path,
        details={"field": field_name},
    )


def _optional_sequence_type(
    value: object,
    *,
    config_path: Path,
    field_name: str,
) -> AlignmentAlphabet | None:
    if value is None:
        return None
    if isinstance(value, str) and value in {"dna", "rna", "protein", "unknown"}:
        return value
    raise _invalid_config_error(
        f"{field_name} must be one of: dna, protein, rna, unknown",
        config_path=config_path,
        details={"field": field_name},
    )


def _resolve_existing_path(
    path_text: str,
    *,
    config_path: Path,
    field_name: str,
) -> Path:
    resolved = _resolve_relative_path(path_text, config_path=config_path)
    if not resolved.exists():
        raise _invalid_config_error(
            f"{field_name} does not exist: {resolved}",
            config_path=config_path,
            code="workflow_config_missing_path",
            details={"field": field_name, "path": str(resolved)},
        )
    return resolved


def _resolve_optional_existing_path(
    value: object,
    *,
    config_path: Path,
    field_name: str,
) -> Path | None:
    path_text = _optional_string(
        value,
        config_path=config_path,
        field_name=field_name,
    )
    if path_text is None:
        return None
    return _resolve_existing_path(
        path_text,
        config_path=config_path,
        field_name=field_name,
    )


def _resolve_output_path(path_text: str, *, config_path: Path) -> Path:
    return _resolve_relative_path(path_text, config_path=config_path)


def _resolve_relative_path(path_text: str, *, config_path: Path) -> Path:
    candidate = Path(path_text)
    if candidate.is_absolute():
        return candidate
    return (config_path.parent / candidate).resolve()


def _invalid_config_error(
    message: str,
    *,
    config_path: Path,
    code: str = "workflow_config_invalid",
    details: dict[str, object] | None = None,
) -> PhylogeneticsError:
    payload = {"config_path": str(config_path)}
    if details is not None:
        payload.update(details)
    return PhylogeneticsError(message, code=code, details=payload)
