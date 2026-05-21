from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, is_dataclass
import json
from pathlib import Path
from typing import Any

from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport
from bijux_phylogenetics.ancestral.presentation.report_rendering import (
    write_ancestral_state_table,
)
from bijux_phylogenetics.comparative.pgls import PGLSResult
from bijux_phylogenetics.compare.topology import TreeComparisonReport
from bijux_phylogenetics.engines.artifacts.bootstrap import (
    build_bootstrap_support_rows,
    write_bootstrap_support_table,
)
from bijux_phylogenetics.engines.inference import (
    FastaToTreeWorkflowReport,
    WorkflowConfigRunReport,
)
from bijux_phylogenetics.engines.workflows.models import EngineWorkflowReport
from bijux_phylogenetics.phylo.alignment import FastaInputValidationReport
from bijux_phylogenetics.reports.service import ReportBuildResult

__all__ = [
    "AlignmentWorkflowResult",
    "AncestralReconstructionWorkflowResult",
    "ComparativeModelWorkflowResult",
    "ConfiguredPhyloWorkflowResult",
    "FastaValidationResult",
    "InferenceWorkflowResult",
    "ReportWorkflowResult",
    "SequenceToTreeWorkflowResult",
    "SupportWorkflowResult",
    "TreeComparisonWorkflowResult",
    "TrimmingWorkflowResult",
]


def _serialize_value(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value):
        return {key: _serialize_value(item) for key, item in asdict(value).items()}
    if isinstance(value, dict):
        return {str(key): _serialize_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serialize_value(item) for item in value]
    return value


def _write_json(path: Path, payload: object) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_serialize_value(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def _stringify_cell(value: object) -> str:
    if value is None:
        return ""
    serialized = _serialize_value(value)
    if isinstance(serialized, bool):
        return "true" if serialized else "false"
    if isinstance(serialized, (dict, list)):
        return json.dumps(serialized, sort_keys=True)
    return str(serialized)


def _write_table(path: Path, rows: list[dict[str, object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return path
    columns = list(rows[0].keys())
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {column: _stringify_cell(row.get(column)) for column in columns}
            )
    return path


@dataclass(slots=True)
class _SerializableWorkflowResult:
    def __getattr__(self, name: str) -> Any:
        report = object.__getattribute__(self, "report")
        return getattr(report, name)

    def to_json_payload(self) -> dict[str, object]:
        return {
            "result_type": self.__class__.__name__,
            "report": _serialize_value(object.__getattribute__(self, "report")),
        }

    def write_json(self, path: Path) -> Path:
        return _write_json(path, self.to_json_payload())


@dataclass(slots=True)
class FastaValidationResult(_SerializableWorkflowResult):
    report: FastaInputValidationReport

    def write_tsv(self, path: Path) -> Path:
        return _write_table(
            path,
            [
                {
                    "path": self.report.path,
                    "sequence_count": self.report.summary.sequence_count,
                    "selected_type": self.report.sequence_type_report.selected_type,
                    "duplicate_identifier_count": len(
                        self.report.duplicate_identifiers
                    ),
                    "illegal_character_count": len(self.report.illegal_characters),
                    "empty_sequence_count": len(self.report.empty_sequences),
                    "length_outlier_count": len(self.report.length_outliers),
                    "warning_count": len(self.report.warnings),
                    "warnings": self.report.warnings,
                }
            ],
        )


@dataclass(slots=True)
class AlignmentWorkflowResult(_SerializableWorkflowResult):
    report: EngineWorkflowReport

    def write_tsv(self, path: Path) -> Path:
        return _write_table(
            path,
            [
                {
                    "workflow": self.report.workflow,
                    "engine_name": self.report.engine_name,
                    "manifest_path": self.report.manifest_path,
                    "alignment_path": self.report.output_paths.get("alignment"),
                    "version": self.report.run.version,
                    "command": self.report.run.command,
                    "runtime_seconds": self.report.run.runtime_seconds,
                }
            ],
        )


@dataclass(slots=True)
class TrimmingWorkflowResult(_SerializableWorkflowResult):
    report: EngineWorkflowReport

    def write_tsv(self, path: Path) -> Path:
        summary = self.report.trimming_summary
        return _write_table(
            path,
            [
                {
                    "workflow": self.report.workflow,
                    "engine_name": self.report.engine_name,
                    "manifest_path": self.report.manifest_path,
                    "trimmed_alignment_path": self.report.output_paths.get(
                        "trimmed_alignment"
                    ),
                    "mode": None if summary is None else summary.mode,
                    "gap_threshold": None if summary is None else summary.gap_threshold,
                    "input_alignment_length": (
                        None if summary is None else summary.input_alignment_length
                    ),
                    "trimmed_alignment_length": (
                        None if summary is None else summary.trimmed_alignment_length
                    ),
                    "removed_site_count": (
                        None if summary is None else summary.removed_site_count
                    ),
                    "runtime_seconds": self.report.run.runtime_seconds,
                }
            ],
        )


@dataclass(slots=True)
class InferenceWorkflowResult(_SerializableWorkflowResult):
    report: EngineWorkflowReport

    def write_tsv(self, path: Path) -> Path:
        return _write_table(
            path,
            [
                {
                    "workflow": self.report.workflow,
                    "engine_name": self.report.engine_name,
                    "manifest_path": self.report.manifest_path,
                    "tree_path": self.report.output_paths.get("tree"),
                    "selected_model": self.report.selected_model,
                    "log_likelihood": self.report.log_likelihood,
                    "runtime_seconds": self.report.run.runtime_seconds,
                }
            ],
        )


@dataclass(slots=True)
class SupportWorkflowResult(_SerializableWorkflowResult):
    report: EngineWorkflowReport

    def write_tsv(self, path: Path) -> Path:
        summary = self.report.bootstrap_support_summary
        if summary is None:
            return _write_table(
                path,
                [
                    {
                        "workflow": self.report.workflow,
                        "engine_name": self.report.engine_name,
                        "support_tree_path": self.report.output_paths.get(
                            "support_tree"
                        ),
                        "supported_node_count": 0,
                    }
                ],
            )
        return write_bootstrap_support_table(
            path, build_bootstrap_support_rows(summary)
        )


@dataclass(slots=True)
class SequenceToTreeWorkflowResult(_SerializableWorkflowResult):
    report: FastaToTreeWorkflowReport

    def write_tsv(self, path: Path) -> Path:
        return _write_table(
            path,
            [
                {
                    "workflow": self.report.workflow,
                    "input_path": self.report.input_path,
                    "sequence_type": self.report.sequence_type,
                    "selected_model": self.report.selected_model,
                    "alignment_path": self.report.output_paths.get("alignment"),
                    "trimmed_alignment_path": self.report.output_paths.get(
                        "trimmed_alignment"
                    ),
                    "tree_path": self.report.output_paths.get("tree"),
                    "support_table_path": self.report.output_paths.get("support_table"),
                    "runtime_seconds": self.report.runtime_seconds,
                }
            ],
        )


@dataclass(slots=True)
class TreeComparisonWorkflowResult(_SerializableWorkflowResult):
    report: TreeComparisonReport

    def write_tsv(self, path: Path) -> Path:
        return _write_table(
            path,
            [
                {
                    "left_path": self.report.left_path,
                    "right_path": self.report.right_path,
                    "shared_taxa_count": len(self.report.shared_taxa),
                    "robinson_foulds_distance": self.report.robinson_foulds_distance,
                    "normalized_robinson_foulds": (
                        self.report.normalized_robinson_foulds
                    ),
                    "rooted_robinson_foulds_distance": (
                        self.report.rooted_robinson_foulds_distance
                    ),
                    "unrooted_robinson_foulds_distance": (
                        self.report.unrooted_robinson_foulds_distance
                    ),
                    "topology_equal": self.report.topology_equal,
                    "same_unrooted_topology": self.report.same_unrooted_topology,
                }
            ],
        )


@dataclass(slots=True)
class ComparativeModelWorkflowResult(_SerializableWorkflowResult):
    report: PGLSResult

    def write_tsv(self, path: Path) -> Path:
        return _write_table(
            path,
            [
                {
                    "coefficient": coefficient.name,
                    "estimate": coefficient.estimate,
                    "standard_error": coefficient.standard_error,
                    "test_statistic": coefficient.test_statistic,
                    "p_value": coefficient.p_value,
                    "lower_95_confidence_interval": (
                        coefficient.lower_95_confidence_interval
                    ),
                    "upper_95_confidence_interval": (
                        coefficient.upper_95_confidence_interval
                    ),
                    "lambda_value": self.report.lambda_value,
                    "log_likelihood": self.report.log_likelihood,
                    "aic": self.report.aic,
                    "taxon_count": self.report.taxon_count,
                }
                for coefficient in self.report.coefficients
            ],
        )


@dataclass(slots=True)
class AncestralReconstructionWorkflowResult(_SerializableWorkflowResult):
    report: DiscreteAncestralReport

    def write_tsv(self, path: Path) -> Path:
        return write_ancestral_state_table(path, self.report)


@dataclass(slots=True)
class ReportWorkflowResult(_SerializableWorkflowResult):
    report: ReportBuildResult

    def write_tsv(self, path: Path) -> Path:
        return _write_table(
            path,
            [
                {
                    "report_kind": self.report.report_kind,
                    "title": self.report.title,
                    "output_path": self.report.output_path,
                    "machine_manifest_path": self.report.machine_manifest_path,
                    "validation_rooted": self.report.validation.rooted,
                    "inspection_binary": self.report.inspection.is_binary,
                    "validation_complete_branch_lengths": (
                        self.report.validation.has_complete_branch_lengths
                    ),
                    "alignment_taxa": (
                        None
                        if self.report.alignment is None
                        else self.report.alignment.sequence_count
                    ),
                    "trait_missing_value_count": (
                        None
                        if self.report.trait_missing_values is None
                        else len(self.report.trait_missing_values.missing_values)
                    ),
                }
            ],
        )


@dataclass(slots=True)
class ConfiguredPhyloWorkflowResult(_SerializableWorkflowResult):
    report: WorkflowConfigRunReport

    def write_tsv(self, path: Path) -> Path:
        return _write_table(
            path,
            [
                {
                    "workflow": self.report.workflow,
                    "config_path": self.report.config_path,
                    "selected_workflow_status": (
                        self.report.selected_workflow_status.readiness_status
                    ),
                    "bundle_manifest_path": (
                        self.report.bundle_report.bundle_manifest_path
                    ),
                    "bundle_valid": self.report.bundle_validation.valid,
                    "tree_path": self.report.fasta_to_tree_report.output_paths.get(
                        "tree"
                    ),
                    "selected_model": self.report.fasta_to_tree_report.selected_model,
                    "warning_count": len(self.report.warnings),
                    "note_count": len(self.report.notes),
                }
            ],
        )
