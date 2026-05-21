from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path

__all__ = [
    "ArtifactSchemaValidationReport",
    "artifact_schema_names",
    "assert_artifact_schema_valid",
    "validate_artifact_schema",
    "validate_clade_table_schema",
    "validate_comparative_report_manifest_schema",
    "validate_comparative_summary_table_schema",
    "validate_comparative_traits_table_schema",
    "validate_fasta_to_tree_manifest_schema",
    "validate_fasta_to_tree_model_table_schema",
    "validate_fasta_to_tree_support_table_schema",
    "validate_geographic_event_table_schema",
    "validate_host_switch_branch_table_schema",
    "validate_run_manifest_schema",
]

_FASTA_TO_TREE_MODEL_COLUMNS = (
    "workflow",
    "engine_name",
    "sequence_type",
    "selected_model",
    "report_selected_model",
    "artifact_selected_model",
    "model_consistent",
    "alignment_path",
    "trimmed_alignment_path",
    "manifest_path",
)
_FASTA_TO_TREE_SUPPORT_COLUMNS = (
    "node",
    "descendant_taxa",
    "support",
    "support_fraction",
    "is_backbone",
)
_CLADE_TABLE_CORE_COLUMNS = (
    "source_path",
    "tree_index",
    "node_kind",
    "clade_id",
    "node_label",
    "taxon_count",
    "taxa",
    "support",
    "support_fraction",
    "branch_length",
    "root_depth",
    "descendant_tip_depth_min",
    "descendant_tip_depth_max",
    "node_age",
)
_HOST_SWITCH_BRANCH_COLUMNS = (
    "branch_id",
    "parent_node",
    "child_node",
    "child_descendant_taxa",
    "branch_length",
    "parent_most_likely_host",
    "child_most_likely_host",
    "parent_host_set",
    "child_host_set",
    "overlapping_hosts",
    "changed",
    "transition",
    "certainty_class",
    "parent_confidence",
    "child_confidence",
    "transition_allowed",
)
_GEOGRAPHIC_EVENT_COLUMNS = (
    "branch_id",
    "parent_node",
    "child_node",
    "child_descendant_taxa",
    "branch_length",
    "parent_depth",
    "child_depth",
    "midpoint_depth",
    "source_region",
    "target_region",
    "support",
    "strongly_supported",
    "confidence_class",
)
_COMPARATIVE_SUMMARY_COLUMNS = (
    "response",
    "formula",
    "predictor_count",
    "analysis_taxa",
    "excluded_taxa",
    "selected_model",
    "pgls_lambda",
    "pgls_log_likelihood",
    "pgls_r_squared",
    "phylogenetic_signal_k",
    "phylogenetic_signal_lambda",
    "independent_contrast_count",
    "better_model_aicc_delta",
)
_COMPARATIVE_TRAITS_COLUMNS = (
    "taxon",
    "host_group",
    "region_group",
    "region_latitude",
    "region_longitude",
)
_RUN_MANIFEST_KEYS = (
    "arguments",
    "command",
    "dependency_versions",
    "host_platform",
    "input_checksums",
    "input_paths",
    "output_checksums",
    "output_paths",
    "package_version",
    "python_version",
    "timestamp_utc",
)


@dataclass(frozen=True, slots=True)
class ArtifactSchemaValidationReport:
    """Validation result for one stable reviewer-facing artifact schema."""

    schema_name: str
    path: Path
    artifact_format: str
    valid: bool
    expected_fields: tuple[str, ...]
    observed_fields: tuple[str, ...]
    missing_fields: tuple[str, ...]
    unexpected_fields: tuple[str, ...]
    order_matches: bool
    notes: tuple[str, ...]

    @property
    def failure_message(self) -> str:
        details = [f"schema `{self.schema_name}` mismatch for {self.path}"]
        if self.missing_fields:
            details.append(f"missing={','.join(self.missing_fields)}")
        if self.unexpected_fields:
            details.append(f"unexpected={','.join(self.unexpected_fields)}")
        if not self.order_matches:
            details.append("field order changed")
        if self.notes:
            details.append("; ".join(self.notes))
        return (
            ": ".join([details[0], " | ".join(details[1:])])
            if len(details) > 1
            else details[0]
        )


def artifact_schema_names() -> tuple[str, ...]:
    """Return the stable schema profile names supported by this module."""
    return tuple(sorted(_SCHEMA_VALIDATORS))


def assert_artifact_schema_valid(report: ArtifactSchemaValidationReport) -> None:
    """Raise one explicit error when a schema report is invalid."""
    if not report.valid:
        raise ValueError(report.failure_message)


def validate_artifact_schema(
    path: Path,
    schema_name: str,
) -> ArtifactSchemaValidationReport:
    """Validate one artifact path against one named schema profile."""
    try:
        validator = _SCHEMA_VALIDATORS[schema_name]
    except KeyError as error:
        raise ValueError(f"unknown artifact schema: {schema_name}") from error
    return validator(path)


def validate_fasta_to_tree_model_table_schema(
    path: Path,
) -> ArtifactSchemaValidationReport:
    """Validate the canonical `.model.tsv` output from FASTA-to-tree workflows."""
    return _validate_exact_tsv_schema(
        path,
        schema_name="fasta_to_tree_model_tsv",
        expected_fields=_FASTA_TO_TREE_MODEL_COLUMNS,
    )


def validate_fasta_to_tree_support_table_schema(
    path: Path,
) -> ArtifactSchemaValidationReport:
    """Validate the canonical `.support.tsv` output from FASTA-to-tree workflows."""
    return _validate_exact_tsv_schema(
        path,
        schema_name="fasta_to_tree_support_tsv",
        expected_fields=_FASTA_TO_TREE_SUPPORT_COLUMNS,
    )


def validate_clade_table_schema(path: Path) -> ArtifactSchemaValidationReport:
    """Validate one clade table, including stable metadata-column triplets."""
    observed_fields = _read_delimited_header(path)
    core_count = len(_CLADE_TABLE_CORE_COLUMNS)
    core_fields = observed_fields[:core_count]
    metadata_fields = observed_fields[core_count:]
    missing_fields = tuple(
        field for field in _CLADE_TABLE_CORE_COLUMNS if field not in core_fields
    )
    unexpected_fields: list[str] = []
    order_matches = tuple(core_fields) == _CLADE_TABLE_CORE_COLUMNS
    notes: list[str] = []
    if len(metadata_fields) % 3 != 0:
        unexpected_fields.extend(metadata_fields)
        notes.append(
            "metadata columns must appear in values/distinct_values/missing_taxa triplets"
        )
    else:
        seen_metadata_columns: set[str] = set()
        for offset in range(0, len(metadata_fields), 3):
            values_field, distinct_field, missing_field = metadata_fields[
                offset : offset + 3
            ]
            if not values_field.endswith("_values"):
                unexpected_fields.extend((values_field, distinct_field, missing_field))
                continue
            column_name = values_field[: -len("_values")]
            expected_triplet = (
                f"{column_name}_values",
                f"{column_name}_distinct_values",
                f"{column_name}_missing_taxa",
            )
            observed_triplet = (values_field, distinct_field, missing_field)
            if observed_triplet != expected_triplet:
                unexpected_fields.extend(observed_triplet)
                continue
            if column_name in seen_metadata_columns:
                unexpected_fields.extend(observed_triplet)
                notes.append(f"duplicate metadata schema triplet for `{column_name}`")
                continue
            seen_metadata_columns.add(column_name)
        notes.append(f"metadata_column_count={len(seen_metadata_columns)}")
    valid = not missing_fields and not unexpected_fields and order_matches
    return ArtifactSchemaValidationReport(
        schema_name="clade_table_tsv",
        path=path,
        artifact_format="tsv",
        valid=valid,
        expected_fields=_CLADE_TABLE_CORE_COLUMNS,
        observed_fields=observed_fields,
        missing_fields=missing_fields,
        unexpected_fields=tuple(unexpected_fields),
        order_matches=order_matches,
        notes=tuple(notes),
    )


def validate_host_switch_branch_table_schema(
    path: Path,
) -> ArtifactSchemaValidationReport:
    """Validate one host-switch branch ledger."""
    return _validate_exact_tsv_schema(
        path,
        schema_name="host_switch_branch_tsv",
        expected_fields=_HOST_SWITCH_BRANCH_COLUMNS,
    )


def validate_geographic_event_table_schema(
    path: Path,
) -> ArtifactSchemaValidationReport:
    """Validate one geographic event ledger."""
    return _validate_exact_tsv_schema(
        path,
        schema_name="geographic_event_tsv",
        expected_fields=_GEOGRAPHIC_EVENT_COLUMNS,
    )


def validate_comparative_summary_table_schema(
    path: Path,
) -> ArtifactSchemaValidationReport:
    """Validate one comparative summary model-output table."""
    return _validate_exact_tsv_schema(
        path,
        schema_name="comparative_summary_tsv",
        expected_fields=_COMPARATIVE_SUMMARY_COLUMNS,
    )


def validate_comparative_traits_table_schema(
    path: Path,
) -> ArtifactSchemaValidationReport:
    """Validate one derived comparative trait table."""
    return _validate_exact_tsv_schema(
        path,
        schema_name="comparative_traits_tsv",
        expected_fields=_COMPARATIVE_TRAITS_COLUMNS,
    )


def validate_run_manifest_schema(path: Path) -> ArtifactSchemaValidationReport:
    """Validate one deterministic run manifest."""
    return _validate_exact_json_schema(
        path,
        schema_name="run_manifest_json",
        expected_fields=_RUN_MANIFEST_KEYS,
    )


def validate_fasta_to_tree_manifest_schema(
    path: Path,
) -> ArtifactSchemaValidationReport:
    """Validate one FASTA-to-tree workflow manifest."""
    return _validate_exact_json_schema(
        path,
        schema_name="fasta_to_tree_manifest_json",
        expected_fields=(
            "alignment_mode",
            "alignment_workflow",
            "bootstrap_replicates",
            "bootstrap_workflow",
            "commands",
            "config",
            "ended_at_utc",
            "engine_artifact_dir",
            "engine_versions",
            "input_checksums",
            "input_path",
            "input_repair",
            "input_validation",
            "iqtree_seed",
            "iqtree_threads",
            "manifest_path",
            "maximum_likelihood_workflow",
            "method_tier",
            "model_rows",
            "model_selection_workflow",
            "model_validation",
            "notes",
            "out_dir",
            "output_checksums",
            "output_paths",
            "prefix",
            "prepared_input_path",
            "repaired_input_validation",
            "run_manifest_path",
            "runtime_seconds",
            "selected_model",
            "sequence_type",
            "stage_fingerprints",
            "started_at_utc",
            "step_manifests",
            "support_rows",
            "support_summary",
            "trimming_mode",
            "trimming_workflow",
            "trim_gap_threshold",
            "warnings",
            "workflow",
        ),
    )


def validate_comparative_report_manifest_schema(
    path: Path,
) -> ArtifactSchemaValidationReport:
    """Validate one comparative report package manifest."""
    return _validate_exact_json_schema(
        path,
        schema_name="comparative_report_manifest_json",
        expected_fields=(
            "input_checksums",
            "input_paths",
            "limitations",
            "metrics",
            "outputs",
            "report_kind",
            "reviewer_audit_checklist",
            "summary",
        ),
    )


def _validate_exact_tsv_schema(
    path: Path,
    *,
    schema_name: str,
    expected_fields: tuple[str, ...],
) -> ArtifactSchemaValidationReport:
    observed_fields = _read_delimited_header(path)
    return _build_exact_report(
        schema_name=schema_name,
        path=path,
        artifact_format="tsv",
        expected_fields=expected_fields,
        observed_fields=observed_fields,
        require_order=True,
    )


def _validate_exact_json_schema(
    path: Path,
    *,
    schema_name: str,
    expected_fields: tuple[str, ...],
) -> ArtifactSchemaValidationReport:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return ArtifactSchemaValidationReport(
            schema_name=schema_name,
            path=path,
            artifact_format="json",
            valid=False,
            expected_fields=expected_fields,
            observed_fields=(),
            missing_fields=tuple(expected_fields),
            unexpected_fields=(),
            order_matches=False,
            notes=("top-level JSON payload must be an object",),
        )
    observed_fields = tuple(payload.keys())
    return _build_exact_report(
        schema_name=schema_name,
        path=path,
        artifact_format="json",
        expected_fields=expected_fields,
        observed_fields=observed_fields,
        require_order=False,
    )


def _build_exact_report(
    *,
    schema_name: str,
    path: Path,
    artifact_format: str,
    expected_fields: tuple[str, ...],
    observed_fields: tuple[str, ...],
    require_order: bool,
) -> ArtifactSchemaValidationReport:
    missing_fields = tuple(
        field for field in expected_fields if field not in observed_fields
    )
    unexpected_fields = tuple(
        field for field in observed_fields if field not in expected_fields
    )
    order_matches = observed_fields == expected_fields if require_order else True
    valid = not missing_fields and not unexpected_fields and order_matches
    return ArtifactSchemaValidationReport(
        schema_name=schema_name,
        path=path,
        artifact_format=artifact_format,
        valid=valid,
        expected_fields=expected_fields,
        observed_fields=observed_fields,
        missing_fields=missing_fields,
        unexpected_fields=unexpected_fields,
        order_matches=order_matches,
        notes=(),
    )


def _read_delimited_header(path: Path) -> tuple[str, ...]:
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError(f"artifact has no header row: {path}") from error
    return tuple(header)


_SCHEMA_VALIDATORS = {
    "clade_table_tsv": validate_clade_table_schema,
    "comparative_report_manifest_json": validate_comparative_report_manifest_schema,
    "comparative_summary_tsv": validate_comparative_summary_table_schema,
    "comparative_traits_tsv": validate_comparative_traits_table_schema,
    "fasta_to_tree_manifest_json": validate_fasta_to_tree_manifest_schema,
    "fasta_to_tree_model_tsv": validate_fasta_to_tree_model_table_schema,
    "fasta_to_tree_support_tsv": validate_fasta_to_tree_support_table_schema,
    "geographic_event_tsv": validate_geographic_event_table_schema,
    "host_switch_branch_tsv": validate_host_switch_branch_table_schema,
    "run_manifest_json": validate_run_manifest_schema,
}
