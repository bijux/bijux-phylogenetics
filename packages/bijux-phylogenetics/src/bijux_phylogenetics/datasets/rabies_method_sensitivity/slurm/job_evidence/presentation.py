from __future__ import annotations

from pathlib import Path
import shlex

from bijux_phylogenetics.render.html import write_html_report

from .contracts import RabiesMethodSensitivitySlurmJobEvidenceRow
from .payloads import _relative_bundle_path
from .serialization import _format_float


def _write_job_html_report(
    *,
    path: Path,
    row: RabiesMethodSensitivitySlurmJobEvidenceRow,
    payload: dict[str, object],
    warnings: list[str],
    output_inventory: list[dict[str, object]],
    execution_record_path: Path,
    workflow_manifest_path: Path,
) -> Path:
    commands = payload["commands"]
    sections = [
        (
            "job",
            "\n".join(
                [
                    f"variant_id: {row.variant_id}",
                    f"label: {row.label}",
                    f"status: {row.status}",
                    f"partition_id: {row.partition_id}",
                    f"array_index: {row.array_index}",
                    f"execution_mode: {row.execution_mode}",
                    f"script_path: {row.script_path}",
                ]
            ),
        ),
        (
            "configuration",
            "\n".join(
                [
                    f"alignment_mode: {row.alignment_mode}",
                    f"trimming_mode: {row.trimming_mode}",
                    f"trim_gap_threshold: {_format_float(row.trim_gap_threshold)}",
                    f"selected_model: {row.selected_model}",
                    f"alignment_length: {row.alignment_length}",
                    f"trimmed_alignment_length: {row.trimmed_alignment_length}",
                ]
            ),
        ),
        (
            "commands",
            "\n".join(
                [
                    f"{name}: {shlex.join(list(command))}"
                    for name, command in dict(commands).items()
                ]
            ),
        ),
        (
            "runtime-and-findings",
            "\n".join(
                [
                    f"total_runtime_seconds: {_format_float(row.total_runtime_seconds)}",
                    (
                        "model_selection_runtime_seconds: "
                        f"{_format_float(row.model_selection_runtime_seconds)}"
                    ),
                    (
                        "iqtree_support_runtime_seconds: "
                        f"{_format_float(row.iqtree_support_runtime_seconds)}"
                    ),
                    (
                        "fasttree_runtime_seconds: "
                        f"{_format_float(row.fasttree_runtime_seconds)}"
                    ),
                    f"serious_conflict_count: {row.serious_conflict_count}",
                    f"rooted_engine_rf_distance: {row.rooted_engine_rf_distance}",
                    (
                        "rooted_same_taxa_different_rooting: "
                        f"{str(row.rooted_same_taxa_different_rooting).lower()}"
                    ),
                ]
            ),
        ),
        (
            "warnings",
            "none" if not warnings else "\n".join(warnings),
        ),
        (
            "output-inventory",
            "\n".join(
                f"{item['relative_path']}: {item['size_bytes']} bytes {item['sha256']}"
                for item in output_inventory
            ),
        ),
    ]
    artifact_links = [
        ("task-log", Path("task.log").as_posix(), None),
        (
            "workflow-manifest",
            _relative_bundle_path(path, workflow_manifest_path),
            None,
        ),
        (
            "execution-record",
            _relative_bundle_path(path, execution_record_path),
            None,
        ),
        ("alignment-manifest", "alignment.manifest.json", None),
        ("trimming-manifest", "trimming.manifest.json", None),
        ("inference-comparison-manifest", "inference-comparison.manifest.json", None),
        ("model-selection-manifest", "model-selection.manifest.json", None),
        ("iqtree-support-manifest", "iqtree-support.manifest.json", None),
        ("fasttree-manifest", "fasttree.manifest.json", None),
    ]
    return write_html_report(
        title=f"Rabies Slurm Job Evidence: {row.variant_id}",
        sections=sections,
        out_path=path,
        embedded_json=payload,
        summary_metrics=[
            ("variant", row.variant_id),
            ("status", row.status),
            ("selected model", row.selected_model),
            ("runtime seconds", _format_float(row.total_runtime_seconds)),
            ("output files", row.output_file_count),
            ("warnings", row.warning_count),
            ("serious conflicts", row.serious_conflict_count),
            ("rooted RF", row.rooted_engine_rf_distance),
        ],
        artifact_links=artifact_links,
    )
