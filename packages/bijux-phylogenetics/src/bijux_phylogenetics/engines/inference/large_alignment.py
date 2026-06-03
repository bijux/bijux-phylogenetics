from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess  # nosec B404
import time
import tracemalloc

from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet
from bijux_phylogenetics.runtime.errors import (
    EngineWorkflowError,
    InvalidAlignmentError,
)

from ..artifacts.fasttree import (
    build_fasttree_low_support_rows,
    build_fasttree_support_histogram_rows,
    build_fasttree_support_rows,
    write_fasttree_support_histogram,
    write_fasttree_support_table,
)
from ..artifacts.support import FastTreeSupportNode, FastTreeSupportSummaryReport
from ..common import (
    EngineIncompleteRunRecord,
    build_file_checksums,
    clear_incomplete_engine_run,
    load_engine_manifest,
    read_engine_version,
    resolve_engine_executable,
    utc_now_text,
    validate_timeout_seconds,
    write_engine_manifest,
    write_incomplete_engine_run,
)
from ..validation import (
    summarize_fasttree_support_distribution,
)
from ..validation.preflight import require_external_engine_surface
from ..workflows.state import (
    _resolve_incomplete_workflow_state,
    _validate_support_value_count,
    _validate_tree_output,
)

__all__ = [
    "LargeAlignmentInputSummary",
    "LargeAlignmentResourceRow",
    "LargeAlignmentInferenceWorkflowReport",
    "run_large_alignment_inference",
    "write_large_alignment_inference_log",
    "write_large_alignment_resource_table",
]

_DNA_CHARS = set("ACGTNRYWSKMBDHV?-")
_RNA_CHARS = set("ACGUNRYWSKMBDHV?-")
_SLEEP_INTERVAL_SECONDS = 0.05


@dataclass(frozen=True, slots=True)
class LargeAlignmentInputSummary:
    """One streamed summary of a large aligned FASTA input."""

    input_path: Path
    input_bytes: int
    sequence_count: int
    alignment_length: int
    min_sequence_length: int
    max_sequence_length: int
    total_site_cells: int
    inferred_sequence_type: AlignmentAlphabet


@dataclass(frozen=True, slots=True)
class LargeAlignmentResourceRow:
    """One stage-level resource observation for the large-alignment workflow."""

    stage: str
    elapsed_seconds: float
    peak_memory_bytes: int | None
    memory_observation_kind: str


@dataclass(slots=True)
class LargeAlignmentInferenceWorkflowReport:
    """One governed large-alignment inference run over an aligned FASTA matrix."""

    workflow: str
    input_path: Path
    out_dir: Path
    prefix: str
    sequence_type: AlignmentAlphabet
    engine_name: str
    approximate_method: bool
    timeout_seconds: float | None
    resumed: bool
    started_at_utc: str
    ended_at_utc: str
    runtime_seconds: float
    manifest_path: Path
    output_paths: dict[str, Path]
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    config: dict[str, object]
    input_summary: LargeAlignmentInputSummary
    resource_rows: list[LargeAlignmentResourceRow]
    engine_version_text: str
    command: list[str]
    support_summary: FastTreeSupportSummaryReport
    warnings: list[str]
    notes: list[str]


def _infer_sequence_type(observed_letters: set[str]) -> AlignmentAlphabet:
    if not observed_letters:
        return "unknown"
    if observed_letters <= _DNA_CHARS:
        return "dna"
    if observed_letters <= _RNA_CHARS:
        return "rna"
    return "protein"


def _scan_alignment_streaming(
    path: Path,
) -> tuple[LargeAlignmentInputSummary, LargeAlignmentResourceRow]:
    started = time.perf_counter()
    tracemalloc.start()
    sequence_count = 0
    min_sequence_length: int | None = None
    max_sequence_length = 0
    current_identifier: str | None = None
    current_length = 0
    identifiers: set[str] = set()
    observed_letters: set[str] = set()

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_identifier is not None:
                    if current_length == 0:
                        raise InvalidAlignmentError(
                            f"alignment contains an empty sequence: {current_identifier}"
                        )
                    sequence_count += 1
                    min_sequence_length = (
                        current_length
                        if min_sequence_length is None
                        else min(min_sequence_length, current_length)
                    )
                    max_sequence_length = max(max_sequence_length, current_length)
                current_identifier = line[1:].strip()
                if not current_identifier:
                    raise InvalidAlignmentError(
                        f"alignment contains an empty FASTA identifier: {path}"
                    )
                if current_identifier in identifiers:
                    raise InvalidAlignmentError(
                        f"alignment contains duplicate sequence ids: {current_identifier}"
                    )
                identifiers.add(current_identifier)
                current_length = 0
                continue
            if current_identifier is None:
                raise InvalidAlignmentError(
                    f"alignment sequence appears before any FASTA header in {path}"
                )
            current_length += len(line)
            observed_letters.update(character.upper() for character in line)

    if current_identifier is not None:
        if current_length == 0:
            raise InvalidAlignmentError(
                f"alignment contains an empty sequence: {current_identifier}"
            )
        sequence_count += 1
        min_sequence_length = (
            current_length
            if min_sequence_length is None
            else min(min_sequence_length, current_length)
        )
        max_sequence_length = max(max_sequence_length, current_length)

    _, peak_memory_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed_seconds = time.perf_counter() - started
    if sequence_count == 0 or min_sequence_length is None:
        raise InvalidAlignmentError(f"alignment contains no FASTA records: {path}")
    if min_sequence_length != max_sequence_length:
        raise InvalidAlignmentError(
            "large-alignment inference requires a rectangular aligned FASTA matrix"
        )

    summary = LargeAlignmentInputSummary(
        input_path=path,
        input_bytes=path.stat().st_size,
        sequence_count=sequence_count,
        alignment_length=max_sequence_length,
        min_sequence_length=min_sequence_length,
        max_sequence_length=max_sequence_length,
        total_site_cells=sequence_count * max_sequence_length,
        inferred_sequence_type=_infer_sequence_type(observed_letters),
    )
    return summary, LargeAlignmentResourceRow(
        stage="preflight-scan",
        elapsed_seconds=elapsed_seconds,
        peak_memory_bytes=peak_memory_bytes,
        memory_observation_kind="python-tracemalloc",
    )


def _fasttree_args(path: Path, sequence_type: AlignmentAlphabet) -> list[str]:
    if sequence_type in {"dna", "rna"}:
        return ["-gtr", "-nt", str(path)]
    if sequence_type == "protein":
        return ["-lg", str(path)]
    return [str(path)]


def _resource_table_path(root: Path) -> Path:
    return root.with_suffix(".resources.tsv")


def _stderr_log_path(root: Path) -> Path:
    return root.with_suffix(".stderr.log")


def _plain_log_path(root: Path) -> Path:
    return root.with_suffix(".log")


def _tree_path(root: Path) -> Path:
    return root.with_suffix(".tree")


def _manifest_path(root: Path) -> Path:
    return root.with_suffix(".manifest.json")


def _sample_process_peak_rss_bytes(pid: int) -> int | None:
    ps_executable = shutil.which("ps")
    if ps_executable is None:
        return None
    result = subprocess.run(  # nosec B603
        [ps_executable, "-o", "rss=", "-p", str(pid)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        return None
    text = result.stdout.strip()
    if not text:
        return None
    try:
        return int(text) * 1024
    except ValueError:
        return None


def _execute_large_fasttree_command(
    *,
    executable: str,
    command_args: list[str],
    out_path: Path,
    stderr_path: Path,
    timeout_seconds: float | None,
) -> tuple[list[str], LargeAlignmentResourceRow]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    command = [executable, *command_args]
    started = time.perf_counter()
    peak_memory_bytes: int | None = None
    with (
        out_path.open("w", encoding="utf-8") as stdout_handle,
        stderr_path.open("w", encoding="utf-8") as stderr_handle,
        subprocess.Popen(  # nosec B603
            command,
            cwd=out_path.parent,
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
        ) as process,
    ):
        while True:
            current_peak = _sample_process_peak_rss_bytes(process.pid)
            if current_peak is not None:
                peak_memory_bytes = (
                    current_peak
                    if peak_memory_bytes is None
                    else max(peak_memory_bytes, current_peak)
                )
            exit_code = process.poll()
            elapsed_seconds = time.perf_counter() - started
            if exit_code is not None:
                break
            if timeout_seconds is not None and elapsed_seconds > timeout_seconds:
                process.kill()
                process.wait()
                raise EngineWorkflowError(
                    "FastTree large-alignment inference timed out after "
                    f"{timeout_seconds:.3f} seconds; stderr log: {stderr_path}"
                )
            time.sleep(_SLEEP_INTERVAL_SECONDS)
    final_peak = _sample_process_peak_rss_bytes(process.pid)
    if final_peak is not None:
        peak_memory_bytes = (
            final_peak
            if peak_memory_bytes is None
            else max(peak_memory_bytes, final_peak)
        )
    stderr_text = (
        stderr_path.read_text(encoding="utf-8") if stderr_path.exists() else ""
    )
    warning_lines = [
        line.strip()
        for line in stderr_text.splitlines()
        if line.strip() and "warn" in line.lower()
    ]
    if exit_code != 0:
        raise EngineWorkflowError(
            f"FastTree large-alignment inference failed with exit code {exit_code}; stderr log: {stderr_path}"
        )
    if not out_path.exists() or not out_path.read_text(encoding="utf-8").strip():
        raise EngineWorkflowError(
            f"FastTree large-alignment inference did not produce a tree: {out_path}"
        )
    return warning_lines, LargeAlignmentResourceRow(
        stage="fasttree-inference",
        elapsed_seconds=time.perf_counter() - started,
        peak_memory_bytes=peak_memory_bytes,
        memory_observation_kind="sampled-process-rss",
    )


def write_large_alignment_resource_table(
    path: Path, rows: list[LargeAlignmentResourceRow]
) -> Path:
    """Write one stage-level resource ledger for large-alignment inference."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "stage",
                "elapsed_seconds",
                "peak_memory_bytes",
                "memory_observation_kind",
            ]
        )
    ]
    for row in rows:
        lines.append(
            "\t".join(
                [
                    row.stage,
                    format(row.elapsed_seconds, ".12g"),
                    "" if row.peak_memory_bytes is None else str(row.peak_memory_bytes),
                    row.memory_observation_kind,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_large_alignment_inference_log(
    path: Path,
    report: LargeAlignmentInferenceWorkflowReport,
) -> Path:
    """Write a reviewer-facing plain-text log for large-alignment inference."""
    lines = [
        "workflow: large-alignment-inference",
        f"input_path: {report.input_path}",
        f"sequence_type: {report.sequence_type}",
        f"engine_name: {report.engine_name}",
        f"approximate_method: {str(report.approximate_method).lower()}",
        f"resumed: {str(report.resumed).lower()}",
        f"timeout_seconds: {'' if report.timeout_seconds is None else format(report.timeout_seconds, '.12g')}",
        "",
        "[input-summary]",
        f"sequence_count: {report.input_summary.sequence_count}",
        f"alignment_length: {report.input_summary.alignment_length}",
        f"total_site_cells: {report.input_summary.total_site_cells}",
        f"input_bytes: {report.input_summary.input_bytes}",
        f"inferred_sequence_type: {report.input_summary.inferred_sequence_type}",
        "",
        "[engine-run]",
        f"command: {' '.join(report.command)}",
        f"version: {report.engine_version_text}",
    ]
    for output_label, output_path in report.output_paths.items():
        lines.append(f"output.{output_label}: {output_path}")
    lines.extend(["", "[resources]"])
    for row in report.resource_rows:
        lines.append(
            f"- {row.stage}: elapsed_seconds={format(row.elapsed_seconds, '.12g')}, "
            f"peak_memory_bytes={'' if row.peak_memory_bytes is None else row.peak_memory_bytes}, "
            f"memory_observation_kind={row.memory_observation_kind}"
        )
    if report.warnings:
        lines.extend(["", "[warnings]"])
        lines.extend(f"- {warning}" for warning in report.warnings)
    if report.notes:
        lines.extend(["", "[notes]"])
        lines.extend(f"- {note}" for note in report.notes)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return path


def _restore_support_summary(
    payload: dict[str, object],
) -> FastTreeSupportSummaryReport:
    return FastTreeSupportSummaryReport(
        tree_path=Path(str(payload["tree_path"])),
        internal_node_count=int(payload["internal_node_count"]),
        annotated_node_count=int(payload["annotated_node_count"]),
        minimum_local_support=(
            None
            if payload["minimum_local_support"] is None
            else float(payload["minimum_local_support"])
        ),
        maximum_local_support=(
            None
            if payload["maximum_local_support"] is None
            else float(payload["maximum_local_support"])
        ),
        median_local_support=(
            None
            if payload["median_local_support"] is None
            else float(payload["median_local_support"])
        ),
        weakly_supported_clade_count=int(payload["weakly_supported_clade_count"]),
        support_histogram={
            str(key): int(value)
            for key, value in dict(payload["support_histogram"]).items()
        },
        approximate_method=bool(payload["approximate_method"]),
        support_label_kind=str(payload["support_label_kind"]),
        support_scale=str(payload["support_scale"]),
        nodes=[
            FastTreeSupportNode(
                node=str(item["node"]),
                descendant_taxa=[str(taxon) for taxon in item["descendant_taxa"]],
                local_support=float(item["local_support"]),
                support_fraction=float(item["support_fraction"]),
                is_backbone=bool(item["is_backbone"]),
            )
            for item in payload["nodes"]
        ],
        warnings=[str(item) for item in payload["warnings"]],
    )


def _restore_large_alignment_report(
    payload: dict[str, object],
) -> LargeAlignmentInferenceWorkflowReport:
    return LargeAlignmentInferenceWorkflowReport(
        workflow=str(payload.get("workflow", "large-alignment-inference")),
        input_path=Path(str(payload["input_path"])),
        out_dir=Path(str(payload["out_dir"])),
        prefix=str(payload["prefix"]),
        sequence_type=str(payload["sequence_type"]),
        engine_name=str(payload["engine_name"]),
        approximate_method=bool(payload["approximate_method"]),
        timeout_seconds=(
            None
            if payload["timeout_seconds"] is None
            else float(payload["timeout_seconds"])
        ),
        resumed=bool(payload.get("resumed", False)),
        started_at_utc=str(payload.get("started_at_utc", "")),
        ended_at_utc=str(payload.get("ended_at_utc", "")),
        runtime_seconds=float(payload.get("runtime_seconds", 0.0)),
        manifest_path=Path(str(payload["manifest_path"])),
        output_paths={
            key: Path(str(value))
            for key, value in dict(payload["output_paths"]).items()
        },
        input_checksums={
            str(key): str(value)
            for key, value in dict(payload["input_checksums"]).items()
        },
        output_checksums={
            str(key): str(value)
            for key, value in dict(payload["output_checksums"]).items()
        },
        config={
            str(key): value for key, value in dict(payload.get("config", {})).items()
        },
        input_summary=LargeAlignmentInputSummary(
            input_path=Path(str(dict(payload["input_summary"])["input_path"])),
            input_bytes=int(dict(payload["input_summary"])["input_bytes"]),
            sequence_count=int(dict(payload["input_summary"])["sequence_count"]),
            alignment_length=int(dict(payload["input_summary"])["alignment_length"]),
            min_sequence_length=int(
                dict(payload["input_summary"])["min_sequence_length"]
            ),
            max_sequence_length=int(
                dict(payload["input_summary"])["max_sequence_length"]
            ),
            total_site_cells=int(dict(payload["input_summary"])["total_site_cells"]),
            inferred_sequence_type=str(
                dict(payload["input_summary"])["inferred_sequence_type"]
            ),
        ),
        resource_rows=[
            LargeAlignmentResourceRow(
                stage=str(item["stage"]),
                elapsed_seconds=float(item["elapsed_seconds"]),
                peak_memory_bytes=(
                    None
                    if item["peak_memory_bytes"] is None
                    else int(item["peak_memory_bytes"])
                ),
                memory_observation_kind=str(item["memory_observation_kind"]),
            )
            for item in payload["resource_rows"]
        ],
        engine_version_text=str(payload["engine_version_text"]),
        command=[str(item) for item in payload["command"]],
        support_summary=_restore_support_summary(dict(payload["support_summary"])),
        warnings=[str(item) for item in payload["warnings"]],
        notes=[str(item) for item in payload["notes"]],
    )


def _resume_existing_large_alignment_report(
    *,
    manifest_path: Path,
    input_path: Path,
    expected_command: list[str],
    expected_timeout_seconds: float | None,
) -> LargeAlignmentInferenceWorkflowReport | None:
    if not manifest_path.exists():
        return None
    payload = load_engine_manifest(manifest_path)
    report = _restore_large_alignment_report(payload)
    if report.command != expected_command:
        return None
    if report.timeout_seconds != expected_timeout_seconds:
        return None
    if report.input_checksums != build_file_checksums([input_path]):
        return None
    required_outputs = [
        report.output_paths["tree"],
        report.output_paths["support_table"],
        report.output_paths["low_support_branches"],
        report.output_paths["support_histogram"],
        report.output_paths["resource_table"],
        report.output_paths["stderr_log"],
        report.output_paths["log"],
    ]
    if any(not path.exists() for path in required_outputs):
        return None
    current_checksums = build_file_checksums(required_outputs)
    if report.output_checksums != current_checksums:
        return None
    clear_incomplete_engine_run(manifest_path)
    report.resumed = True
    return report


def run_large_alignment_inference(
    input_path: Path,
    *,
    out_dir: Path,
    prefix: str = "large-alignment-inference",
    sequence_type: AlignmentAlphabet | None = None,
    executable: str | Path = "FastTree",
    timeout_seconds: float | None = None,
    resume: bool = False,
    incomplete_run_policy: str = "reject",
) -> LargeAlignmentInferenceWorkflowReport:
    """Run governed large-alignment inference with streaming preflight and resource reporting."""
    validate_timeout_seconds(timeout_seconds)
    require_external_engine_surface(
        workflow_id="large-alignment-inference",
        summary="FastTree large-alignment inference workflow.",
        required_engines=("fasttree",),
        executables={"fasttree": executable},
        preserve_missing_error=True,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    root = out_dir / prefix
    tree_path = _tree_path(root)
    support_table_path = root.with_suffix(".support.tsv")
    low_support_branches_path = root.with_suffix(".low-support.tsv")
    support_histogram_path = root.with_suffix(".support-histogram.tsv")
    resource_table_path = _resource_table_path(root)
    stderr_log_path = _stderr_log_path(root)
    plain_log_path = _plain_log_path(root)
    manifest_path = _manifest_path(root)

    input_summary, preflight_resource = _scan_alignment_streaming(input_path)
    effective_sequence_type = (
        input_summary.inferred_sequence_type if sequence_type is None else sequence_type
    )
    version = read_engine_version("FastTree", executable, version_args=("-help",))
    resolved = resolve_engine_executable(executable)
    command = [resolved, *_fasttree_args(input_path.resolve(), effective_sequence_type)]
    if resume:
        resumed = _resume_existing_large_alignment_report(
            manifest_path=manifest_path,
            input_path=input_path,
            expected_command=command,
            expected_timeout_seconds=timeout_seconds,
        )
        if resumed is not None:
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )

    incomplete_record = EngineIncompleteRunRecord(
        engine_name="FastTree",
        workflow="large-alignment-inference",
        executable=resolved,
        working_directory=out_dir,
        manifest_path=manifest_path,
        command=command,
        stdout_path=tree_path,
        stderr_path=stderr_log_path,
        output_paths={
            "tree": tree_path,
            "support_table": support_table_path,
            "low_support_branches": low_support_branches_path,
            "support_histogram": support_histogram_path,
            "resource_table": resource_table_path,
            "log": plain_log_path,
        },
        started_at_utc=utc_now_text(),
        ended_at_utc=None,
        timeout_seconds=timeout_seconds,
        timed_out=False,
        exit_code=None,
        failure_reason=None,
        failure_message=None,
    )
    write_incomplete_engine_run(incomplete_record)

    started_at_utc = incomplete_record.started_at_utc
    started = time.perf_counter()
    try:
        warning_lines, inference_resource = _execute_large_fasttree_command(
            executable=resolved,
            command_args=command[1:],
            out_path=tree_path,
            stderr_path=stderr_log_path,
            timeout_seconds=timeout_seconds,
        )
        _validate_tree_output(
            tree_path,
            engine_name="FastTree",
            workflow="large-alignment-inference",
            output_name="tree",
            artifact_kind="large-alignment-tree",
        )
    except EngineWorkflowError as error:
        incomplete_record.ended_at_utc = utc_now_text()
        incomplete_record.timed_out = "timed out after" in str(error)
        incomplete_record.failure_message = str(error)
        write_incomplete_engine_run(incomplete_record)
        raise
    support_summary = summarize_fasttree_support_distribution(tree_path)
    _validate_support_value_count(
        engine_name="FastTree",
        workflow="large-alignment-inference",
        path=tree_path,
        output_name="tree",
        artifact_kind="large-alignment-tree",
        support_value_count=support_summary.annotated_node_count,
        support_kind="FastTree local support",
    )
    write_fasttree_support_table(
        support_table_path,
        build_fasttree_support_rows(support_summary),
    )
    write_fasttree_support_table(
        low_support_branches_path,
        build_fasttree_low_support_rows(support_summary),
    )
    write_fasttree_support_histogram(
        support_histogram_path,
        build_fasttree_support_histogram_rows(support_summary),
    )
    resource_rows = [preflight_resource, inference_resource]
    write_large_alignment_resource_table(resource_table_path, resource_rows)

    warnings = list(dict.fromkeys(warning_lines + support_summary.warnings))
    notes = [
        "input alignment was scanned linearly before inference so the workflow can validate size and rectangularity without materializing the full matrix in Python memory",
        "large-alignment inference reuses the provided aligned FASTA in place and does not create an intermediate alignment copy before FastTree runs",
        "resource ledger records streamed preflight allocations and sampled FastTree process RSS during inference",
        "FastTree is an approximately maximum-likelihood engine and should be reviewed as exploratory or large-alignment evidence",
        "timeout control applies only to the FastTree inference step and a completed manifest can be resumed when outputs still match the recorded checksums",
        *incomplete_notes,
    ]
    report = LargeAlignmentInferenceWorkflowReport(
        workflow="large-alignment-inference",
        input_path=input_path,
        out_dir=out_dir,
        prefix=prefix,
        sequence_type=effective_sequence_type,
        engine_name="FastTree",
        approximate_method=True,
        timeout_seconds=timeout_seconds,
        resumed=False,
        started_at_utc=started_at_utc,
        ended_at_utc=utc_now_text(),
        runtime_seconds=max(0.0, round(time.perf_counter() - started, 6)),
        manifest_path=manifest_path,
        output_paths={
            "tree": tree_path,
            "support_table": support_table_path,
            "low_support_branches": low_support_branches_path,
            "support_histogram": support_histogram_path,
            "resource_table": resource_table_path,
            "stderr_log": stderr_log_path,
            "log": plain_log_path,
            "manifest": manifest_path,
        },
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        config={
            "sequence_type": effective_sequence_type,
            "timeout_seconds": timeout_seconds,
            "resume": resume,
            "incomplete_run_policy": incomplete_run_policy,
        },
        input_summary=input_summary,
        resource_rows=resource_rows,
        engine_version_text=version.text,
        command=command,
        support_summary=support_summary,
        warnings=warnings,
        notes=notes,
    )
    write_large_alignment_inference_log(plain_log_path, report)
    report.output_checksums = build_file_checksums(
        [
            tree_path,
            support_table_path,
            low_support_branches_path,
            support_histogram_path,
            resource_table_path,
            stderr_log_path,
            plain_log_path,
        ]
    )
    clear_incomplete_engine_run(manifest_path)
    write_engine_manifest(manifest_path, report)
    return report
