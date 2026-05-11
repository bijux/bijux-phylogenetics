from __future__ import annotations

from dataclasses import dataclass, field
import gzip
from pathlib import Path
import re

from bijux_phylogenetics.compare.reports import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.core.alignment import AlignmentAlphabet, AlignmentSummary
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.errors import EngineWorkflowError
from bijux_phylogenetics.io.fasta import (
    infer_alignment_alphabet,
    load_fasta_alignment,
    summarise_fasta,
)
from bijux_phylogenetics.io.newick import loads_newick

from .common import (
    EngineRunReport,
    EngineVersionInfo,
    build_file_checksums,
    execute_engine_command,
    load_engine_manifest,
    load_unaligned_fasta,
    read_engine_version,
    resolve_engine_executable,
    write_engine_manifest,
)

_BEST_MODEL_PATTERN = re.compile(
    r"(?:best-fit model(?: according to [A-Z0-9]+)?|best model)\s*[:=]\s*(?P<model>[A-Za-z0-9+._-]+)",
    re.IGNORECASE,
)
_MAFFT_ALIGNMENT_MODE_ARGUMENTS: dict[str, tuple[str, ...]] = {
    "auto": ("--auto",),
    "linsi": ("--localpair", "--maxiterate", "1000"),
    "ginsi": ("--globalpair", "--maxiterate", "1000"),
    "einsi": ("--ep", "0", "--genafpair", "--maxiterate", "1000"),
    "fast": ("--retree", "2", "--maxiterate", "0"),
}
_TRIMAL_TRIMMING_MODES: tuple[str, ...] = (
    "gap-threshold",
    "gappyout",
    "strict",
    "strictplus",
    "automated1",
)


@dataclass(slots=True)
class EngineWorkflowReport:
    workflow: str
    engine_name: str
    input_paths: list[Path]
    output_paths: dict[str, Path]
    run: EngineRunReport
    manifest_path: Path
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    selected_model: str | None = None
    trimming_summary: AlignmentTrimmingSummary | None = None
    resumed: bool = False
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class AlignmentTrimmingSummary:
    mode: str
    gap_threshold: float | None
    input_alignment_length: int
    trimmed_alignment_length: int
    retained_site_count: int
    removed_site_count: int
    retained_site_fraction: float
    removed_site_fraction: float
    input_gap_fraction: float
    trimmed_gap_fraction: float
    input_gap_percentage: float
    trimmed_gap_percentage: float


@dataclass(slots=True)
class ExternalTreeComparisonReport:
    fast_tree_path: Path
    ml_tree_path: Path
    comparison_report: ComparisonReportBuildResult


def list_mafft_alignment_modes() -> tuple[str, ...]:
    """Return the supported named MAFFT alignment strategies."""
    return tuple(_MAFFT_ALIGNMENT_MODE_ARGUMENTS)


def resolve_mafft_alignment_mode(mode: str) -> tuple[str, ...]:
    """Resolve one named MAFFT alignment strategy into explicit engine arguments."""
    try:
        return _MAFFT_ALIGNMENT_MODE_ARGUMENTS[mode]
    except KeyError as error:
        available = ", ".join(sorted(_MAFFT_ALIGNMENT_MODE_ARGUMENTS))
        raise ValueError(
            f"unsupported mafft alignment mode '{mode}', expected one of: {available}"
        ) from error


def list_trimal_trimming_modes() -> tuple[str, ...]:
    """Return the supported named trimAl trimming strategies."""
    return _TRIMAL_TRIMMING_MODES


def resolve_trimal_trimming_mode(
    mode: str,
    *,
    gap_threshold: float,
) -> tuple[str, ...]:
    """Resolve one named trimAl strategy into explicit engine arguments."""
    if mode == "gap-threshold":
        return ("-gt", f"{gap_threshold:.6f}")
    if mode in {"gappyout", "strict", "strictplus", "automated1"}:
        return (f"-{mode}",)
    available = ", ".join(sorted(_TRIMAL_TRIMMING_MODES))
    raise ValueError(
        f"unsupported trimAl trimming mode '{mode}', expected one of: {available}"
    )


def _sidecar(path: Path, label: str) -> Path:
    return path.parent / f"{path.name}.{label}"


def _prefix_path(out_dir: Path, prefix: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / prefix


def _manifest_path_from_output(path: Path) -> Path:
    return _sidecar(path, "manifest.json")


def _validate_alignment_output(path: Path) -> None:
    records = load_fasta_alignment(path)
    if not records or len(records[0].sequence) < 1:
        raise EngineWorkflowError(
            f"inference alignment is empty after filtering: {path}"
        )
    summarise_fasta(path)


def _build_alignment_trimming_summary(
    *,
    mode: str,
    gap_threshold: float,
    input_summary: AlignmentSummary,
    trimmed_summary: AlignmentSummary,
) -> AlignmentTrimmingSummary:
    if trimmed_summary.alignment_length > input_summary.alignment_length:
        raise EngineWorkflowError(
            "trimmed alignment is longer than the input alignment, which is not a valid trimAl result"
        )
    retained_site_count = trimmed_summary.alignment_length
    removed_site_count = input_summary.alignment_length - trimmed_summary.alignment_length
    retained_site_fraction = retained_site_count / input_summary.alignment_length
    removed_site_fraction = removed_site_count / input_summary.alignment_length
    return AlignmentTrimmingSummary(
        mode=mode,
        gap_threshold=gap_threshold if mode == "gap-threshold" else None,
        input_alignment_length=input_summary.alignment_length,
        trimmed_alignment_length=trimmed_summary.alignment_length,
        retained_site_count=retained_site_count,
        removed_site_count=removed_site_count,
        retained_site_fraction=retained_site_fraction,
        removed_site_fraction=removed_site_fraction,
        input_gap_fraction=input_summary.gap_fraction,
        trimmed_gap_fraction=trimmed_summary.gap_fraction,
        input_gap_percentage=input_summary.gap_fraction * 100.0,
        trimmed_gap_percentage=trimmed_summary.gap_fraction * 100.0,
    )


def _validate_tree_output(path: Path) -> None:
    loads_newick(path.read_text(encoding="utf-8"))
    validate_tree_path(path)


def _ensure_inference_ready_alignment(path: Path) -> None:
    records = load_fasta_alignment(path)
    if not records or len(records[0].sequence) < 1:
        raise EngineWorkflowError(
            f"inference alignment is empty after filtering: {path}"
        )
    summarise_fasta(path)


def _restore_workflow_report(payload: dict[str, object]) -> EngineWorkflowReport:
    run_payload = payload["run"]
    if not isinstance(run_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid run payload")
    version_payload = run_payload["version"]
    if not isinstance(version_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid version payload")
    version = EngineVersionInfo(
        engine_name=str(version_payload["engine_name"]),
        executable=str(version_payload["executable"]),
        command=[str(item) for item in version_payload["command"]],
        text=str(version_payload["text"]),
    )
    run = EngineRunReport(
        engine_name=str(run_payload["engine_name"]),
        workflow=str(run_payload["workflow"]),
        executable=str(run_payload["executable"]),
        working_directory=Path(run_payload["working_directory"]),
        version=version,
        command=[str(item) for item in run_payload["command"]],
        exit_code=int(run_payload["exit_code"]),
        stdout_path=Path(run_payload["stdout_path"]),
        stderr_path=Path(run_payload["stderr_path"]),
        output_paths={
            str(key): Path(value)
            for key, value in dict(run_payload["output_paths"]).items()
        },
        warning_lines=[str(item) for item in run_payload["warning_lines"]],
    )
    return EngineWorkflowReport(
        workflow=str(payload["workflow"]),
        engine_name=str(payload["engine_name"]),
        input_paths=[Path(item) for item in payload["input_paths"]],
        output_paths={
            str(key): Path(value)
            for key, value in dict(payload["output_paths"]).items()
        },
        run=run,
        manifest_path=Path(payload["manifest_path"]),
        input_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("input_checksums", {})).items()
        },
        output_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("output_checksums", {})).items()
        },
        selected_model=None
        if payload.get("selected_model") is None
        else str(payload["selected_model"]),
        trimming_summary=None
        if payload.get("trimming_summary") is None
        else AlignmentTrimmingSummary(
            mode=str(dict(payload["trimming_summary"])["mode"]),
            gap_threshold=(
                None
                if dict(payload["trimming_summary"]).get("gap_threshold") is None
                else float(dict(payload["trimming_summary"])["gap_threshold"])
            ),
            input_alignment_length=int(
                dict(payload["trimming_summary"])["input_alignment_length"]
            ),
            trimmed_alignment_length=int(
                dict(payload["trimming_summary"])["trimmed_alignment_length"]
            ),
            retained_site_count=int(
                dict(payload["trimming_summary"])["retained_site_count"]
            ),
            removed_site_count=int(
                dict(payload["trimming_summary"])["removed_site_count"]
            ),
            retained_site_fraction=float(
                dict(payload["trimming_summary"])["retained_site_fraction"]
            ),
            removed_site_fraction=float(
                dict(payload["trimming_summary"])["removed_site_fraction"]
            ),
            input_gap_fraction=float(
                dict(payload["trimming_summary"])["input_gap_fraction"]
            ),
            trimmed_gap_fraction=float(
                dict(payload["trimming_summary"])["trimmed_gap_fraction"]
            ),
            input_gap_percentage=float(
                dict(payload["trimming_summary"])["input_gap_percentage"]
            ),
            trimmed_gap_percentage=float(
                dict(payload["trimming_summary"])["trimmed_gap_percentage"]
            ),
        ),
        resumed=bool(payload.get("resumed", False)),
        notes=[str(item) for item in payload.get("notes", [])],
    )


def _resume_existing_workflow(
    *,
    manifest_path: Path,
    input_paths: list[Path],
    expected_command: list[str],
) -> EngineWorkflowReport | None:
    if not manifest_path.exists():
        return None
    payload = load_engine_manifest(manifest_path)
    report = _restore_workflow_report(payload)
    if report.run.command != expected_command:
        return None
    current_input_checksums = build_file_checksums(input_paths)
    if report.input_checksums != current_input_checksums:
        return None
    if any(not path.exists() for path in report.output_paths.values()):
        return None
    current_output_checksums = build_file_checksums(list(report.output_paths.values()))
    if report.output_checksums != current_output_checksums:
        return None
    report.resumed = True
    return report


def _persist_workflow_report(report: EngineWorkflowReport) -> EngineWorkflowReport:
    report.output_checksums = build_file_checksums(list(report.output_paths.values()))
    write_engine_manifest(report.manifest_path, report)
    return report


def _iqtree_sequence_type_flag(
    path: Path, sequence_type: AlignmentAlphabet | None
) -> list[str]:
    detected = sequence_type
    if detected is None:
        detected = infer_alignment_alphabet(load_fasta_alignment(path))
    if detected in {"dna", "rna"}:
        return ["-st", "DNA"]
    if detected == "protein":
        return ["-st", "AA"]
    return []


def _fasttree_args(path: Path, sequence_type: AlignmentAlphabet | None) -> list[str]:
    detected = sequence_type
    if detected is None:
        detected = infer_alignment_alphabet(load_fasta_alignment(path))
    if detected in {"dna", "rna"}:
        return ["-gtr", "-nt", str(path)]
    if detected == "protein":
        return ["-lg", str(path)]
    return [str(path)]


def _parse_best_model(iqtree_report_path: Path) -> str | None:
    if not iqtree_report_path.exists():
        return None
    match = _BEST_MODEL_PATTERN.search(iqtree_report_path.read_text(encoding="utf-8"))
    if match is not None:
        return match.group("model")
    return None


def _parse_best_model_artifact(prefix_path: Path) -> str | None:
    for candidate in (
        prefix_path.with_suffix(".iqtree"),
        prefix_path.with_suffix(".model"),
    ):
        model = _parse_best_model(candidate)
        if model is not None:
            return model
    gz_candidate = prefix_path.with_suffix(".model.gz")
    if gz_candidate.exists():
        text = gzip.decompress(gz_candidate.read_bytes()).decode(
            "utf-8", errors="replace"
        )
        match = _BEST_MODEL_PATTERN.search(text)
        if match is not None:
            return match.group("model")
    return None


def run_multiple_sequence_alignment(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "mafft",
    mode: str = "auto",
    extra_args: tuple[str, ...] = (),
) -> EngineWorkflowReport:
    """Run a multiple-sequence alignment engine against an unaligned FASTA file."""
    load_unaligned_fasta(input_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mode_args = resolve_mafft_alignment_mode(mode)
    version = read_engine_version("mafft", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    run = execute_engine_command(
        engine_name="mafft",
        workflow="multiple-sequence-alignment",
        executable=resolved,
        version=version,
        command_args=[*mode_args, *extra_args, str(input_path.resolve())],
        work_dir=out_path.parent,
        stdout_path=out_path,
        stderr_path=_sidecar(out_path, "stderr.log"),
        output_paths={"alignment": out_path},
    )
    _validate_alignment_output(out_path)
    manifest_path = _manifest_path_from_output(out_path)
    report = EngineWorkflowReport(
        workflow="multiple-sequence-alignment",
        engine_name="mafft",
        input_paths=[input_path],
        output_paths={"alignment": out_path},
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        notes=[
            f"mafft alignment mode: {mode}",
            "alignment output validated as deterministic equal-length FASTA",
        ],
    )
    return _persist_workflow_report(report)


def run_alignment_trimming(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "trimal",
    mode: str = "gap-threshold",
    gap_threshold: float = 0.1,
) -> EngineWorkflowReport:
    """Run an external alignment trimming engine against an aligned FASTA file."""
    load_fasta_alignment(input_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    input_summary = summarise_fasta(input_path)
    mode_args = resolve_trimal_trimming_mode(mode, gap_threshold=gap_threshold)
    version = read_engine_version("trimal", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    run = execute_engine_command(
        engine_name="trimal",
        workflow="alignment-trimming",
        executable=resolved,
        version=version,
        command_args=[
            "-in",
            str(input_path.resolve()),
            "-out",
            str(out_path.resolve()),
            *mode_args,
        ],
        work_dir=out_path.parent,
        stdout_path=_sidecar(out_path, "stdout.log"),
        stderr_path=_sidecar(out_path, "stderr.log"),
        output_paths={"trimmed_alignment": out_path},
    )
    _validate_alignment_output(out_path)
    trimmed_summary = summarise_fasta(out_path)
    trimming_summary = _build_alignment_trimming_summary(
        mode=mode,
        gap_threshold=gap_threshold,
        input_summary=input_summary,
        trimmed_summary=trimmed_summary,
    )
    manifest_path = _manifest_path_from_output(out_path)
    report = EngineWorkflowReport(
        workflow="alignment-trimming",
        engine_name="trimal",
        input_paths=[input_path],
        output_paths={"trimmed_alignment": out_path},
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        trimming_summary=trimming_summary,
        notes=[
            f"trimal trimming mode: {mode}",
            f"retained sites: {trimming_summary.retained_site_count} of {trimming_summary.input_alignment_length}",
            f"gap percentage: {trimming_summary.input_gap_percentage:.3f} -> {trimming_summary.trimmed_gap_percentage:.3f}",
            "trimmed alignment validated as nonempty equal-length FASTA",
        ],
    )
    return _persist_workflow_report(report)


def run_model_selection(
    input_path: Path,
    *,
    out_dir: Path,
    prefix: str = "model-selection",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
) -> EngineWorkflowReport:
    """Run a model-selection workflow on an aligned FASTA file."""
    _ensure_inference_ready_alignment(input_path)
    prefix_path = _prefix_path(out_dir, prefix)
    version = read_engine_version("iqtree", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    iqtree_report_path = prefix_path.with_suffix(".iqtree")
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="model-selection",
        executable=resolved,
        version=version,
        command_args=[
            "-s",
            str(input_path.resolve()),
            *(_iqtree_sequence_type_flag(input_path, sequence_type)),
            "-m",
            "MF",
            "-pre",
            str(prefix_path.resolve()),
        ],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={"iqtree_report": iqtree_report_path},
    )
    selected_model = _parse_best_model_artifact(prefix_path)
    if selected_model is None:
        raise EngineWorkflowError(
            f"iqtree model-selection did not expose a parsable best-fit model in {iqtree_report_path}"
        )
    selected_model_path = prefix_path.with_suffix(".selected-model.txt")
    selected_model_path.write_text(selected_model + "\n", encoding="utf-8")
    manifest_path = prefix_path.with_suffix(".manifest.json")
    report = EngineWorkflowReport(
        workflow="model-selection",
        engine_name="iqtree",
        input_paths=[input_path],
        output_paths={
            "iqtree_report": iqtree_report_path,
            "selected_model": selected_model_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        selected_model=selected_model,
        notes=["best-fit substitution model parsed from engine output"],
    )
    return _persist_workflow_report(report)


def run_maximum_likelihood_tree_inference(
    input_path: Path,
    *,
    out_dir: Path,
    model: str,
    prefix: str = "maximum-likelihood",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    resume: bool = False,
) -> EngineWorkflowReport:
    """Run an external maximum-likelihood tree inference workflow."""
    _ensure_inference_ready_alignment(input_path)
    prefix_path = _prefix_path(out_dir, prefix)
    version = read_engine_version("iqtree", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    tree_path = prefix_path.with_suffix(".treefile")
    report_path = prefix_path.with_suffix(".iqtree")
    manifest_path = prefix_path.with_suffix(".manifest.json")
    command = [
        resolved,
        "-s",
        str(input_path.resolve()),
        *(_iqtree_sequence_type_flag(input_path, sequence_type)),
        "-m",
        model,
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[input_path],
            expected_command=command,
        )
        if resumed is not None:
            return resumed
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="maximum-likelihood-tree",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "tree": tree_path,
            "iqtree_report": report_path,
        },
    )
    _validate_tree_output(tree_path)
    report = EngineWorkflowReport(
        workflow="maximum-likelihood-tree",
        engine_name="iqtree",
        input_paths=[input_path],
        output_paths={
            "tree": tree_path,
            "iqtree_report": report_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        selected_model=model,
        notes=["maximum-likelihood tree validated as parseable Newick output"],
    )
    return _persist_workflow_report(report)


def run_bootstrap_support_estimation(
    input_path: Path,
    *,
    out_dir: Path,
    model: str,
    replicates: int = 1000,
    prefix: str = "bootstrap-support",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    resume: bool = False,
) -> EngineWorkflowReport:
    """Run external bootstrap support estimation and retain bootstrap trees."""
    if replicates < 1:
        raise ValueError(f"replicates must be positive, got {replicates}")
    _ensure_inference_ready_alignment(input_path)
    prefix_path = _prefix_path(out_dir, prefix)
    version = read_engine_version("iqtree", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    support_tree_path = prefix_path.with_suffix(".treefile")
    bootstrap_tree_path = prefix_path.with_suffix(".ufboot")
    report_path = prefix_path.with_suffix(".iqtree")
    manifest_path = prefix_path.with_suffix(".manifest.json")
    command = [
        resolved,
        "-s",
        str(input_path.resolve()),
        *(_iqtree_sequence_type_flag(input_path, sequence_type)),
        "-m",
        model,
        "-bb",
        str(replicates),
        "-wbt",
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[input_path],
            expected_command=command,
        )
        if resumed is not None:
            return resumed
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="bootstrap-support",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "support_tree": support_tree_path,
            "bootstrap_trees": bootstrap_tree_path,
            "iqtree_report": report_path,
        },
    )
    _validate_tree_output(support_tree_path)
    if not bootstrap_tree_path.read_text(encoding="utf-8").strip():
        raise EngineWorkflowError(f"bootstrap tree set is empty: {bootstrap_tree_path}")
    report = EngineWorkflowReport(
        workflow="bootstrap-support",
        engine_name="iqtree",
        input_paths=[input_path],
        output_paths={
            "support_tree": support_tree_path,
            "bootstrap_trees": bootstrap_tree_path,
            "iqtree_report": report_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        selected_model=model,
        notes=["bootstrap tree set retained for downstream consensus construction"],
    )
    return _persist_workflow_report(report)


def run_bootstrap_consensus_tree(
    bootstrap_trees_path: Path,
    *,
    out_dir: Path,
    prefix: str = "bootstrap-consensus",
    executable: str | Path = "iqtree2",
    minimum_support: float = 0.5,
) -> EngineWorkflowReport:
    """Construct a consensus tree from bootstrap trees."""
    if not 0.0 <= minimum_support <= 1.0:
        raise ValueError(
            f"minimum_support must be between 0 and 1 inclusive, got {minimum_support}"
        )
    if not bootstrap_trees_path.exists():
        raise FileNotFoundError(bootstrap_trees_path)
    prefix_path = _prefix_path(out_dir, prefix)
    version = read_engine_version("iqtree", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    consensus_tree_path = prefix_path.with_suffix(".contree")
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="bootstrap-consensus",
        executable=resolved,
        version=version,
        command_args=[
            "-t",
            str(bootstrap_trees_path.resolve()),
            "-con",
            "-minsup",
            str(minimum_support),
            "-pre",
            str(prefix_path.resolve()),
        ],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={"consensus_tree": consensus_tree_path},
    )
    _validate_tree_output(consensus_tree_path)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    report = EngineWorkflowReport(
        workflow="bootstrap-consensus",
        engine_name="iqtree",
        input_paths=[bootstrap_trees_path],
        output_paths={"consensus_tree": consensus_tree_path},
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([bootstrap_trees_path]),
        output_checksums={},
        notes=["consensus tree validated as parseable Newick output"],
    )
    return _persist_workflow_report(report)


def run_fast_tree_inference(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "FastTree",
    sequence_type: AlignmentAlphabet | None = None,
    resume: bool = False,
) -> EngineWorkflowReport:
    """Run a fast approximate tree inference engine against an aligned FASTA file."""
    _ensure_inference_ready_alignment(input_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    version = read_engine_version("FastTree", executable, version_args=("-help",))
    resolved = resolve_engine_executable(executable)
    manifest_path = _manifest_path_from_output(out_path)
    command = [resolved, *_fasttree_args(input_path.resolve(), sequence_type)]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[input_path],
            expected_command=command,
        )
        if resumed is not None:
            return resumed
    run = execute_engine_command(
        engine_name="FastTree",
        workflow="fast-approximate-tree",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_path.parent,
        stdout_path=out_path,
        stderr_path=_sidecar(out_path, "stderr.log"),
        output_paths={"tree": out_path},
    )
    _validate_tree_output(out_path)
    report = EngineWorkflowReport(
        workflow="fast-approximate-tree",
        engine_name="FastTree",
        input_paths=[input_path],
        output_paths={"tree": out_path},
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        notes=["fast approximate tree validated as parseable Newick output"],
    )
    return _persist_workflow_report(report)


def compare_fast_and_ml_trees(
    fast_tree_path: Path,
    ml_tree_path: Path,
    *,
    out_path: Path,
) -> ExternalTreeComparisonReport:
    """Compare a fast approximate tree against a maximum-likelihood tree."""
    comparison_report = build_tree_comparison_report(
        fast_tree_path, ml_tree_path, out_path=out_path
    )
    return ExternalTreeComparisonReport(
        fast_tree_path=fast_tree_path,
        ml_tree_path=ml_tree_path,
        comparison_report=comparison_report,
    )
