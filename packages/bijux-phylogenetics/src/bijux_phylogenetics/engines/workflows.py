from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

from bijux_phylogenetics.compare.reports import ComparisonReportBuildResult, build_tree_comparison_report
from bijux_phylogenetics.core.alignment import AlignmentAlphabet
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.errors import EngineWorkflowError
from bijux_phylogenetics.io.fasta import infer_alignment_alphabet, load_fasta_alignment, summarise_fasta
from bijux_phylogenetics.io.newick import loads_newick

from .common import (
    EngineRunReport,
    execute_engine_command,
    load_unaligned_fasta,
    read_engine_version,
    resolve_engine_executable,
    write_engine_manifest,
)

_BEST_MODEL_PATTERN = re.compile(
    r"(?:best-fit model(?: according to [A-Z0-9]+)?|best model)\s*[:=]\s*(?P<model>[A-Za-z0-9+._-]+)",
    re.IGNORECASE,
)


@dataclass(slots=True)
class EngineWorkflowReport:
    workflow: str
    engine_name: str
    input_paths: list[Path]
    output_paths: dict[str, Path]
    run: EngineRunReport
    manifest_path: Path
    selected_model: str | None = None
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ExternalTreeComparisonReport:
    fast_tree_path: Path
    ml_tree_path: Path
    comparison_report: ComparisonReportBuildResult


def _sidecar(path: Path, label: str) -> Path:
    return path.parent / f"{path.name}.{label}"


def _prefix_path(out_dir: Path, prefix: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / prefix


def _manifest_path_from_output(path: Path) -> Path:
    return _sidecar(path, "manifest.json")


def _validate_alignment_output(path: Path) -> None:
    summarise_fasta(path)


def _validate_tree_output(path: Path) -> None:
    loads_newick(path.read_text(encoding="utf-8"))
    validate_tree_path(path)


def _iqtree_sequence_type_flag(path: Path, sequence_type: AlignmentAlphabet | None) -> list[str]:
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
    for candidate in (prefix_path.with_suffix(".iqtree"), prefix_path.with_suffix(".model")):
        model = _parse_best_model(candidate)
        if model is not None:
            return model
    gz_candidate = prefix_path.with_suffix(".model.gz")
    if gz_candidate.exists():
        text = gzip.decompress(gz_candidate.read_bytes()).decode("utf-8", errors="replace")
        match = _BEST_MODEL_PATTERN.search(text)
        if match is not None:
            return match.group("model")
    return None


def run_multiple_sequence_alignment(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "mafft",
    extra_args: tuple[str, ...] = (),
) -> EngineWorkflowReport:
    """Run a multiple-sequence alignment engine against an unaligned FASTA file."""
    load_unaligned_fasta(input_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    version = read_engine_version("mafft", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    run = execute_engine_command(
        engine_name="mafft",
        workflow="multiple-sequence-alignment",
        executable=resolved,
        version=version,
        command_args=["--auto", *extra_args, str(input_path.resolve())],
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
        notes=["alignment output validated as deterministic equal-length FASTA"],
    )
    write_engine_manifest(manifest_path, report)
    return report


def run_alignment_trimming(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "trimal",
    gap_threshold: float = 0.1,
) -> EngineWorkflowReport:
    """Run an external alignment trimming engine against an aligned FASTA file."""
    load_fasta_alignment(input_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
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
            "-gt",
            f"{gap_threshold:.6f}",
        ],
        work_dir=out_path.parent,
        stdout_path=_sidecar(out_path, "stdout.log"),
        stderr_path=_sidecar(out_path, "stderr.log"),
        output_paths={"trimmed_alignment": out_path},
    )
    _validate_alignment_output(out_path)
    manifest_path = _manifest_path_from_output(out_path)
    report = EngineWorkflowReport(
        workflow="alignment-trimming",
        engine_name="trimal",
        input_paths=[input_path],
        output_paths={"trimmed_alignment": out_path},
        run=run,
        manifest_path=manifest_path,
        notes=["trimmed alignment validated as nonempty equal-length FASTA"],
    )
    write_engine_manifest(manifest_path, report)
    return report


def run_model_selection(
    input_path: Path,
    *,
    out_dir: Path,
    prefix: str = "model-selection",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
) -> EngineWorkflowReport:
    """Run a model-selection workflow on an aligned FASTA file."""
    load_fasta_alignment(input_path)
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
        raise EngineWorkflowError(f"iqtree model-selection did not expose a parsable best-fit model in {iqtree_report_path}")
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
        selected_model=selected_model,
        notes=["best-fit substitution model parsed from engine output"],
    )
    write_engine_manifest(manifest_path, report)
    return report


def run_maximum_likelihood_tree_inference(
    input_path: Path,
    *,
    out_dir: Path,
    model: str,
    prefix: str = "maximum-likelihood",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
) -> EngineWorkflowReport:
    """Run an external maximum-likelihood tree inference workflow."""
    load_fasta_alignment(input_path)
    prefix_path = _prefix_path(out_dir, prefix)
    version = read_engine_version("iqtree", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    tree_path = prefix_path.with_suffix(".treefile")
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="maximum-likelihood-tree",
        executable=resolved,
        version=version,
        command_args=[
            "-s",
            str(input_path.resolve()),
            *(_iqtree_sequence_type_flag(input_path, sequence_type)),
            "-m",
            model,
            "-pre",
            str(prefix_path.resolve()),
        ],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "tree": tree_path,
            "iqtree_report": prefix_path.with_suffix(".iqtree"),
        },
    )
    _validate_tree_output(tree_path)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    report = EngineWorkflowReport(
        workflow="maximum-likelihood-tree",
        engine_name="iqtree",
        input_paths=[input_path],
        output_paths={
            "tree": tree_path,
            "iqtree_report": prefix_path.with_suffix(".iqtree"),
        },
        run=run,
        manifest_path=manifest_path,
        selected_model=model,
        notes=["maximum-likelihood tree validated as parseable Newick output"],
    )
    write_engine_manifest(manifest_path, report)
    return report


def run_bootstrap_support_estimation(
    input_path: Path,
    *,
    out_dir: Path,
    model: str,
    replicates: int = 1000,
    prefix: str = "bootstrap-support",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
) -> EngineWorkflowReport:
    """Run external bootstrap support estimation and retain bootstrap trees."""
    if replicates < 1:
        raise ValueError(f"replicates must be positive, got {replicates}")
    load_fasta_alignment(input_path)
    prefix_path = _prefix_path(out_dir, prefix)
    version = read_engine_version("iqtree", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    support_tree_path = prefix_path.with_suffix(".treefile")
    bootstrap_tree_path = prefix_path.with_suffix(".ufboot")
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="bootstrap-support",
        executable=resolved,
        version=version,
        command_args=[
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
        ],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "support_tree": support_tree_path,
            "bootstrap_trees": bootstrap_tree_path,
            "iqtree_report": prefix_path.with_suffix(".iqtree"),
        },
    )
    _validate_tree_output(support_tree_path)
    if not bootstrap_tree_path.read_text(encoding="utf-8").strip():
        raise EngineWorkflowError(f"bootstrap tree set is empty: {bootstrap_tree_path}")
    manifest_path = prefix_path.with_suffix(".manifest.json")
    report = EngineWorkflowReport(
        workflow="bootstrap-support",
        engine_name="iqtree",
        input_paths=[input_path],
        output_paths={
            "support_tree": support_tree_path,
            "bootstrap_trees": bootstrap_tree_path,
            "iqtree_report": prefix_path.with_suffix(".iqtree"),
        },
        run=run,
        manifest_path=manifest_path,
        selected_model=model,
        notes=["bootstrap tree set retained for downstream consensus construction"],
    )
    write_engine_manifest(manifest_path, report)
    return report


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
        raise ValueError(f"minimum_support must be between 0 and 1 inclusive, got {minimum_support}")
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
        notes=["consensus tree validated as parseable Newick output"],
    )
    write_engine_manifest(manifest_path, report)
    return report


def run_fast_tree_inference(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "FastTree",
    sequence_type: AlignmentAlphabet | None = None,
) -> EngineWorkflowReport:
    """Run a fast approximate tree inference engine against an aligned FASTA file."""
    load_fasta_alignment(input_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    version = read_engine_version("FastTree", executable, version_args=("-help",))
    resolved = resolve_engine_executable(executable)
    run = execute_engine_command(
        engine_name="FastTree",
        workflow="fast-approximate-tree",
        executable=resolved,
        version=version,
        command_args=_fasttree_args(input_path.resolve(), sequence_type),
        work_dir=out_path.parent,
        stdout_path=out_path,
        stderr_path=_sidecar(out_path, "stderr.log"),
        output_paths={"tree": out_path},
    )
    _validate_tree_output(out_path)
    manifest_path = _manifest_path_from_output(out_path)
    report = EngineWorkflowReport(
        workflow="fast-approximate-tree",
        engine_name="FastTree",
        input_paths=[input_path],
        output_paths={"tree": out_path},
        run=run,
        manifest_path=manifest_path,
        notes=["fast approximate tree validated as parseable Newick output"],
    )
    write_engine_manifest(manifest_path, report)
    return report


def compare_fast_and_ml_trees(
    fast_tree_path: Path,
    ml_tree_path: Path,
    *,
    out_path: Path,
) -> ExternalTreeComparisonReport:
    """Compare a fast approximate tree against a maximum-likelihood tree."""
    comparison_report = build_tree_comparison_report(fast_tree_path, ml_tree_path, out_path=out_path)
    return ExternalTreeComparisonReport(
        fast_tree_path=fast_tree_path,
        ml_tree_path=ml_tree_path,
        comparison_report=comparison_report,
    )
