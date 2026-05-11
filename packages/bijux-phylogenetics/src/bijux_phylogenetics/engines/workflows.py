from __future__ import annotations

from dataclasses import dataclass, field
import gzip
import hashlib
from pathlib import Path
import re

from bijux_phylogenetics.compare.reports import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.core.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    AlignmentSummary,
    CodingSequenceExclusion,
)
from bijux_phylogenetics.core.partitions import (
    LocusPartition,
    PartitionSummaryReport,
    build_partition_summary_report,
    normalize_partition_data_type,
    parse_locus_partitions,
    slice_partition_sequence,
    write_locus_partitions,
    write_partition_summary_table,
)
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.errors import EngineWorkflowError
from bijux_phylogenetics.io.fasta import (
    back_translate_aligned_coding_sequences,
    infer_alignment_alphabet,
    load_fasta_alignment,
    prepare_coding_sequences_for_alignment,
    summarise_fasta,
    translate_prepared_coding_sequences,
    write_fasta_alignment,
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
from .validation import summarize_bootstrap_support_distribution

_BEST_MODEL_PATTERN = re.compile(
    r"(?:best-fit model(?: according to [A-Z0-9]+)?|best model)\s*[:=]\s*(?P<model>[A-Za-z0-9+._-]+)",
    re.IGNORECASE,
)
_LOG_LIKELIHOOD_PATTERNS = (
    re.compile(
        r"(?:log-likelihood(?: of the tree)?|log likelihood)\s*[:=]\s*(?P<value>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)",
        re.IGNORECASE,
    ),
    re.compile(
        r"best score found\s*[:=]\s*(?P<value>-?[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)",
        re.IGNORECASE,
    ),
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
_MINIMUM_UFBOOT_REPLICATES = 1000


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
    log_likelihood: float | None = None
    iqtree_summary: IqtreeWorkflowSummary | None = None
    trimming_summary: AlignmentTrimmingSummary | None = None
    resumed: bool = False
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class IqtreeSupportValue:
    node: str
    descendant_taxa: list[str]
    support: float
    support_fraction: float
    is_backbone: bool


@dataclass(slots=True)
class IqtreeWorkflowSummary:
    selected_model: str | None
    log_likelihood: float | None
    support_value_count: int
    minimum_support: float | None
    maximum_support: float | None
    support_values: list[IqtreeSupportValue] = field(default_factory=list)


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


@dataclass(slots=True)
class CodonAwareAlignmentWorkflowReport:
    workflow: str
    engine_name: str
    input_path: Path
    output_paths: dict[str, Path]
    run: EngineRunReport
    manifest_path: Path
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    sequence_type: AlignmentAlphabet
    accepted_sequence_count: int
    excluded_sequences: list[CodingSequenceExclusion]
    terminal_stop_sequence_count: int
    notes: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


@dataclass(frozen=True, slots=True)
class _PreparedIqtreePartitions:
    command_args: list[str]
    summary: PartitionSummaryReport
    output_paths: dict[str, Path]
    notes: list[str]
    mixed_data_types: bool


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


def _partition_support_path(prefix_path: Path, suffix: str) -> Path:
    return prefix_path.parent / f"{prefix_path.name}.{suffix}"


def _validate_alignment_output(path: Path) -> None:
    records = load_fasta_alignment(path)
    if not records or len(records[0].sequence) < 1:
        raise EngineWorkflowError(
            f"inference alignment is empty after filtering: {path}"
        )
    summarise_fasta(path)


def _write_coding_exclusion_table(
    path: Path, exclusions: list[CodingSequenceExclusion]
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "identifier",
                "comparable_length",
                "reason",
                "premature_stop_count",
                "terminal_stop_count",
                "trailing_bases",
                "note",
            ]
        )
    ]
    lines.extend(
        "\t".join(
            [
                row.identifier,
                str(row.comparable_length),
                row.reason,
                str(row.premature_stop_count),
                str(row.terminal_stop_count),
                str(row.trailing_bases),
                row.note,
            ]
        )
        for row in exclusions
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


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


def _partition_alignment_file_name(partition: LocusPartition) -> str:
    normalized = re.sub(r"[^A-Za-z0-9._-]+", "-", partition.name.strip().lower())
    normalized = normalized.strip("-._") or "partition"
    digest = hashlib.sha1(partition.name.encode("utf-8")).hexdigest()[:8]  # nosec B324
    return f"{normalized}-{digest}.fasta"


def _iqtree_partition_supports_fixed_model(
    *,
    model: str,
    mixed_data_types: bool,
) -> bool:
    if not mixed_data_types:
        return True
    normalized = model.strip().upper()
    return normalized in {
        "TEST",
        "TESTONLY",
        "TESTNEWONLY",
        "MF",
        "MFP",
        "TESTMERGE",
        "TESTMERGEONLY",
        "MF+MERGE",
        "MFP+MERGE",
    }


def _prepare_iqtree_partitions(
    input_path: Path,
    partition_path: Path,
    *,
    prefix_path: Path,
) -> _PreparedIqtreePartitions:
    records = load_fasta_alignment(input_path)
    alignment_summary = summarise_fasta(input_path)
    partitions = parse_locus_partitions(partition_path)
    summary = build_partition_summary_report(
        partitions,
        alignment_length=alignment_summary.alignment_length,
    )
    summary_path = _partition_support_path(prefix_path, "partition-summary.tsv")
    write_partition_summary_table(summary_path, summary)
    notes = [
        f"validated {summary.partition_count} partitions across {summary.assigned_site_count} assigned sites",
    ]
    output_paths: dict[str, Path] = {
        "partition_summary": summary_path,
    }

    declared_types = {
        normalize_partition_data_type(partition.data_type)
        for partition in partitions
        if partition.data_type is not None
    }
    if len(declared_types) <= 1:
        normalized_partition_path = _partition_support_path(
            prefix_path, "partition-scheme.partitions"
        )
        write_locus_partitions(normalized_partition_path, partitions)
        output_paths["partition_scheme"] = normalized_partition_path
        notes.append("prepared a normalized partition scheme for single-alignment IQ-TREE analysis")
        return _PreparedIqtreePartitions(
            command_args=[
                "-s",
                str(input_path.resolve()),
                "-p",
                str(normalized_partition_path.resolve()),
            ],
            summary=summary,
            output_paths=output_paths,
            notes=notes,
            mixed_data_types=False,
        )

    if any(partition.data_type is None for partition in partitions):
        raise EngineWorkflowError(
            "mixed partition analyses require every partition to declare a data_type"
        )
    unsupported_types = sorted(
        {
            data_type
            for data_type in declared_types
            if data_type not in {"DNA", "RNA", "PROTEIN"}
        }
    )
    if unsupported_types:
        raise EngineWorkflowError(
            "mixed partition analyses currently support only DNA, RNA, and PROTEIN datatypes; "
            f"got: {', '.join(unsupported_types)}"
        )

    partition_alignment_dir = _partition_support_path(prefix_path, "partition-alignments")
    partition_alignment_dir.mkdir(parents=True, exist_ok=True)
    lines = ["#nexus", "begin sets;"]
    for partition in partitions:
        partition_alignment_path = (
            partition_alignment_dir / _partition_alignment_file_name(partition)
        )
        write_fasta_alignment(
            partition_alignment_path,
            [
                AlignmentRecord(
                    identifier=record.identifier,
                    sequence=slice_partition_sequence(record.sequence, partition),
                )
                for record in records
            ],
        )
        output_paths[f"partition_alignment_{partition.name}"] = partition_alignment_path
        lines.append(
            f"    charset {partition.name} = {partition_alignment_path.name}: *;"
        )
    lines.append("end;")
    mixed_partition_path = _partition_support_path(prefix_path, "partition-scheme.nex")
    mixed_partition_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    output_paths["partition_scheme"] = mixed_partition_path
    notes.append(
        "prepared a mixed-datatype NEXUS partition scheme with one extracted alignment per partition"
    )
    return _PreparedIqtreePartitions(
        command_args=["-p", str(mixed_partition_path.resolve())],
        summary=summary,
        output_paths=output_paths,
        notes=notes,
        mixed_data_types=True,
    )


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
        log_likelihood=(
            None
            if payload.get("log_likelihood") is None
            else float(payload["log_likelihood"])
        ),
        iqtree_summary=(
            None
            if payload.get("iqtree_summary") is None
            else IqtreeWorkflowSummary(
                selected_model=(
                    None
                    if dict(payload["iqtree_summary"]).get("selected_model") is None
                    else str(dict(payload["iqtree_summary"])["selected_model"])
                ),
                log_likelihood=(
                    None
                    if dict(payload["iqtree_summary"]).get("log_likelihood") is None
                    else float(dict(payload["iqtree_summary"])["log_likelihood"])
                ),
                support_value_count=int(
                    dict(payload["iqtree_summary"])["support_value_count"]
                ),
                minimum_support=(
                    None
                    if dict(payload["iqtree_summary"]).get("minimum_support") is None
                    else float(dict(payload["iqtree_summary"])["minimum_support"])
                ),
                maximum_support=(
                    None
                    if dict(payload["iqtree_summary"]).get("maximum_support") is None
                    else float(dict(payload["iqtree_summary"])["maximum_support"])
                ),
                support_values=[
                    IqtreeSupportValue(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon)
                            for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        support=float(dict(item)["support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(dict(payload["iqtree_summary"])["support_values"])
                ],
            )
        ),
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


def _iqtree_execution_controls(*, seed: int, threads: int) -> list[str]:
    if seed < 1:
        raise ValueError(f"iqtree seed must be positive, got {seed}")
    if threads < 1:
        raise ValueError(f"iqtree threads must be positive, got {threads}")
    return ["-seed", str(seed), "-nt", str(threads)]


def _validate_ufboot_replicates(replicates: int) -> None:
    if replicates < _MINIMUM_UFBOOT_REPLICATES:
        raise EngineWorkflowError(
            "iqtree ultrafast bootstrap requires at least "
            f"{_MINIMUM_UFBOOT_REPLICATES} replicates, got {replicates}"
        )


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


def _parse_log_likelihood_text(text: str) -> float | None:
    for pattern in _LOG_LIKELIHOOD_PATTERNS:
        match = pattern.search(text)
        if match is None:
            continue
        try:
            return float(match.group("value"))
        except ValueError:
            continue
    return None


def _parse_log_likelihood_artifact(prefix_path: Path) -> float | None:
    for candidate in (
        prefix_path.with_suffix(".iqtree"),
        prefix_path.with_suffix(".log"),
    ):
        if not candidate.exists():
            continue
        log_likelihood = _parse_log_likelihood_text(
            candidate.read_text(encoding="utf-8", errors="replace")
        )
        if log_likelihood is not None:
            return log_likelihood
    return None


def _existing_iqtree_outputs(
    prefix_path: Path,
    *,
    include_tree: bool = False,
    include_bootstrap: bool = False,
    include_consensus: bool = False,
) -> dict[str, Path]:
    outputs: dict[str, Path] = {}
    for key, candidate in (
        ("iqtree_report", prefix_path.with_suffix(".iqtree")),
        ("iqtree_log", prefix_path.with_suffix(".log")),
    ):
        if candidate.exists():
            outputs[key] = candidate
    tree_candidate = prefix_path.with_suffix(".treefile")
    if include_tree and tree_candidate.exists():
        outputs["tree"] = tree_candidate
    bootstrap_candidate = prefix_path.with_suffix(".ufboot")
    if include_bootstrap and bootstrap_candidate.exists():
        outputs["bootstrap_trees"] = bootstrap_candidate
    consensus_candidate = prefix_path.with_suffix(".contree")
    if include_consensus and consensus_candidate.exists():
        outputs["consensus_tree"] = consensus_candidate
    return outputs


def _build_iqtree_summary(
    prefix_path: Path,
    *,
    default_selected_model: str | None,
    support_tree_path: Path | None = None,
) -> IqtreeWorkflowSummary:
    selected_model = _parse_best_model_artifact(prefix_path) or default_selected_model
    log_likelihood = _parse_log_likelihood_artifact(prefix_path)
    support_values: list[IqtreeSupportValue] = []
    minimum_support: float | None = None
    maximum_support: float | None = None
    if support_tree_path is not None and support_tree_path.exists():
        support_summary = summarize_bootstrap_support_distribution(support_tree_path)
        support_values = [
            IqtreeSupportValue(
                node=node.node,
                descendant_taxa=list(node.descendant_taxa),
                support=node.support,
                support_fraction=node.support_fraction,
                is_backbone=node.is_backbone,
            )
            for node in support_summary.nodes
        ]
        minimum_support = support_summary.minimum_support
        maximum_support = support_summary.maximum_support
    return IqtreeWorkflowSummary(
        selected_model=selected_model,
        log_likelihood=log_likelihood,
        support_value_count=len(support_values),
        minimum_support=minimum_support,
        maximum_support=maximum_support,
        support_values=support_values,
    )


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


def run_codon_aware_multiple_sequence_alignment(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "mafft",
    mode: str = "auto",
    sequence_type: AlignmentAlphabet | None = None,
) -> CodonAwareAlignmentWorkflowReport:
    """Align coding nucleotide sequences through a translated amino-acid guide."""
    prepared_records, preparation = prepare_coding_sequences_for_alignment(
        input_path,
        sequence_type=sequence_type,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    guide_input_path = _sidecar(out_path, "guide-input.fasta")
    guide_alignment_path = _sidecar(out_path, "guide-alignment.fasta")
    exclusion_report_path = _sidecar(out_path, "excluded.tsv")
    guide_records = translate_prepared_coding_sequences(prepared_records)
    write_fasta_alignment(guide_input_path, guide_records)
    mode_args = resolve_mafft_alignment_mode(mode)
    version = read_engine_version("mafft", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    run = execute_engine_command(
        engine_name="mafft",
        workflow="codon-aware-multiple-sequence-alignment",
        executable=resolved,
        version=version,
        command_args=[*mode_args, str(guide_input_path.resolve())],
        work_dir=out_path.parent,
        stdout_path=guide_alignment_path,
        stderr_path=_sidecar(out_path, "stderr.log"),
        output_paths={"guide_alignment": guide_alignment_path},
    )
    aligned_guide = load_fasta_alignment(guide_alignment_path)
    codon_records = back_translate_aligned_coding_sequences(
        aligned_guide,
        coding_records=prepared_records,
    )
    write_fasta_alignment(out_path, codon_records)
    _validate_alignment_output(out_path)
    codon_summary = summarise_fasta(out_path)
    if codon_summary.alignment_length % 3 != 0:
        raise EngineWorkflowError(
            "codon-aware alignment produced an alignment length that is not divisible by three"
        )
    _write_coding_exclusion_table(exclusion_report_path, preparation.excluded_sequences)
    output_paths = {
        "alignment": out_path,
        "guide_input": guide_input_path,
        "guide_alignment": guide_alignment_path,
        "excluded_sequences": exclusion_report_path,
    }
    manifest_path = _manifest_path_from_output(out_path)
    report = CodonAwareAlignmentWorkflowReport(
        workflow="codon-aware-multiple-sequence-alignment",
        engine_name="mafft",
        input_path=input_path,
        output_paths=output_paths,
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        sequence_type=preparation.sequence_type,
        accepted_sequence_count=preparation.accepted_sequence_count,
        excluded_sequences=preparation.excluded_sequences,
        terminal_stop_sequence_count=preparation.terminal_stop_sequence_count,
        notes=[
            "codon-aware alignment preserved nucleotide codon triplets through amino-acid guide alignment",
            f"mafft alignment mode: {mode}",
            f"accepted coding sequences: {preparation.accepted_sequence_count} of {preparation.input_sequence_count}",
            f"retained nucleotide alignment length: {codon_summary.alignment_length}",
        ],
        warnings=list(dict.fromkeys(run.warning_lines + preparation.warnings)),
    )
    report.output_checksums = build_file_checksums(list(output_paths.values()))
    write_engine_manifest(manifest_path, report)
    return report


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
    partition_path: Path | None = None,
    seed: int = 1,
    threads: int = 1,
) -> EngineWorkflowReport:
    """Run a model-selection workflow on an aligned FASTA file."""
    _ensure_inference_ready_alignment(input_path)
    prefix_path = _prefix_path(out_dir, prefix)
    version = read_engine_version("iqtree", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    iqtree_report_path = prefix_path.with_suffix(".iqtree")
    iqtree_log_path = prefix_path.with_suffix(".log")
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="model-selection",
        executable=resolved,
        version=version,
        command_args=[
            *(
                prepared_partitions.command_args
                if prepared_partitions is not None
                else ["-s", str(input_path.resolve())]
            ),
            *(
                []
                if prepared_partitions is not None and prepared_partitions.mixed_data_types
                else _iqtree_sequence_type_flag(input_path, sequence_type)
            ),
            *_iqtree_execution_controls(seed=seed, threads=threads),
            "-m",
            "MF",
            "-pre",
            str(prefix_path.resolve()),
        ],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "iqtree_report": iqtree_report_path,
            "iqtree_log": iqtree_log_path,
        },
    )
    iqtree_summary = _build_iqtree_summary(
        prefix_path,
        default_selected_model=None,
    )
    if iqtree_summary.selected_model is None:
        raise EngineWorkflowError(
            f"iqtree model-selection did not expose a parsable best-fit model in {iqtree_report_path}"
        )
    selected_model = iqtree_summary.selected_model
    selected_model_path = prefix_path.with_suffix(".selected-model.txt")
    selected_model_path.write_text(selected_model + "\n", encoding="utf-8")
    manifest_path = prefix_path.with_suffix(".manifest.json")
    report = EngineWorkflowReport(
        workflow="model-selection",
        engine_name="iqtree",
        input_paths=(
            [input_path]
            if partition_path is None
            else [input_path, partition_path]
        ),
        output_paths={
            **(
                {}
            if prepared_partitions is None
            else prepared_partitions.output_paths
        ),
            **_existing_iqtree_outputs(prefix_path),
            "selected_model": selected_model_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        selected_model=selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        notes=[
            *(
                []
                if prepared_partitions is None
                else prepared_partitions.notes
            ),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            "best-fit substitution model parsed from engine output",
            *(
                []
                if iqtree_summary.log_likelihood is None
                else [
                    "model-selection workflow exposed a parsable log-likelihood score"
                ]
            ),
        ],
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
    partition_path: Path | None = None,
    resume: bool = False,
    seed: int = 1,
    threads: int = 1,
) -> EngineWorkflowReport:
    """Run an external maximum-likelihood tree inference workflow."""
    _ensure_inference_ready_alignment(input_path)
    prefix_path = _prefix_path(out_dir, prefix)
    version = read_engine_version("iqtree", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
    if (
        prepared_partitions is not None
        and not _iqtree_partition_supports_fixed_model(
            model=model,
            mixed_data_types=prepared_partitions.mixed_data_types,
        )
    ):
        raise EngineWorkflowError(
            "mixed DNA/protein partition analyses require a model-selection keyword such as MF, MFP, TEST, or TESTMERGE"
        )
    tree_path = prefix_path.with_suffix(".treefile")
    report_path = prefix_path.with_suffix(".iqtree")
    log_path = prefix_path.with_suffix(".log")
    manifest_path = prefix_path.with_suffix(".manifest.json")
    command = [
        resolved,
        *(
            prepared_partitions.command_args
            if prepared_partitions is not None
            else ["-s", str(input_path.resolve())]
        ),
        *(
            []
            if prepared_partitions is not None and prepared_partitions.mixed_data_types
            else _iqtree_sequence_type_flag(input_path, sequence_type)
        ),
        *_iqtree_execution_controls(seed=seed, threads=threads),
        "-m",
        model,
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=(
                [input_path]
                if partition_path is None
                else [input_path, partition_path]
            ),
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
            "iqtree_log": log_path,
        },
    )
    _validate_tree_output(tree_path)
    iqtree_summary = _build_iqtree_summary(
        prefix_path,
        default_selected_model=model,
        support_tree_path=tree_path,
    )
    report = EngineWorkflowReport(
        workflow="maximum-likelihood-tree",
        engine_name="iqtree",
        input_paths=(
            [input_path]
            if partition_path is None
            else [input_path, partition_path]
        ),
        output_paths={
            **(
                {}
            if prepared_partitions is None
            else prepared_partitions.output_paths
        ),
            **_existing_iqtree_outputs(prefix_path, include_tree=True),
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        selected_model=iqtree_summary.selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        notes=[
            *(
                []
                if prepared_partitions is None
                else prepared_partitions.notes
            ),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            "maximum-likelihood tree validated as parseable Newick output",
            *(
                []
                if iqtree_summary.log_likelihood is None
                else ["log-likelihood parsed from iqtree inference artifacts"]
            ),
            *(
                []
                if iqtree_summary.support_value_count == 0
                else [
                    "support values parsed from the inferred maximum-likelihood tree"
                ]
            ),
        ],
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
    partition_path: Path | None = None,
    resume: bool = False,
    seed: int = 1,
    threads: int = 1,
) -> EngineWorkflowReport:
    """Run external bootstrap support estimation and retain bootstrap trees."""
    if replicates < 1:
        raise ValueError(f"replicates must be positive, got {replicates}")
    _validate_ufboot_replicates(replicates)
    _ensure_inference_ready_alignment(input_path)
    prefix_path = _prefix_path(out_dir, prefix)
    version = read_engine_version("iqtree", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
    if (
        prepared_partitions is not None
        and not _iqtree_partition_supports_fixed_model(
            model=model,
            mixed_data_types=prepared_partitions.mixed_data_types,
        )
    ):
        raise EngineWorkflowError(
            "mixed DNA/protein partition analyses require a model-selection keyword such as MF, MFP, TEST, or TESTMERGE"
        )
    support_tree_path = prefix_path.with_suffix(".treefile")
    bootstrap_tree_path = prefix_path.with_suffix(".ufboot")
    report_path = prefix_path.with_suffix(".iqtree")
    log_path = prefix_path.with_suffix(".log")
    manifest_path = prefix_path.with_suffix(".manifest.json")
    command = [
        resolved,
        *(
            prepared_partitions.command_args
            if prepared_partitions is not None
            else ["-s", str(input_path.resolve())]
        ),
        *(
            []
            if prepared_partitions is not None and prepared_partitions.mixed_data_types
            else _iqtree_sequence_type_flag(input_path, sequence_type)
        ),
        *_iqtree_execution_controls(seed=seed, threads=threads),
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
            input_paths=(
                [input_path]
                if partition_path is None
                else [input_path, partition_path]
            ),
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
            "iqtree_log": log_path,
        },
    )
    _validate_tree_output(support_tree_path)
    if not bootstrap_tree_path.read_text(encoding="utf-8").strip():
        raise EngineWorkflowError(f"bootstrap tree set is empty: {bootstrap_tree_path}")
    iqtree_summary = _build_iqtree_summary(
        prefix_path,
        default_selected_model=model,
        support_tree_path=support_tree_path,
    )
    report = EngineWorkflowReport(
        workflow="bootstrap-support",
        engine_name="iqtree",
        input_paths=(
            [input_path]
            if partition_path is None
            else [input_path, partition_path]
        ),
        output_paths={
            **(
                {}
            if prepared_partitions is None
            else prepared_partitions.output_paths
        ),
            **_existing_iqtree_outputs(
                prefix_path,
                include_tree=False,
                include_bootstrap=True,
                include_consensus=True,
            ),
            "support_tree": support_tree_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        selected_model=iqtree_summary.selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        notes=[
            *(
                []
                if prepared_partitions is None
                else prepared_partitions.notes
            ),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            "bootstrap tree set retained for downstream consensus construction",
            *(
                []
                if iqtree_summary.log_likelihood is None
                else ["log-likelihood parsed from iqtree bootstrap inference artifacts"]
            ),
            *(
                []
                if iqtree_summary.support_value_count == 0
                else [
                    "support values parsed from the bootstrap-supported tree artifact"
                ]
            ),
        ],
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
    log_path = prefix_path.with_suffix(".log")
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
        output_paths={
            "consensus_tree": consensus_tree_path,
            "iqtree_log": log_path,
        },
    )
    _validate_tree_output(consensus_tree_path)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    iqtree_summary = _build_iqtree_summary(
        prefix_path,
        default_selected_model=None,
        support_tree_path=consensus_tree_path,
    )
    report = EngineWorkflowReport(
        workflow="bootstrap-consensus",
        engine_name="iqtree",
        input_paths=[bootstrap_trees_path],
        output_paths=_existing_iqtree_outputs(
            prefix_path, include_consensus=True
        ),
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([bootstrap_trees_path]),
        output_checksums={},
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        notes=[
            "consensus tree validated as parseable Newick output",
            *(
                []
                if iqtree_summary.support_value_count == 0
                else ["support values parsed from the bootstrap consensus tree"]
            ),
        ],
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
