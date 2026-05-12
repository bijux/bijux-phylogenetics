from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from Bio import Phylo

from bijux_phylogenetics.bayesian.diagnostics import (
    TraceConvergenceReport,
    summarize_trace_convergence,
)
from bijux_phylogenetics.core.alignment import AlignmentAlphabet
from bijux_phylogenetics.core.partitions import (
    LocusPartition,
    build_partition_summary_report,
    normalize_partition_data_type,
    parse_locus_partitions,
    partition_coordinate_text,
)
from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.engines.common import (
    build_file_checksums,
    execute_engine_command,
    read_engine_version,
    resolve_engine_executable,
)
from bijux_phylogenetics.engines.workflows import (
    EngineWorkflowReport,
    _ensure_inference_ready_alignment,
    _persist_workflow_report,
    _resume_existing_workflow,
)
from bijux_phylogenetics.errors import EngineWorkflowError
from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.fasta import infer_alignment_alphabet, load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.tree_set import (
    compute_clade_frequency_table,
    compute_consensus_tree,
    load_tree_set,
)


@dataclass(slots=True)
class MrBayesPreparationReport:
    alignment_path: Path
    nexus_path: Path
    partition_path: Path | None
    taxon_count: int
    character_count: int
    inferred_alphabet: AlignmentAlphabet
    partition_count: int
    partition_names: list[str]
    partition_data_types: list[str]
    partition_warnings: list[str]
    model: str
    rates: str
    ngen: int
    nchains: int
    samplefreq: int
    printfreq: int
    burnin_fraction: float


@dataclass(slots=True)
class MrBayesTraceRow:
    generation: int
    values: dict[str, float]


@dataclass(slots=True)
class MrBayesTraceReport:
    path: Path
    row_count: int
    columns: list[str]
    rows: list[MrBayesTraceRow]


@dataclass(slots=True)
class EffectiveSampleSize:
    parameter: str
    sample_count: int
    effective_sample_size: float


@dataclass(slots=True)
class MrBayesESSReport:
    path: Path
    sample_count: int
    effective_sample_sizes: list[EffectiveSampleSize]


@dataclass(slots=True)
class MrBayesConvergenceReport:
    path: Path
    sample_count: int
    converged: bool
    ess_threshold: float
    mean_shift_threshold: float
    warnings: list[dict[str, object]]
    parameter_summaries: list[dict[str, object]]


@dataclass(slots=True)
class MrBayesPosteriorSummaryReport:
    source_path: Path
    filtered_tree_set_path: Path
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    rooted_topology_count: int
    shared_taxa: list[str]
    consensus_newick: str
    clade_frequency_count: int


def _mrbayes_datatype(alphabet: AlignmentAlphabet) -> str:
    if alphabet in {"dna", "rna"}:
        return "dna"
    if alphabet == "protein":
        return "protein"
    raise EngineWorkflowError(
        f"MrBayes preparation requires a recognized alignment alphabet, got {alphabet}"
    )


def _mrbayes_model_commands(
    *,
    alphabet: AlignmentAlphabet,
    model: str,
    rates: str,
    partition_count: int,
) -> list[str]:
    normalized_model = model.lower()
    normalized_rates = rates.lower()
    applyto_prefix = " applyto=(all)" if partition_count > 1 else ""
    if alphabet in {"dna", "rna"}:
        if normalized_model == "jc69":
            base = [f"lset{applyto_prefix} nst=1 rates=equal;"]
        elif normalized_model == "hky":
            base = [f"lset{applyto_prefix} nst=2 rates=equal;"]
        elif normalized_model == "gtr":
            base = [f"lset{applyto_prefix} nst=6 rates=equal;"]
        else:
            raise EngineWorkflowError(f"unsupported nucleotide MrBayes model: {model}")
    else:
        if normalized_model not in {"wag", "jones", "dayhoff", "poisson"}:
            raise EngineWorkflowError(f"unsupported protein MrBayes model: {model}")
        base = [
            f"prset{applyto_prefix} aamodelpr=fixed({normalized_model});",
            f"lset{applyto_prefix} rates=equal;",
        ]

    if normalized_rates == "equal":
        return base
    if normalized_rates in {"gamma", "invgamma", "propinv"}:
        updated: list[str] = []
        for line in base:
            if line.startswith("lset "):
                updated.append(line.replace("rates=equal", f"rates={normalized_rates}"))
            else:
                updated.append(line)
        return updated
    raise EngineWorkflowError(f"unsupported MrBayes rate model: {rates}")


def _mrbayes_partition_datatype(alphabet: AlignmentAlphabet) -> str:
    if alphabet in {"dna", "rna"}:
        return "DNA"
    if alphabet == "protein":
        return "PROTEIN"
    raise EngineWorkflowError(
        f"MrBayes preparation requires a recognized alignment alphabet, got {alphabet}"
    )


def _validate_mrbayes_partitions(
    partitions: tuple[LocusPartition, ...],
    *,
    alphabet: AlignmentAlphabet,
    alignment_length: int,
) -> tuple[list[str], list[str]]:
    summary = build_partition_summary_report(
        partitions, alignment_length=alignment_length
    )
    expected_data_type = _mrbayes_partition_datatype(alphabet)
    declared_data_types = [
        normalize_partition_data_type(partition.data_type)
        for partition in partitions
        if partition.data_type is not None
    ]
    mismatched = sorted(
        {
            data_type
            for data_type in declared_data_types
            if data_type is not None and data_type != expected_data_type
        }
    )
    if mismatched:
        raise EngineWorkflowError(
            "MrBayes preparation requires partition datatypes that match the "
            f"alignment alphabet {expected_data_type}, got {', '.join(mismatched)}"
        )
    return summary.declared_data_types, summary.warnings


def _mrbayes_charset_commands(partitions: tuple[LocusPartition, ...]) -> list[str]:
    return [
        f"charset {partition.name} = {partition_coordinate_text(partition)};"
        for partition in partitions
    ]


def _mrbayes_partition_commands(partitions: tuple[LocusPartition, ...]) -> list[str]:
    if len(partitions) <= 1:
        return []
    partition_name = "loci"
    joined_names = ", ".join(partition.name for partition in partitions)
    return [
        f"partition {partition_name} = {len(partitions)}: {joined_names};",
        f"set partition={partition_name};",
        "prset applyto=(all) ratepr=variable;",
    ]


def prepare_mrbayes_analysis(
    alignment_path: Path,
    nexus_path: Path,
    *,
    partition_path: Path | None = None,
    model: str = "gtr",
    rates: str = "gamma",
    ngen: int = 10000,
    nchains: int = 4,
    samplefreq: int = 100,
    printfreq: int = 100,
    burnin_fraction: float = 0.25,
) -> MrBayesPreparationReport:
    """Prepare a MrBayes NEXUS analysis specification from an aligned FASTA file."""
    _ensure_inference_ready_alignment(alignment_path)
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    records = load_fasta_alignment(alignment_path)
    alphabet = infer_alignment_alphabet(records)
    datatype = _mrbayes_datatype(alphabet)
    partitions: tuple[LocusPartition, ...] = ()
    partition_data_types: list[str] = []
    partition_warnings: list[str] = []
    if partition_path is not None:
        partitions = parse_locus_partitions(partition_path)
        partition_data_types, partition_warnings = _validate_mrbayes_partitions(
            partitions,
            alphabet=alphabet,
            alignment_length=len(records[0].sequence),
        )
    model_commands = _mrbayes_model_commands(
        alphabet=alphabet,
        model=model,
        rates=rates,
        partition_count=len(partitions),
    )
    matrix_lines = "\n".join(
        f"{record.identifier} {record.sequence}" for record in records
    )
    command_block = "\n".join(
        [
            "begin mrbayes;",
            "  set autoclose=yes nowarn=yes;",
            *[f"  {line}" for line in _mrbayes_charset_commands(partitions)],
            *[f"  {line}" for line in _mrbayes_partition_commands(partitions)],
            *[f"  {line}" for line in model_commands],
            f"  mcmcp ngen={ngen} nchains={nchains} samplefreq={samplefreq} printfreq={printfreq};",
            "  mcmc;",
            f"  sump burninfrac={burnin_fraction:.6f};",
            f"  sumt burninfrac={burnin_fraction:.6f};",
            "end;",
        ]
    )
    nexus_text = "\n".join(
        [
            "#NEXUS",
            "",
            "begin data;",
            f"  dimensions ntax={len(records)} nchar={len(records[0].sequence)};",
            f"  format datatype={datatype} missing=? gap=-;",
            "  matrix",
            matrix_lines,
            "  ;",
            "end;",
            "",
            command_block,
            "",
        ]
    )
    nexus_path.parent.mkdir(parents=True, exist_ok=True)
    nexus_path.write_text(nexus_text, encoding="utf-8")
    return MrBayesPreparationReport(
        alignment_path=alignment_path,
        nexus_path=nexus_path,
        partition_path=partition_path,
        taxon_count=len(records),
        character_count=len(records[0].sequence),
        inferred_alphabet=alphabet,
        partition_count=max(len(partitions), 1),
        partition_names=[partition.name for partition in partitions],
        partition_data_types=partition_data_types,
        partition_warnings=partition_warnings,
        model=model,
        rates=rates,
        ngen=ngen,
        nchains=nchains,
        samplefreq=samplefreq,
        printfreq=printfreq,
        burnin_fraction=burnin_fraction,
    )


def run_mrbayes_posterior_inference(
    nexus_path: Path,
    *,
    executable: str | Path = "mb",
    resume: bool = False,
) -> EngineWorkflowReport:
    """Run a MrBayes posterior tree inference workflow from a prepared NEXUS file."""
    if not nexus_path.exists():
        raise FileNotFoundError(nexus_path)
    version = read_engine_version("MrBayes", executable, version_args=("--version",))
    resolved = resolve_engine_executable(executable)
    prefix_path = nexus_path.with_suffix("")
    trace_path = Path(f"{nexus_path}.run1.p")
    tree_path = Path(f"{nexus_path}.run1.t")
    manifest_path = prefix_path.with_suffix(".manifest.json")
    command = [resolved, nexus_path.name]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[nexus_path],
            expected_command=command,
        )
        if resumed is not None:
            return resumed
    run = execute_engine_command(
        engine_name="MrBayes",
        workflow="posterior-tree-inference",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=nexus_path.parent,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "posterior_trees": tree_path,
            "parameter_traces": trace_path,
        },
    )
    parse_mrbayes_parameter_traces(trace_path)
    summarize_mrbayes_posterior_trees(tree_path, burnin_fraction=0.25)
    report = EngineWorkflowReport(
        workflow="posterior-tree-inference",
        engine_name="MrBayes",
        input_paths=[nexus_path],
        output_paths={
            "posterior_trees": tree_path,
            "parameter_traces": trace_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([nexus_path]),
        output_checksums={},
        notes=[
            "MrBayes posterior trees and parameter traces validated after engine execution"
        ],
    )
    return _persist_workflow_report(report)


def parse_mrbayes_parameter_traces(path: Path) -> MrBayesTraceReport:
    """Parse a MrBayes parameter trace table into deterministic numeric rows."""
    rows: list[MrBayesTraceRow] = []
    with path.open(encoding="utf-8", newline="") as handle:
        filtered_lines = [
            line
            for line in handle
            if line.strip() and not line.lstrip().startswith("[")
        ]
    reader = csv.DictReader(filtered_lines, delimiter="\t")
    if reader.fieldnames is None:
        raise EngineWorkflowError(f"MrBayes trace file contains no header: {path}")
    columns = [field for field in reader.fieldnames if field and field != "Gen"]
    for raw_row in reader:
        generation_text = raw_row.get("Gen") or raw_row.get("gen")
        if generation_text is None:
            raise EngineWorkflowError(f"MrBayes trace file lacks a Gen column: {path}")
        values = {
            column: float(raw_row[column])
            for column in columns
            if raw_row.get(column) not in {None, ""}
        }
        rows.append(
            MrBayesTraceRow(generation=int(float(generation_text)), values=values)
        )
    if not rows:
        raise EngineWorkflowError(
            f"MrBayes trace file contains no sampled rows: {path}"
        )
    return MrBayesTraceReport(
        path=path, row_count=len(rows), columns=columns, rows=rows
    )


def compute_mrbayes_effective_sample_sizes(path: Path) -> MrBayesESSReport:
    """Compute per-parameter effective sample sizes from a MrBayes trace file."""
    report = parse_mrbayes_parameter_traces(path)
    convergence = summarize_trace_convergence(
        path=path,
        rows=[row.values for row in report.rows],
        columns=report.columns,
    )
    effective_sample_sizes = [
        EffectiveSampleSize(
            parameter=summary.parameter,
            sample_count=summary.sample_count,
            effective_sample_size=summary.effective_sample_size,
        )
        for summary in convergence.series
    ]
    return MrBayesESSReport(
        path=path,
        sample_count=report.row_count,
        effective_sample_sizes=effective_sample_sizes,
    )


def assess_mrbayes_convergence(
    path: Path,
    *,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
) -> MrBayesConvergenceReport:
    """Flag low-ESS or unstable MrBayes trace parameters."""
    report = parse_mrbayes_parameter_traces(path)
    convergence = summarize_trace_convergence(
        path=path,
        rows=[row.values for row in report.rows],
        columns=report.columns,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
    )
    return _build_mrbayes_convergence_report(convergence)


def _build_mrbayes_convergence_report(
    convergence: TraceConvergenceReport,
) -> MrBayesConvergenceReport:
    return MrBayesConvergenceReport(
        path=convergence.path,
        sample_count=convergence.sample_count,
        converged=convergence.converged,
        ess_threshold=convergence.ess_threshold,
        mean_shift_threshold=convergence.mean_shift_threshold,
        warnings=[
            {
                "parameter": warning.parameter,
                "code": warning.code,
                "message": warning.message,
                "observed_value": warning.observed_value,
                "threshold": warning.threshold,
            }
            for warning in convergence.warnings
        ],
        parameter_summaries=[
            {
                "parameter": summary.parameter,
                "sample_count": summary.sample_count,
                "effective_sample_size": summary.effective_sample_size,
                "mean": summary.mean,
                "minimum": summary.minimum,
                "maximum": summary.maximum,
                "first_half_mean": summary.first_half_mean,
                "second_half_mean": summary.second_half_mean,
                "standardized_mean_shift": summary.standardized_mean_shift,
            }
            for summary in convergence.series
        ],
    )


def summarize_mrbayes_posterior_trees(
    tree_set_path: Path,
    *,
    burnin_fraction: float = 0.25,
) -> tuple[PhyloTree, MrBayesPosteriorSummaryReport]:
    """Summarize MrBayes posterior trees after discarding a burn-in fraction."""
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    tree_format = detect_tree_format(tree_set_path)
    bio_trees = list(Phylo.parse(tree_set_path, tree_format))
    if not bio_trees:
        raise EngineWorkflowError(
            f"MrBayes posterior tree file contains no trees: {tree_set_path}"
        )
    burnin_tree_count = int(len(bio_trees) * burnin_fraction)
    kept_trees = bio_trees[burnin_tree_count:]
    if not kept_trees:
        raise EngineWorkflowError(
            f"MrBayes posterior tree file is empty after burn-in filtering: {tree_set_path}"
        )
    filtered_tree_set_path = tree_set_path.with_suffix(".postburnin.nwk")
    filtered_tree_set_path.write_text(
        "".join(
            dumps_newick(tree_from_biophylo(tree, source_format="newick")) + "\n"
            for tree in kept_trees
        ),
        encoding="utf-8",
    )
    summary = load_tree_set(filtered_tree_set_path)
    consensus_tree, consensus = compute_consensus_tree(filtered_tree_set_path)
    clade_frequencies = compute_clade_frequency_table(filtered_tree_set_path)
    return consensus_tree, MrBayesPosteriorSummaryReport(
        source_path=tree_set_path,
        filtered_tree_set_path=filtered_tree_set_path,
        total_tree_count=len(bio_trees),
        burnin_tree_count=burnin_tree_count,
        kept_tree_count=len(kept_trees),
        rooted_topology_count=summary.rooted_topology_count,
        shared_taxa=summary.shared_taxa,
        consensus_newick=consensus.consensus_newick,
        clade_frequency_count=len(clade_frequencies.clade_frequencies),
    )
