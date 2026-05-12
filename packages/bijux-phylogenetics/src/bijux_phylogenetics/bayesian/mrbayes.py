from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re

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
from bijux_phylogenetics.io.biopython import loads_biophylo
from bijux_phylogenetics.io.fasta import infer_alignment_alphabet, load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.tree_set import (
    compute_clade_frequency_table,
    compute_consensus_tree,
    load_tree_set,
)

_MRBAYES_TREE_PATTERN = re.compile(
    r"tree\s+([^\s=]+)\s*=\s*(.+?);", flags=re.IGNORECASE | re.DOTALL
)
_MRBAYES_TREE_GENERATION_PATTERN = re.compile(r"(\d+)$")
_MRBAYES_PROBABILITY_PATTERN = re.compile(r"prob=([0-9.eE+-]+)")
_MRBAYES_PROBABILITY_PERCENT_PATTERN = re.compile(
    r'prob\(percent\)="([0-9.eE+-]+)"'
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


@dataclass(slots=True)
class MrBayesPosteriorTreeSample:
    tree_name: str
    generation: int | None
    rooted: bool | None
    tip_names: list[str]
    newick: str


@dataclass(slots=True)
class MrBayesPosteriorTreeSetReport:
    path: Path
    tree_count: int
    rooted_tree_count: int
    sampled_generations: list[int]
    tip_names: list[str]
    trees: list[MrBayesPosteriorTreeSample]


@dataclass(slots=True)
class MrBayesMcmcRow:
    generation: int
    values: dict[str, float | None]


@dataclass(slots=True)
class MrBayesMcmcReport:
    path: Path
    row_count: int
    columns: list[str]
    comment_lines: list[str]
    rows: list[MrBayesMcmcRow]


@dataclass(slots=True)
class MrBayesConsensusTreeReport:
    path: Path
    tree_name: str
    rooted: bool | None
    tip_names: list[str]
    consensus_newick: str
    annotated_node_count: int
    minimum_posterior_probability: float | None
    maximum_posterior_probability: float | None
    minimum_posterior_probability_percent: float | None
    maximum_posterior_probability_percent: float | None


def _extract_mrbayes_tree_entries(text: str) -> list[tuple[str, str]]:
    entries = [
        (match.group(1), match.group(2).strip())
        for match in _MRBAYES_TREE_PATTERN.finditer(text)
    ]
    if not entries:
        raise EngineWorkflowError("MrBayes tree file contains no tree entries")
    return entries


def _split_nexus_translate_entries(raw_block: str) -> list[str]:
    entries: list[str] = []
    current: list[str] = []
    in_single_quote = False
    for character in raw_block:
        if character == "'":
            in_single_quote = not in_single_quote
        if character == "," and not in_single_quote:
            candidate = "".join(current).strip()
            if candidate:
                entries.append(candidate)
            current = []
            continue
        current.append(character)
    tail = "".join(current).strip()
    if tail:
        entries.append(tail)
    return entries


def _parse_nexus_translate_map(text: str) -> dict[str, str]:
    lowered = text.lower()
    marker = "translate"
    start = lowered.find(marker)
    if start == -1:
        return {}
    remainder = text[start + len(marker) :]
    end = remainder.find(";")
    if end == -1:
        raise EngineWorkflowError("MrBayes tree file has an unterminated translate block")
    block = remainder[:end]
    mapping: dict[str, str] = {}
    for entry in _split_nexus_translate_entries(block):
        parts = entry.split(None, 1)
        if len(parts) != 2:
            continue
        key, value = parts
        mapping[key.strip()] = value.strip().strip("'")
    return mapping


def _strip_square_bracket_comments(text: str) -> str:
    stripped: list[str] = []
    depth = 0
    for character in text:
        if character == "[":
            depth += 1
            continue
        if character == "]" and depth:
            depth -= 1
            continue
        if depth == 0:
            stripped.append(character)
    return "".join(stripped)


def _translate_mrbayes_tip_labels(newick: str, mapping: dict[str, str]) -> str:
    if not mapping:
        return newick

    def replace(match: re.Match[str]) -> str:
        token = match.group(1)
        translated = mapping.get(token, token)
        return match.group(0).replace(token, translated)

    return re.sub(r"(?<=[(,])\s*([A-Za-z0-9_.-]+)(?=\s*[:),])", replace, newick)


def _detect_mrbayes_rooted_flag(tree_text: str) -> bool | None:
    prefix = tree_text.lstrip()
    if prefix.startswith("[&R]"):
        return True
    if prefix.startswith("[&U]"):
        return False
    return None


def _parse_mrbayes_tree_generation(tree_name: str) -> int | None:
    match = _MRBAYES_TREE_GENERATION_PATTERN.search(tree_name)
    return None if match is None else int(match.group(1))


def _parse_mrbayes_tree_text(
    tree_text: str, *, translation: dict[str, str]
) -> tuple[str, PhyloTree, bool | None]:
    rooted = _detect_mrbayes_rooted_flag(tree_text)
    stripped = _strip_square_bracket_comments(tree_text).strip()
    translated = _translate_mrbayes_tip_labels(stripped, translation)
    tree = loads_biophylo(f"{translated};", source_format="newick")
    return dumps_newick(tree), tree, rooted


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
    mcmc_path = Path(f"{nexus_path}.mcmc")
    consensus_path = Path(f"{nexus_path}.con.tre")
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
            "mcmc_diagnostics": mcmc_path,
            "consensus_tree": consensus_path,
        },
    )
    parse_mrbayes_parameter_traces(trace_path)
    parse_mrbayes_mcmc_diagnostics(mcmc_path)
    parse_mrbayes_consensus_tree(consensus_path)
    summarize_mrbayes_posterior_trees(tree_path, burnin_fraction=0.25)
    report = EngineWorkflowReport(
        workflow="posterior-tree-inference",
        engine_name="MrBayes",
        input_paths=[nexus_path],
        output_paths={
            "posterior_trees": tree_path,
            "parameter_traces": trace_path,
            "mcmc_diagnostics": mcmc_path,
            "consensus_tree": consensus_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([nexus_path]),
        output_checksums={},
        notes=[
            "MrBayes posterior trees, parameter traces, consensus tree, and MCMC diagnostics validated after engine execution"
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


def parse_mrbayes_mcmc_diagnostics(path: Path) -> MrBayesMcmcReport:
    """Parse a MrBayes .mcmc diagnostics table into deterministic rows."""
    comment_lines: list[str] = []
    table_lines: list[str] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("["):
                comment_lines.append(stripped)
                continue
            table_lines.append(line)
    reader = csv.DictReader(table_lines, delimiter="\t")
    if reader.fieldnames is None:
        raise EngineWorkflowError(
            f"MrBayes MCMC diagnostics file contains no header: {path}"
        )
    columns = [field for field in reader.fieldnames if field and field != "Gen"]
    rows: list[MrBayesMcmcRow] = []
    for raw_row in reader:
        generation_text = raw_row.get("Gen") or raw_row.get("gen")
        if generation_text is None:
            raise EngineWorkflowError(
                f"MrBayes MCMC diagnostics file lacks a Gen column: {path}"
            )
        values: dict[str, float | None] = {}
        for column in columns:
            raw_value = raw_row.get(column)
            if raw_value in {None, ""}:
                continue
            normalized = raw_value.strip()
            if normalized.lower() in {"na", "nan"}:
                values[column] = None
            else:
                values[column] = float(normalized)
        rows.append(
            MrBayesMcmcRow(generation=int(float(generation_text)), values=values)
        )
    if not rows:
        raise EngineWorkflowError(
            f"MrBayes MCMC diagnostics file contains no sampled rows: {path}"
        )
    return MrBayesMcmcReport(
        path=path,
        row_count=len(rows),
        columns=columns,
        comment_lines=comment_lines,
        rows=rows,
    )


def parse_mrbayes_posterior_tree_samples(path: Path) -> MrBayesPosteriorTreeSetReport:
    """Parse a MrBayes posterior tree set into generation-tagged samples."""
    text = path.read_text(encoding="utf-8")
    translation = _parse_nexus_translate_map(text)
    samples: list[MrBayesPosteriorTreeSample] = []
    for tree_name, tree_text in _extract_mrbayes_tree_entries(text):
        newick, tree, rooted = _parse_mrbayes_tree_text(
            tree_text, translation=translation
        )
        samples.append(
            MrBayesPosteriorTreeSample(
                tree_name=tree_name,
                generation=_parse_mrbayes_tree_generation(tree_name),
                rooted=rooted if rooted is not None else tree.rooted,
                tip_names=tree.tip_names,
                newick=newick,
            )
        )
    if not samples:
        raise EngineWorkflowError(
            f"MrBayes posterior tree file contains no trees: {path}"
        )
    rooted_tree_count = sum(1 for sample in samples if sample.rooted)
    sampled_generations = [
        generation
        for generation in (sample.generation for sample in samples)
        if generation is not None
    ]
    return MrBayesPosteriorTreeSetReport(
        path=path,
        tree_count=len(samples),
        rooted_tree_count=rooted_tree_count,
        sampled_generations=sampled_generations,
        tip_names=samples[0].tip_names,
        trees=samples,
    )


def parse_mrbayes_consensus_tree(path: Path) -> tuple[PhyloTree, MrBayesConsensusTreeReport]:
    """Parse a MrBayes consensus tree with posterior-probability annotations."""
    text = path.read_text(encoding="utf-8")
    translation = _parse_nexus_translate_map(text)
    entries = _extract_mrbayes_tree_entries(text)
    if len(entries) != 1:
        raise EngineWorkflowError(
            f"MrBayes consensus tree file must contain exactly one tree: {path}"
        )
    tree_name, tree_text = entries[0]
    consensus_newick, tree, rooted = _parse_mrbayes_tree_text(
        tree_text, translation=translation
    )
    posterior_probabilities = [
        float(match.group(1)) for match in _MRBAYES_PROBABILITY_PATTERN.finditer(tree_text)
    ]
    posterior_probability_percents = [
        float(match.group(1))
        for match in _MRBAYES_PROBABILITY_PERCENT_PATTERN.finditer(tree_text)
    ]
    report = MrBayesConsensusTreeReport(
        path=path,
        tree_name=tree_name,
        rooted=rooted if rooted is not None else tree.rooted,
        tip_names=tree.tip_names,
        consensus_newick=consensus_newick,
        annotated_node_count=len(posterior_probabilities),
        minimum_posterior_probability=(
            None if not posterior_probabilities else min(posterior_probabilities)
        ),
        maximum_posterior_probability=(
            None if not posterior_probabilities else max(posterior_probabilities)
        ),
        minimum_posterior_probability_percent=(
            None
            if not posterior_probability_percents
            else min(posterior_probability_percents)
        ),
        maximum_posterior_probability_percent=(
            None
            if not posterior_probability_percents
            else max(posterior_probability_percents)
        ),
    )
    return tree, report


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
    tree_set_report = parse_mrbayes_posterior_tree_samples(tree_set_path)
    burnin_tree_count = int(tree_set_report.tree_count * burnin_fraction)
    kept_trees = tree_set_report.trees[burnin_tree_count:]
    if not kept_trees:
        raise EngineWorkflowError(
            f"MrBayes posterior tree file is empty after burn-in filtering: {tree_set_path}"
        )
    filtered_tree_set_path = tree_set_path.with_suffix(".postburnin.nwk")
    filtered_tree_set_path.write_text(
        "".join(f"{sample.newick}\n" for sample in kept_trees),
        encoding="utf-8",
    )
    summary = load_tree_set(filtered_tree_set_path)
    consensus_tree, consensus = compute_consensus_tree(filtered_tree_set_path)
    clade_frequencies = compute_clade_frequency_table(filtered_tree_set_path)
    return consensus_tree, MrBayesPosteriorSummaryReport(
        source_path=tree_set_path,
        filtered_tree_set_path=filtered_tree_set_path,
        total_tree_count=tree_set_report.tree_count,
        burnin_tree_count=burnin_tree_count,
        kept_tree_count=len(kept_trees),
        rooted_topology_count=summary.rooted_topology_count,
        shared_taxa=summary.shared_taxa,
        consensus_newick=consensus.consensus_newick,
        clade_frequency_count=len(clade_frequencies.clade_frequencies),
    )
