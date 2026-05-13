from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re

from bijux_phylogenetics.bayesian.burnin import (
    DEFAULT_BURNIN_FRACTIONS,
    BurninSensitivityCladeShift,
    BurninSensitivityParameterShift,
    normalize_burnin_fractions,
    summarize_burnin_clade_shifts,
    summarize_burnin_parameter_shifts,
)
from bijux_phylogenetics.bayesian.diagnostics import (
    TraceConvergenceReport,
    summarize_trace_convergence,
    summarize_trace_parameters,
)
from bijux_phylogenetics.core.alignment import AlignmentAlphabet
from bijux_phylogenetics.core.metadata import write_taxon_rows
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
    validate_timeout_seconds,
)
from bijux_phylogenetics.engines.workflows import (
    EngineWorkflowReport,
    _ensure_inference_ready_alignment,
    _persist_workflow_report,
    _resolve_incomplete_workflow_state,
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
_MRBAYES_PROBABILITY_PERCENT_PATTERN = re.compile(r'prob\(percent\)="([0-9.eE+-]+)"')


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
class MrBayesParameterSummary:
    parameter: str
    sample_count: int
    effective_sample_size: float
    mean: float
    median: float
    standard_deviation: float
    minimum: float
    maximum: float
    hpd_95_lower: float
    hpd_95_upper: float
    first_half_mean: float
    second_half_mean: float
    standardized_mean_shift: float


@dataclass(slots=True)
class MrBayesParameterDiagnosticsReport:
    path: Path
    burnin_fraction: float
    burnin_row_count: int
    kept_row_count: int
    first_kept_generation: int
    last_kept_generation: int
    parameter_summaries: list[MrBayesParameterSummary]


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
class MrBayesBurninSensitivitySlice:
    burnin_fraction: float
    burnin_tree_count: int
    kept_tree_count: int
    rooted_topology_count: int
    clade_frequency_count: int
    consensus_newick: str
    kept_row_count: int | None
    first_kept_generation: int | None
    last_kept_generation: int | None
    lnl_mean: float | None
    tree_length_mean: float | None


@dataclass(slots=True)
class MrBayesBurninSensitivityReport:
    posterior_tree_path: Path
    trace_path: Path | None
    slices: list[MrBayesBurninSensitivitySlice]
    changed_consensus_count: int
    parameter_shifts: list[BurninSensitivityParameterShift]
    clade_shifts: list[BurninSensitivityCladeShift]
    unstable_parameter_count: int
    unstable_clade_count: int
    warnings: list[str]


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
        raise EngineWorkflowError(
            "MrBayes tree file has an unterminated translate block"
        )
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
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run a MrBayes posterior tree inference workflow from a prepared NEXUS file."""
    if not nexus_path.exists():
        raise FileNotFoundError(nexus_path)
    validate_timeout_seconds(timeout_seconds)
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
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    version = read_engine_version(
        "MrBayes",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
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
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
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
            "MrBayes posterior trees, parameter traces, consensus tree, and MCMC diagnostics validated after engine execution",
            *incomplete_notes,
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


def parse_mrbayes_consensus_tree(
    path: Path,
) -> tuple[PhyloTree, MrBayesConsensusTreeReport]:
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
        float(match.group(1))
        for match in _MRBAYES_PROBABILITY_PATTERN.finditer(tree_text)
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


def summarize_mrbayes_parameter_diagnostics(
    path: Path,
    *,
    burnin_fraction: float = 0.0,
) -> MrBayesParameterDiagnosticsReport:
    """Summarize burn-in-aware posterior parameter diagnostics from MrBayes traces."""
    report = parse_mrbayes_parameter_traces(path)
    burnin_row_count, kept_rows = _split_mrbayes_trace_rows(
        report, burnin_fraction=burnin_fraction
    )
    diagnostics = summarize_trace_parameters(
        path=path,
        rows=[row.values for row in kept_rows],
        columns=report.columns,
    )
    return MrBayesParameterDiagnosticsReport(
        path=path,
        burnin_fraction=burnin_fraction,
        burnin_row_count=burnin_row_count,
        kept_row_count=len(kept_rows),
        first_kept_generation=kept_rows[0].generation,
        last_kept_generation=kept_rows[-1].generation,
        parameter_summaries=[
            MrBayesParameterSummary(
                parameter=summary.parameter,
                sample_count=summary.sample_count,
                effective_sample_size=summary.effective_sample_size,
                mean=summary.mean,
                median=summary.median,
                standard_deviation=summary.standard_deviation,
                minimum=summary.minimum,
                maximum=summary.maximum,
                hpd_95_lower=summary.hpd_95_lower,
                hpd_95_upper=summary.hpd_95_upper,
                first_half_mean=summary.first_half_mean,
                second_half_mean=summary.second_half_mean,
                standardized_mean_shift=summary.standardized_mean_shift,
            )
            for summary in diagnostics.series
        ],
    )


def write_mrbayes_parameter_summary_table(
    path: Path,
    report: MrBayesParameterDiagnosticsReport,
) -> Path:
    """Write a reviewer-facing TSV summary of MrBayes posterior parameter diagnostics."""
    return write_taxon_rows(
        path,
        columns=[
            "parameter",
            "sample_count",
            "effective_sample_size",
            "mean",
            "median",
            "standard_deviation",
            "minimum",
            "maximum",
            "hpd_95_lower",
            "hpd_95_upper",
            "first_half_mean",
            "second_half_mean",
            "standardized_mean_shift",
            "burnin_fraction",
            "burnin_row_count",
            "kept_row_count",
            "first_kept_generation",
            "last_kept_generation",
        ],
        rows=[
            {
                "parameter": summary.parameter,
                "sample_count": str(summary.sample_count),
                "effective_sample_size": format(summary.effective_sample_size, ".15g"),
                "mean": format(summary.mean, ".15g"),
                "median": format(summary.median, ".15g"),
                "standard_deviation": format(summary.standard_deviation, ".15g"),
                "minimum": format(summary.minimum, ".15g"),
                "maximum": format(summary.maximum, ".15g"),
                "hpd_95_lower": format(summary.hpd_95_lower, ".15g"),
                "hpd_95_upper": format(summary.hpd_95_upper, ".15g"),
                "first_half_mean": format(summary.first_half_mean, ".15g"),
                "second_half_mean": format(summary.second_half_mean, ".15g"),
                "standardized_mean_shift": format(
                    summary.standardized_mean_shift, ".15g"
                ),
                "burnin_fraction": format(report.burnin_fraction, ".15g"),
                "burnin_row_count": str(report.burnin_row_count),
                "kept_row_count": str(report.kept_row_count),
                "first_kept_generation": str(report.first_kept_generation),
                "last_kept_generation": str(report.last_kept_generation),
            }
            for summary in report.parameter_summaries
        ],
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
                "median": summary.median,
                "standard_deviation": summary.standard_deviation,
                "minimum": summary.minimum,
                "maximum": summary.maximum,
                "hpd_95_lower": summary.hpd_95_lower,
                "hpd_95_upper": summary.hpd_95_upper,
                "first_half_mean": summary.first_half_mean,
                "second_half_mean": summary.second_half_mean,
                "standardized_mean_shift": summary.standardized_mean_shift,
            }
            for summary in convergence.series
        ],
    )


def _split_mrbayes_trace_rows(
    report: MrBayesTraceReport,
    *,
    burnin_fraction: float,
) -> tuple[int, list[MrBayesTraceRow]]:
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    burnin_row_count = int(report.row_count * burnin_fraction)
    kept_rows = report.rows[burnin_row_count:]
    if not kept_rows:
        raise ValueError(
            "burnin_fraction discards every MrBayes trace row; reduce the burn-in"
        )
    return burnin_row_count, kept_rows


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


def assess_mrbayes_burnin_sensitivity(
    posterior_tree_path: Path,
    *,
    trace_path: Path | None = None,
    burnin_fractions: tuple[float, ...] = DEFAULT_BURNIN_FRACTIONS,
) -> MrBayesBurninSensitivityReport:
    """Compare MrBayes posterior summaries across multiple burn-in fractions."""
    ordered_fractions = normalize_burnin_fractions(burnin_fractions)
    slices: list[MrBayesBurninSensitivitySlice] = []
    previous_consensus: str | None = None
    changed_consensus_count = 0
    parameter_summaries_by_fraction: dict[float, list[MrBayesParameterSummary]] = {}
    clade_frequencies_by_fraction: dict[float, list[object]] = {}
    for fraction in ordered_fractions:
        _, posterior_summary = summarize_mrbayes_posterior_trees(
            posterior_tree_path,
            burnin_fraction=fraction,
        )
        clade_report = compute_clade_frequency_table(
            posterior_summary.filtered_tree_set_path
        )
        kept_row_count = None
        first_kept_generation = None
        last_kept_generation = None
        lnl_mean = None
        tree_length_mean = None
        if trace_path is not None:
            trace_summary = summarize_mrbayes_parameter_diagnostics(
                trace_path,
                burnin_fraction=fraction,
            )
            parameter_summaries_by_fraction[fraction] = (
                trace_summary.parameter_summaries
            )
            kept_row_count = trace_summary.kept_row_count
            first_kept_generation = trace_summary.first_kept_generation
            last_kept_generation = trace_summary.last_kept_generation
            lnl_mean = _mean_mrbayes_parameter(trace_summary, "LnL")
            tree_length_mean = _mean_mrbayes_parameter(trace_summary, "TL")
            if tree_length_mean is None:
                tree_length_mean = _mean_mrbayes_parameter(trace_summary, "TL{all}")
        clade_frequencies_by_fraction[fraction] = list(clade_report.clade_frequencies)
        slices.append(
            MrBayesBurninSensitivitySlice(
                burnin_fraction=fraction,
                burnin_tree_count=posterior_summary.burnin_tree_count,
                kept_tree_count=posterior_summary.kept_tree_count,
                rooted_topology_count=posterior_summary.rooted_topology_count,
                clade_frequency_count=posterior_summary.clade_frequency_count,
                consensus_newick=posterior_summary.consensus_newick,
                kept_row_count=kept_row_count,
                first_kept_generation=first_kept_generation,
                last_kept_generation=last_kept_generation,
                lnl_mean=lnl_mean,
                tree_length_mean=tree_length_mean,
            )
        )
        if (
            previous_consensus is not None
            and previous_consensus != posterior_summary.consensus_newick
        ):
            changed_consensus_count += 1
        previous_consensus = posterior_summary.consensus_newick
    parameter_shifts = summarize_burnin_parameter_shifts(
        parameter_summaries_by_fraction
    )
    clade_shifts = summarize_burnin_clade_shifts(clade_frequencies_by_fraction)
    warnings: list[str] = []
    if changed_consensus_count:
        warnings.append(
            "majority-rule consensus topology changes across tested burn-in fractions"
        )
    if any(shift.unstable for shift in parameter_shifts):
        warnings.append(
            "one or more posterior parameter 95% HPD intervals do not overlap across tested burn-in fractions"
        )
    if any(shift.unstable for shift in clade_shifts):
        warnings.append(
            "one or more posterior clade probabilities cross the majority-rule threshold across tested burn-in fractions"
        )
    return MrBayesBurninSensitivityReport(
        posterior_tree_path=posterior_tree_path,
        trace_path=trace_path,
        slices=slices,
        changed_consensus_count=changed_consensus_count,
        parameter_shifts=parameter_shifts,
        clade_shifts=clade_shifts,
        unstable_parameter_count=sum(1 for shift in parameter_shifts if shift.unstable),
        unstable_clade_count=sum(1 for shift in clade_shifts if shift.unstable),
        warnings=warnings,
    )


def write_mrbayes_burnin_sensitivity_slice_table(
    path: Path,
    report: MrBayesBurninSensitivityReport,
) -> Path:
    """Write one row per tested MrBayes burn-in fraction."""
    return write_taxon_rows(
        path,
        columns=[
            "burnin_fraction",
            "burnin_tree_count",
            "kept_tree_count",
            "rooted_topology_count",
            "clade_frequency_count",
            "kept_row_count",
            "first_kept_generation",
            "last_kept_generation",
            "lnl_mean",
            "tree_length_mean",
            "consensus_newick",
        ],
        rows=[
            {
                "burnin_fraction": format(row.burnin_fraction, ".15g"),
                "burnin_tree_count": str(row.burnin_tree_count),
                "kept_tree_count": str(row.kept_tree_count),
                "rooted_topology_count": str(row.rooted_topology_count),
                "clade_frequency_count": str(row.clade_frequency_count),
                "kept_row_count": ""
                if row.kept_row_count is None
                else str(row.kept_row_count),
                "first_kept_generation": ""
                if row.first_kept_generation is None
                else str(row.first_kept_generation),
                "last_kept_generation": ""
                if row.last_kept_generation is None
                else str(row.last_kept_generation),
                "lnl_mean": ""
                if row.lnl_mean is None
                else format(row.lnl_mean, ".15g"),
                "tree_length_mean": ""
                if row.tree_length_mean is None
                else format(row.tree_length_mean, ".15g"),
                "consensus_newick": row.consensus_newick,
            }
            for row in report.slices
        ],
    )


def _mean_mrbayes_parameter(
    report: MrBayesParameterDiagnosticsReport,
    parameter: str,
) -> float | None:
    for summary in report.parameter_summaries:
        if summary.parameter == parameter:
            return summary.mean
    return None
