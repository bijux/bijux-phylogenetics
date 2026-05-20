from __future__ import annotations

from pathlib import Path
import re

from bijux_phylogenetics.bayesian.burnin import (
    DEFAULT_BURNIN_FRACTIONS,
    normalize_burnin_fractions,
    summarize_burnin_clade_shifts,
    summarize_burnin_parameter_shifts,
)
from bijux_phylogenetics.bayesian.diagnostics import (
    TraceConvergenceReport,
    summarize_trace_convergence,
    summarize_trace_parameters,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.engines.common import (
    build_file_checksums,
    execute_engine_command,
    read_engine_version,
    resolve_engine_executable,
    validate_timeout_seconds,
)
from bijux_phylogenetics.engines.workflows.models import EngineWorkflowReport
from bijux_phylogenetics.engines.workflows.state import (
    _persist_workflow_report,
    _record_output_validation_failure,
    _resolve_incomplete_workflow_state,
    _resume_existing_workflow,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError, PhylogeneticsError
from bijux_phylogenetics.io.biopython import loads_biophylo
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.trees import (
    compute_clade_frequency_table,
    compute_consensus_tree,
    load_tree_set,
)
from .models import (
    EffectiveSampleSize,
    MrBayesBurninSensitivityReport,
    MrBayesBurninSensitivitySlice,
    MrBayesConsensusTreeReport,
    MrBayesConvergenceReport,
    MrBayesESSReport,
    MrBayesMcmcReport,
    MrBayesMcmcRow,
    MrBayesParameterDiagnosticsReport,
    MrBayesParameterSummary,
    MrBayesPosteriorSummaryReport,
    MrBayesPosteriorTreeSample,
    MrBayesPosteriorTreeSetReport,
    MrBayesPreparationReport,
    MrBayesTraceReport,
    MrBayesTraceRow,
)
from .preparation import prepare_mrbayes_analysis
from .tabular import (
    parse_mrbayes_mcmc_diagnostics,
    parse_mrbayes_parameter_traces,
)

_MRBAYES_TREE_PATTERN = re.compile(
    r"tree\s+([^\s=]+)\s*=\s*(.+?);", flags=re.IGNORECASE | re.DOTALL
)
_MRBAYES_TREE_GENERATION_PATTERN = re.compile(r"(\d+)$")
_MRBAYES_PROBABILITY_PATTERN = re.compile(r"prob=([0-9.eE+-]+)")
_MRBAYES_PROBABILITY_PERCENT_PATTERN = re.compile(r'prob\(percent\)="([0-9.eE+-]+)"')



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
        raise _mrbayes_artifact_error(
            f"MrBayes analysis NEXUS file was not found: {nexus_path}",
            code="mrbayes_analysis_missing_file",
            path=nexus_path,
            artifact_kind="mrbayes-analysis-nexus",
            details={"expected_section": "analysis nexus file"},
        )
    validate_timeout_seconds(timeout_seconds)
    resolved = resolve_engine_executable(executable)
    prefix_path = nexus_path.with_suffix("")
    trace_path = Path(f"{nexus_path}.run1.p")
    tree_path = Path(f"{nexus_path}.run1.t")
    mcmc_path = Path(f"{nexus_path}.mcmc")
    consensus_path = Path(f"{nexus_path}.con.tre")
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "MrBayes",
        executable,
        version_args=("-v",),
        timeout_seconds=timeout_seconds,
    )
    command = [resolved, nexus_path.name]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[nexus_path],
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None:
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
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
    try:
        parse_mrbayes_parameter_traces(trace_path)
        parse_mrbayes_mcmc_diagnostics(mcmc_path)
        parse_mrbayes_consensus_tree(consensus_path)
        summarize_mrbayes_posterior_trees(tree_path, burnin_fraction=0.25)
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
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
        config={
            "timeout_seconds": timeout_seconds,
        },
        notes=[
            "MrBayes posterior trees, parameter traces, consensus tree, and MCMC diagnostics validated after engine execution",
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


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
