from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import (
    dumps_newick,
    loads_newick,
    write_newick,
    write_newick_tree_set,
)
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.fixed_topology_branch_lengths import (
    optimize_fixed_topology_nucleotide_branch_lengths,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodMultiStartRunSummary,
    NucleotideLikelihoodNniSearchReport,
    NucleotideLikelihoodSprSearchReport,
    NucleotideLikelihoodTbrSearchReport,
    NucleotideLikelihoodTreeInferenceReport,
    SubstitutionModelSelectionReport,
    SubstitutionModelSelectionRow,
)
from bijux_phylogenetics.phylo.likelihood.multi_start_search import (
    build_likelihood_multi_start_candidates,
    rank_likelihood_multi_start_runs,
)
from bijux_phylogenetics.phylo.likelihood.nni_search import (
    search_nucleotide_likelihood_nni,
    write_nucleotide_likelihood_nni_candidate_table,
    write_nucleotide_likelihood_nni_trace_table,
)
from bijux_phylogenetics.phylo.likelihood.search_artifacts import (
    write_nucleotide_likelihood_best_tree_set,
)
from bijux_phylogenetics.phylo.likelihood.spr_search import (
    search_nucleotide_likelihood_spr,
    write_nucleotide_likelihood_spr_trace_table,
)
from bijux_phylogenetics.phylo.likelihood.stepwise_addition import (
    build_likelihood_stepwise_addition_tree,
)
from bijux_phylogenetics.phylo.likelihood.substitution_model_selection import (
    compare_nucleotide_substitution_models,
)
from bijux_phylogenetics.phylo.likelihood.tbr_search import (
    search_nucleotide_likelihood_tbr,
    write_nucleotide_likelihood_tbr_trace_table,
)
from bijux_phylogenetics.phylo.topology import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_SUPPORTED_NATIVE_TREE_INFERENCE_MODEL_NAMES = frozenset(
    {"auto", "jc69", "k80", "f81", "hky85", "gtr"}
)
_SUPPORTED_NATIVE_TREE_INFERENCE_FIXED_MODELS = frozenset(
    {"jc69", "k80", "f81", "hky85", "gtr"}
)
_SUPPORTED_NATIVE_TREE_INFERENCE_MODEL_SELECTION_CRITERIA = frozenset(
    {"aic", "aicc", "bic"}
)
_SUPPORTED_NATIVE_TREE_INFERENCE_SEARCH_METHODS = frozenset({"nni", "spr", "tbr"})
_SUPPORTED_NATIVE_TREE_INFERENCE_MODEL_CANDIDATES = (
    "JC69",
    "K80",
    "F81",
    "HKY85",
    "GTR",
)


@dataclass(frozen=True, slots=True)
class _PreparedTreeInferenceModel:
    model_name: str
    kappa: float | None
    base_frequencies: dict[str, float] | None
    exchangeabilities: dict[str, float] | None


@dataclass(frozen=True, slots=True)
class _TreeInferenceStartCandidate:
    source_kind: str
    source_label: str
    generation_seed: int | None
    tree: PhyloTree


def default_nucleotide_likelihood_tree_inference_model_candidates() -> tuple[str, ...]:
    """Return the native fixed-rate nucleotide candidates supported by ML tree inference."""
    return _SUPPORTED_NATIVE_TREE_INFERENCE_MODEL_CANDIDATES


def validate_nucleotide_likelihood_tree_inference_model_name(model_name: str) -> str:
    """Validate one native ML tree-inference model selection mode."""
    normalized_model_name = model_name.strip().lower()
    if normalized_model_name not in _SUPPORTED_NATIVE_TREE_INFERENCE_MODEL_NAMES:
        raise ValueError(
            "nucleotide likelihood tree inference model_name must be one of "
            + ", ".join(sorted(_SUPPORTED_NATIVE_TREE_INFERENCE_MODEL_NAMES))
        )
    return normalized_model_name


def validate_nucleotide_likelihood_tree_inference_model_selection_criterion(
    criterion: str,
) -> str:
    """Validate one native ML tree-inference model ranking criterion."""
    normalized_criterion = criterion.strip().lower()
    if (
        normalized_criterion
        not in _SUPPORTED_NATIVE_TREE_INFERENCE_MODEL_SELECTION_CRITERIA
    ):
        raise ValueError(
            "nucleotide likelihood tree inference model_selection_criterion must be one of "
            + ", ".join(
                sorted(_SUPPORTED_NATIVE_TREE_INFERENCE_MODEL_SELECTION_CRITERIA)
            )
        )
    return normalized_criterion


def validate_nucleotide_likelihood_tree_inference_search_method(
    search_method: str,
) -> str:
    """Validate one native ML tree-inference local topology search method."""
    normalized_search_method = search_method.strip().lower()
    if normalized_search_method not in _SUPPORTED_NATIVE_TREE_INFERENCE_SEARCH_METHODS:
        raise ValueError(
            "nucleotide likelihood tree inference search_method must be one of "
            + ", ".join(sorted(_SUPPORTED_NATIVE_TREE_INFERENCE_SEARCH_METHODS))
        )
    return normalized_search_method


def validate_nucleotide_likelihood_tree_inference_start_tree_count(
    start_tree_count: int,
) -> int:
    """Validate the requested native ML tree-inference start-tree count."""
    if start_tree_count < 1:
        raise ValueError("start_tree_count must be at least one")
    return start_tree_count


def infer_nucleotide_likelihood_tree(
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str = "auto",
    model_selection_criterion: str = "aic",
    search_method: str = "nni",
    start_tree_count: int = 4,
    start_tree_seed: int = 1,
    branch_reoptimization_policy: str = "coordinate-branch-lengths",
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodTreeInferenceReport:
    """Infer one rooted nucleotide ML tree natively from an alignment."""
    resolved_records, resolved_alignment_path = _resolve_tree_inference_records(records)
    normalized_model_name = validate_nucleotide_likelihood_tree_inference_model_name(
        model_name
    )
    normalized_model_selection_criterion = (
        validate_nucleotide_likelihood_tree_inference_model_selection_criterion(
            model_selection_criterion
        )
    )
    normalized_search_method = (
        validate_nucleotide_likelihood_tree_inference_search_method(search_method)
    )
    validated_start_tree_count = (
        validate_nucleotide_likelihood_tree_inference_start_tree_count(start_tree_count)
    )

    stepwise_tree, stepwise_report = build_likelihood_stepwise_addition_tree(
        resolved_records,
        model_name="jc69",
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    model_selection_report = _build_tree_inference_model_selection_report(
        stepwise_tree,
        resolved_records,
        normalized_model_name=normalized_model_name,
        max_coordinate_passes=max_coordinate_passes,
        improvement_tolerance=improvement_tolerance,
    )
    selected_model_row = _select_tree_inference_model_row(
        model_selection_report,
        normalized_model_name=normalized_model_name,
        model_selection_criterion=normalized_model_selection_criterion,
    )
    selected_model = _prepare_tree_inference_model(selected_model_row)
    optimized_start_tree = _optimize_tree_inference_start_tree(
        stepwise_tree,
        resolved_records,
        selected_model=selected_model,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    start_tree_candidates = _build_tree_inference_start_candidates(
        optimized_start_tree,
        start_tree_count=validated_start_tree_count,
        start_tree_seed=start_tree_seed,
    )
    local_reports = [
        _run_tree_inference_local_search(
            candidate.tree,
            resolved_records,
            selected_model=selected_model,
            search_method=normalized_search_method,
            branch_reoptimization_policy=branch_reoptimization_policy,
            lower_branch_length_bound=lower_branch_length_bound,
            upper_branch_length_bound=upper_branch_length_bound,
            improvement_tolerance=improvement_tolerance,
            max_coordinate_passes=max_coordinate_passes,
        )
        for candidate in start_tree_candidates
    ]
    run_summaries = [
        _build_tree_inference_run_summary(candidate, local_report)
        for candidate, local_report in zip(
            start_tree_candidates,
            local_reports,
            strict=True,
        )
    ]
    ranked_run_indices = rank_likelihood_multi_start_runs(run_summaries)
    for final_likelihood_rank, ranked_run_index in enumerate(
        ranked_run_indices,
        start=1,
    ):
        run_summaries[ranked_run_index].final_likelihood_rank = final_likelihood_rank
    best_run_index = ranked_run_indices[0]
    run_summaries[best_run_index].best_run = True
    best_run = run_summaries[best_run_index]
    best_search_report = local_reports[best_run_index]
    return NucleotideLikelihoodTreeInferenceReport(
        algorithm="nucleotide-likelihood-tree-inference",
        alignment_path=None
        if resolved_alignment_path is None
        else str(resolved_alignment_path),
        taxon_count=len(model_selection_report.taxa),
        site_count=model_selection_report.site_count,
        pattern_count=model_selection_report.pattern_count,
        stepwise_addition_model_name="JC69",
        stepwise_addition_tree_newick=dumps_newick(stepwise_tree),
        stepwise_addition_final_score=stepwise_report.final_score,
        start_tree_source_policy="stepwise-addition-tree-plus-random-tree",
        random_start_tree_count=max(0, validated_start_tree_count - 1),
        start_tree_seed=start_tree_seed,
        model_selection_strategy=(
            "model-selection" if normalized_model_name == "auto" else "fixed-model"
        ),
        model_selection_criterion=(
            normalized_model_selection_criterion
            if normalized_model_name == "auto"
            else None
        ),
        model_selection_tree_newick=model_selection_report.tree_newick,
        selected_model_name=selected_model_row.model_name,
        search_method=normalized_search_method,
        branch_reoptimization_policy=branch_reoptimization_policy,
        run_summaries=run_summaries,
        best_run_source_label=best_run.start_tree_source_label,
        best_final_tree_newick=best_run.final_tree_newick,
        best_final_log_likelihood=best_run.final_log_likelihood,
        best_final_topology_fingerprint=best_run.final_topology_fingerprint,
        model_selection_report=model_selection_report,
        best_search_report=best_search_report,
    )


def infer_nucleotide_likelihood_tree_from_alignment(
    alignment_path: Path,
    *,
    model_name: str = "auto",
    model_selection_criterion: str = "aic",
    search_method: str = "nni",
    start_tree_count: int = 4,
    start_tree_seed: int = 1,
    branch_reoptimization_policy: str = "coordinate-branch-lengths",
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodTreeInferenceReport:
    """Infer one rooted nucleotide ML tree natively from one FASTA alignment path."""
    return infer_nucleotide_likelihood_tree(
        alignment_path,
        model_name=model_name,
        model_selection_criterion=model_selection_criterion,
        search_method=search_method,
        start_tree_count=start_tree_count,
        start_tree_seed=start_tree_seed,
        branch_reoptimization_policy=branch_reoptimization_policy,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def write_nucleotide_likelihood_tree_inference_model_table(
    path: Path,
    report: NucleotideLikelihoodTreeInferenceReport,
) -> Path:
    """Write one deterministic native ML tree-inference model selection ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "model_name",
        "base_model_name",
        "rate_heterogeneity_model",
        "fit_succeeded",
        "parameter_count",
        "log_likelihood",
        "aic",
        "aicc",
        "bic",
        "delta_aic",
        "akaike_weight",
        "rank",
        "selected_by_aic",
        "selected_by_aicc",
        "selected_by_bic",
        "parameter_values",
        "warnings",
    ]
    rows = ["\t".join(columns)]
    for row in report.model_selection_report.rows:
        payload = [
            row.model_name,
            row.base_model_name,
            row.rate_heterogeneity_model,
            row.fit_succeeded,
            row.parameter_count,
            row.log_likelihood,
            row.aic,
            row.aicc,
            row.bic,
            row.delta_aic,
            row.akaike_weight,
            row.rank,
            row.selected_by_aic,
            row.selected_by_aicc,
            row.selected_by_bic,
            json.dumps(row.parameter_values, sort_keys=True),
            "|".join(row.warnings),
        ]
        rows.append(_format_tsv_row(payload))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_tree_inference_likelihood_table(
    path: Path,
    report: NucleotideLikelihoodTreeInferenceReport,
) -> Path:
    """Write one deterministic native ML tree-inference run-likelihood ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "start_tree_source_kind",
        "start_tree_source_label",
        "start_tree_generation_seed",
        "search_algorithm",
        "start_log_likelihood",
        "final_log_likelihood",
        "final_likelihood_rank",
        "final_topology_fingerprint",
        "search_iteration_count",
        "accepted_move_count",
        "evaluated_neighbor_count",
        "branch_reoptimization_policy",
        "substitution_parameter_policy",
        "stopping_reason",
        "best_run",
        "start_tree_newick",
        "final_tree_newick",
    ]
    rows = ["\t".join(columns)]
    for row in report.run_summaries:
        payload = [
            row.start_tree_source_kind,
            row.start_tree_source_label,
            row.start_tree_generation_seed,
            row.search_algorithm,
            row.start_log_likelihood,
            row.final_log_likelihood,
            row.final_likelihood_rank,
            row.final_topology_fingerprint,
            row.search_iteration_count,
            row.accepted_move_count,
            row.evaluated_neighbor_count,
            row.branch_reoptimization_policy,
            row.substitution_parameter_policy,
            row.stopping_reason,
            row.best_run,
            row.start_tree_newick,
            row.final_tree_newick,
        ]
        rows.append(_format_tsv_row(payload))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_tree_inference_run_json(
    path: Path,
    report: NucleotideLikelihoodTreeInferenceReport,
) -> Path:
    """Write one governed native ML tree-inference manifest payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_nucleotide_likelihood_tree_inference_artifacts(
    out_dir: Path,
    report: NucleotideLikelihoodTreeInferenceReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one native ML tree-inference run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    stepwise_tree_path = write_newick(
        out_dir / "stepwise_start_tree.nwk",
        loads_newick(report.stepwise_addition_tree_newick),
    )
    final_tree_path = write_newick(
        out_dir / "final_tree.nwk",
        loads_newick(report.best_final_tree_newick),
    )
    best_tree_set_path = write_nucleotide_likelihood_best_tree_set(
        out_dir / "best_trees.nwk",
        report.best_search_report.equal_best_tree_report,
    )
    start_tree_path = write_newick_tree_set(
        out_dir / "start_trees.nwk",
        [loads_newick(row.start_tree_newick) for row in report.run_summaries],
    )
    likelihood_table_path = write_nucleotide_likelihood_tree_inference_likelihood_table(
        out_dir / "likelihood_table.tsv",
        report,
    )
    model_table_path = write_nucleotide_likelihood_tree_inference_model_table(
        out_dir / "model_table.tsv",
        report,
    )
    search_trace_path = _write_tree_inference_search_trace(
        out_dir / "search_trace.tsv",
        report.best_search_report,
    )
    run_json_path = write_nucleotide_likelihood_tree_inference_run_json(
        out_dir / "run.json",
        report,
    )
    artifact_paths = {
        "stepwise_tree_path": stepwise_tree_path,
        "start_tree_path": start_tree_path,
        "final_tree_path": final_tree_path,
        "best_tree_set_path": best_tree_set_path,
        "likelihood_table_path": likelihood_table_path,
        "model_table_path": model_table_path,
        "search_trace_path": search_trace_path,
        "run_json_path": run_json_path,
    }
    if isinstance(report.best_search_report, NucleotideLikelihoodNniSearchReport):
        artifact_paths["candidate_table_path"] = (
            write_nucleotide_likelihood_nni_candidate_table(
                out_dir / "candidate_table.tsv",
                report.best_search_report,
            )
        )
    return artifact_paths


def _resolve_tree_inference_records(
    records: list[AlignmentRecord] | Path,
) -> tuple[list[AlignmentRecord], Path | None]:
    if isinstance(records, Path):
        return load_fasta_alignment(records), records
    return records, None


def _build_tree_inference_model_selection_report(
    stepwise_tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    normalized_model_name: str,
    max_coordinate_passes: int,
    improvement_tolerance: float,
) -> SubstitutionModelSelectionReport:
    candidate_models = (
        default_nucleotide_likelihood_tree_inference_model_candidates()
        if normalized_model_name == "auto"
        else (normalized_model_name.upper(),)
    )
    return compare_nucleotide_substitution_models(
        stepwise_tree,
        records,
        candidate_models=candidate_models,
        max_coordinate_passes=max_coordinate_passes,
        improvement_tolerance=improvement_tolerance,
    )


def _select_tree_inference_model_row(
    report: SubstitutionModelSelectionReport,
    *,
    normalized_model_name: str,
    model_selection_criterion: str,
) -> SubstitutionModelSelectionRow:
    if normalized_model_name != "auto":
        selected_row = report.rows[0]
        if not selected_row.fit_succeeded:
            raise ValueError(
                f"fixed model {selected_row.model_name} could not be optimized for native tree inference"
            )
        return selected_row
    for row in report.rows:
        if not row.fit_succeeded:
            continue
        if model_selection_criterion == "aic" and row.selected_by_aic:
            return row
        if model_selection_criterion == "aicc" and row.selected_by_aicc:
            return row
        if model_selection_criterion == "bic" and row.selected_by_bic:
            return row
    raise ValueError(
        "native tree inference could not select one comparable substitution model"
    )


def _prepare_tree_inference_model(
    row: SubstitutionModelSelectionRow,
) -> _PreparedTreeInferenceModel:
    base_frequency_keys = ("a", "c", "g", "t")
    base_frequencies: dict[str, float] | None = None
    if any(
        f"base_frequency_{state}" in row.parameter_values
        for state in base_frequency_keys
    ):
        base_frequencies = {
            state.upper(): row.parameter_values[f"base_frequency_{state}"]
            for state in base_frequency_keys
        }
    exchangeabilities = {
        key.removeprefix("exchangeability_").upper(): value
        for key, value in row.parameter_values.items()
        if key.startswith("exchangeability_")
    }
    return _PreparedTreeInferenceModel(
        model_name=row.model_name,
        kappa=row.parameter_values.get("kappa"),
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities or None,
    )


def _optimize_tree_inference_start_tree(
    stepwise_tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    selected_model: _PreparedTreeInferenceModel,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
    improvement_tolerance: float,
    max_coordinate_passes: int,
) -> PhyloTree:
    optimization_report = optimize_fixed_topology_nucleotide_branch_lengths(
        stepwise_tree,
        records,
        model_name=selected_model.model_name.lower(),
        kappa=selected_model.kappa,
        base_frequencies=selected_model.base_frequencies,
        exchangeabilities=selected_model.exchangeabilities,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    return loads_newick(optimization_report.optimized_tree_newick)


def _build_tree_inference_start_candidates(
    optimized_start_tree: PhyloTree,
    *,
    start_tree_count: int,
    start_tree_seed: int,
) -> list[_TreeInferenceStartCandidate]:
    return [
        _TreeInferenceStartCandidate(
            source_kind=(
                "stepwise-addition-tree"
                if candidate.source_kind == "input-tree"
                else candidate.source_kind
            ),
            source_label=(
                "stepwise-addition-tree"
                if candidate.source_kind == "input-tree"
                else candidate.source_label
            ),
            generation_seed=candidate.generation_seed,
            tree=candidate.tree,
        )
        for candidate in build_likelihood_multi_start_candidates(
            optimized_start_tree,
            start_tree_count=start_tree_count,
            start_tree_source_policy="input-tree-plus-random-tree",
            start_tree_seed=start_tree_seed,
        )
    ]


def _run_tree_inference_local_search(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    selected_model: _PreparedTreeInferenceModel,
    search_method: str,
    branch_reoptimization_policy: str,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
    improvement_tolerance: float,
    max_coordinate_passes: int,
) -> (
    NucleotideLikelihoodNniSearchReport
    | NucleotideLikelihoodSprSearchReport
    | NucleotideLikelihoodTbrSearchReport
):
    if search_method == "nni":
        return search_nucleotide_likelihood_nni(
            tree,
            records,
            model_name=selected_model.model_name.lower(),
            branch_reoptimization_policy=branch_reoptimization_policy,
            kappa=selected_model.kappa,
            base_frequencies=selected_model.base_frequencies,
            exchangeabilities=selected_model.exchangeabilities,
            lower_branch_length_bound=lower_branch_length_bound,
            upper_branch_length_bound=upper_branch_length_bound,
            improvement_tolerance=improvement_tolerance,
            max_coordinate_passes=max_coordinate_passes,
        )
    if search_method == "spr":
        return search_nucleotide_likelihood_spr(
            tree,
            records,
            model_name=selected_model.model_name.lower(),
            branch_reoptimization_policy=branch_reoptimization_policy,
            kappa=selected_model.kappa,
            base_frequencies=selected_model.base_frequencies,
            exchangeabilities=selected_model.exchangeabilities,
            lower_branch_length_bound=lower_branch_length_bound,
            upper_branch_length_bound=upper_branch_length_bound,
            improvement_tolerance=improvement_tolerance,
            max_coordinate_passes=max_coordinate_passes,
        )
    return search_nucleotide_likelihood_tbr(
        tree,
        records,
        model_name=selected_model.model_name.lower(),
        branch_reoptimization_policy=branch_reoptimization_policy,
        kappa=selected_model.kappa,
        base_frequencies=selected_model.base_frequencies,
        exchangeabilities=selected_model.exchangeabilities,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def _build_tree_inference_run_summary(
    candidate: _TreeInferenceStartCandidate,
    local_report: (
        NucleotideLikelihoodNniSearchReport
        | NucleotideLikelihoodSprSearchReport
        | NucleotideLikelihoodTbrSearchReport
    ),
) -> NucleotideLikelihoodMultiStartRunSummary:
    return NucleotideLikelihoodMultiStartRunSummary(
        search_algorithm=local_report.algorithm,
        start_tree_source_kind=candidate.source_kind,
        start_tree_source_label=candidate.source_label,
        start_tree_generation_seed=candidate.generation_seed,
        start_tree_newick=local_report.start_tree_newick,
        start_log_likelihood=local_report.start_log_likelihood,
        final_tree_newick=local_report.final_tree_newick,
        final_log_likelihood=local_report.final_log_likelihood,
        final_topology_fingerprint=rooted_topology_fingerprint(
            loads_newick(local_report.final_tree_newick)
        ),
        search_iteration_count=_resolve_tree_inference_iteration_count(local_report),
        accepted_move_count=local_report.accepted_move_count,
        evaluated_neighbor_count=local_report.evaluated_neighbor_count,
        final_likelihood_rank=0,
        branch_reoptimization_policy=local_report.branch_reoptimization_policy,
        substitution_parameter_policy=local_report.substitution_parameter_policy,
        substitution_parameter_values=dict(local_report.substitution_parameter_values),
        substitution_parameter_warnings=list(
            local_report.substitution_parameter_warnings
        ),
        total_branch_optimization_pass_count=local_report.total_branch_optimization_pass_count,
        total_branch_function_evaluation_count=local_report.total_branch_function_evaluation_count,
        stopping_reason=local_report.stopping_reason,
        best_run=False,
    )


def _resolve_tree_inference_iteration_count(
    local_report: (
        NucleotideLikelihoodNniSearchReport
        | NucleotideLikelihoodSprSearchReport
        | NucleotideLikelihoodTbrSearchReport
    ),
) -> int:
    if hasattr(local_report, "iteration_count"):
        return int(local_report.iteration_count)
    return max((row.iteration for row in local_report.trace_rows), default=0)


def _write_tree_inference_search_trace(
    path: Path,
    report: (
        NucleotideLikelihoodNniSearchReport
        | NucleotideLikelihoodSprSearchReport
        | NucleotideLikelihoodTbrSearchReport
    ),
) -> Path:
    if isinstance(report, NucleotideLikelihoodNniSearchReport):
        return write_nucleotide_likelihood_nni_trace_table(path, report)
    if isinstance(report, NucleotideLikelihoodSprSearchReport):
        return write_nucleotide_likelihood_spr_trace_table(path, report)
    return write_nucleotide_likelihood_tbr_trace_table(path, report)


def _format_tsv_row(values: list[object]) -> str:
    return "\t".join(
        format(value, ".15g") if isinstance(value, float) else str(value)
        for value in values
    )
