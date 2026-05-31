from __future__ import annotations

from pathlib import Path

# Bandit flags subprocess imports broadly; this benchmark intentionally shells
# out so it can measure external-process behavior in a governed way.
import subprocess  # nosec B404

from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.comparative.evolutionary_modes import (
    transform_tree_for_evolutionary_mode,
)
from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.simulation import (
    simulate_brownian_traits,
    simulate_discrete_traits,
    write_continuous_trait_table,
    write_discrete_trait_table,
)

from ..geiger_reference import (
    GEIGER_LARGE_TREE_MODEL_FITTING_REFERENCE_PAYLOADS,
)
from .case_definitions import (
    ContinuousCaseDefinition,
    DiscreteCaseDefinition,
)
from .contracts import (
    LargeTreeModelFittingObservation,
    LargeTreeModelFittingThreshold,
)
from .measurement import (
    measure,
    measure_continuous_fit,
    optional_float,
    optional_int,
    optional_string,
)


def run_case(
    case_definition: ContinuousCaseDefinition | DiscreteCaseDefinition,
    working_root: Path,
) -> LargeTreeModelFittingObservation:
    if isinstance(case_definition, ContinuousCaseDefinition):
        return run_continuous_case(case_definition, working_root)
    return run_discrete_case(case_definition, working_root)


def run_continuous_case(
    case_definition: ContinuousCaseDefinition,
    working_root: Path,
) -> LargeTreeModelFittingObservation:
    tree_path = write_balanced_benchmark_tree(
        working_root / "benchmark-tree.nwk",
        benchmark_taxa(case_definition.taxon_count),
    )
    if case_definition.fit_mode == "pagel-lambda":
        transformed_tree = transform_tree_for_evolutionary_mode(
            load_tree(tree_path),
            mode="pagel-lambda",
            parameter_value=case_definition.lambda_value or 0.0,
        )
        transformed_tree_path = working_root / "simulation-tree.nwk"
        transformed_tree_path.write_text(
            dumps_newick(transformed_tree) + "\n",
            encoding="utf-8",
        )
        simulation = simulate_brownian_traits(
            transformed_tree_path,
            root_state=case_definition.root_state,
            sigma=case_definition.sigma,
            seed=case_definition.seed,
        )
    else:
        simulation = simulate_brownian_traits(
            tree_path,
            root_state=case_definition.root_state,
            sigma=case_definition.sigma,
            seed=case_definition.seed,
        )
    traits_path = write_continuous_trait_table(working_root / "traits.tsv", simulation)

    if case_definition.fit_mode == "pagel-lambda":
        notes = [
            "simulates Brownian trait values on a lambda-transformed rooted tree and refits Pagel-lambda on the original branch-length surface",
            "uses governed reduced profile-search budgets so the benchmark measures large-tree throughput rather than full dense parity search",
        ]
    else:
        notes = [
            "simulates Brownian trait values on a rooted 512-taxon balanced tree and refits the Brownian continuous likelihood surface directly on that large tree",
            "the heavy tier uses a parameter-free continuous fit so the repository still proves 500-plus-taxon continuous fitting even when transformed optimizers remain materially slower",
        ]
    try:
        fit_report, runtime_seconds, peak_memory_bytes = measure_continuous_fit(
            case_definition=case_definition,
            tree_path=tree_path,
            traits_path=traits_path,
        )
    except (TimeoutError, subprocess.TimeoutExpired):
        geiger_timeout_fields = continuous_geiger_fields(
            GEIGER_LARGE_TREE_MODEL_FITTING_REFERENCE_PAYLOADS.get(
                case_definition.case_id
            )
        )
        return LargeTreeModelFittingObservation(
            case_id=case_definition.case_id,
            tier=case_definition.tier,
            trait_kind="continuous",
            fit_surface=case_definition.fit_surface,
            taxon_count=case_definition.taxon_count,
            status="timeout-review",
            runtime_seconds=case_definition.timeout_seconds,
            peak_memory_bytes=None,
            memory_observation_kind=None,
            optimizer_name=None,
            optimizer_iteration_count=None,
            optimizer_function_evaluation_count=None,
            converged=None,
            hit_lower_parameter_boundary=None,
            hit_upper_parameter_boundary=None,
            unstable_review=False,
            too_slow_review=True,
            stable_conclusion_supported=None,
            threshold=case_definition.threshold,
            runtime_within_threshold=False,
            peak_memory_within_threshold=None,
            optimizer_step_within_threshold=None,
            performance_threshold_passed=False,
            geiger_reference_available=(
                case_definition.case_id
                in GEIGER_LARGE_TREE_MODEL_FITTING_REFERENCE_PAYLOADS
            ),
            geiger_runtime_seconds=geiger_timeout_fields[0],
            geiger_optimizer_step_count=geiger_timeout_fields[1],
            geiger_parameter_name=geiger_timeout_fields[2],
            geiger_parameter_value=geiger_timeout_fields[3],
            geiger_rate=geiger_timeout_fields[4],
            geiger_log_likelihood=geiger_timeout_fields[5],
            geiger_aic=geiger_timeout_fields[6],
            geiger_aicc=geiger_timeout_fields[7],
            parameter_delta=None,
            rate_delta=None,
            log_likelihood_delta=None,
            aic_delta=None,
            geiger_match_tolerance=case_definition.geiger_match_tolerance,
            matches_geiger_reference=None,
            notes=[
                *notes,
                "owned fit exceeded the governed heavy-tier runtime budget before returning a continuous-model result",
            ],
        )
    except Exception as error:
        return failed_observation(
            case_id=case_definition.case_id,
            tier=case_definition.tier,
            trait_kind="continuous",
            fit_surface=case_definition.fit_surface,
            taxon_count=case_definition.taxon_count,
            threshold=case_definition.threshold,
            geiger_match_tolerance=case_definition.geiger_match_tolerance,
            notes=[*notes, f"owned fit raised {type(error).__name__}: {error}"],
        )

    optimizer_diagnostics = fit_report.optimizer_diagnostics
    optimizer_step_count = (
        None
        if optimizer_diagnostics is None
        else optimizer_diagnostics.function_evaluation_count
    )
    geiger_reference = GEIGER_LARGE_TREE_MODEL_FITTING_REFERENCE_PAYLOADS.get(
        case_definition.case_id
    )
    (
        geiger_runtime_seconds,
        geiger_optimizer_step_count,
        geiger_parameter_name,
        geiger_parameter_value,
        geiger_rate,
        geiger_log_likelihood,
        geiger_aic,
        geiger_aicc,
    ) = continuous_geiger_fields(geiger_reference)
    parameter_delta = optional_delta(fit_report.parameter_value, geiger_parameter_value)
    rate_delta = optional_delta(fit_report.rate, geiger_rate)
    log_likelihood_delta = optional_delta(
        fit_report.log_likelihood,
        geiger_log_likelihood,
    )
    aic_delta = optional_delta(fit_report.aic, geiger_aic)
    matches_geiger_reference = match_summary(
        [
            parameter_delta,
            rate_delta,
            log_likelihood_delta,
            aic_delta,
        ],
        tolerance=case_definition.geiger_match_tolerance,
        reference_available=geiger_reference is not None,
    )
    (
        runtime_within_threshold,
        peak_memory_within_threshold,
        optimizer_step_within_threshold,
        performance_threshold_passed,
    ) = evaluate_threshold(
        threshold=case_definition.threshold,
        runtime_seconds=runtime_seconds,
        peak_memory_bytes=peak_memory_bytes,
        optimizer_step_count=optimizer_step_count,
    )
    boundary_assessment = fit_report.boundary_assessment
    unstable_review = (
        boundary_assessment is not None
        and boundary_assessment.stable_conclusion_supported is False
    ) or (optimizer_diagnostics is not None and not optimizer_diagnostics.converged)
    too_slow_review = performance_threshold_passed is False and (
        runtime_within_threshold is False
    )
    if unstable_review:
        notes.append(
            "owned fit triggered identifiability or boundary review on the large-tree surface"
        )
    if too_slow_review:
        notes.append(
            "owned fit completed but exceeded the governed runtime budget for this benchmark case"
        )
    return LargeTreeModelFittingObservation(
        case_id=case_definition.case_id,
        tier=case_definition.tier,
        trait_kind="continuous",
        fit_surface=case_definition.fit_surface,
        taxon_count=case_definition.taxon_count,
        status="ok",
        runtime_seconds=runtime_seconds,
        peak_memory_bytes=peak_memory_bytes,
        memory_observation_kind="python-tracemalloc",
        optimizer_name=(
            None
            if optimizer_diagnostics is None
            else optimizer_diagnostics.optimizer_name
        ),
        optimizer_iteration_count=optimizer_step_count,
        optimizer_function_evaluation_count=optimizer_step_count,
        converged=(
            None if optimizer_diagnostics is None else optimizer_diagnostics.converged
        ),
        hit_lower_parameter_boundary=(
            None
            if optimizer_diagnostics is None
            else optimizer_diagnostics.hit_lower_boundary
        ),
        hit_upper_parameter_boundary=(
            None
            if optimizer_diagnostics is None
            else optimizer_diagnostics.hit_upper_boundary
        ),
        unstable_review=unstable_review,
        too_slow_review=too_slow_review,
        stable_conclusion_supported=(
            None
            if boundary_assessment is None
            else boundary_assessment.stable_conclusion_supported
        ),
        threshold=case_definition.threshold,
        runtime_within_threshold=runtime_within_threshold,
        peak_memory_within_threshold=peak_memory_within_threshold,
        optimizer_step_within_threshold=optimizer_step_within_threshold,
        performance_threshold_passed=performance_threshold_passed,
        geiger_reference_available=geiger_reference is not None,
        geiger_runtime_seconds=geiger_runtime_seconds,
        geiger_optimizer_step_count=geiger_optimizer_step_count,
        geiger_parameter_name=geiger_parameter_name,
        geiger_parameter_value=geiger_parameter_value,
        geiger_rate=geiger_rate,
        geiger_log_likelihood=geiger_log_likelihood,
        geiger_aic=geiger_aic,
        geiger_aicc=geiger_aicc,
        parameter_delta=parameter_delta,
        rate_delta=rate_delta,
        log_likelihood_delta=log_likelihood_delta,
        aic_delta=aic_delta,
        geiger_match_tolerance=case_definition.geiger_match_tolerance,
        matches_geiger_reference=matches_geiger_reference,
        notes=notes,
    )


def run_discrete_case(
    case_definition: DiscreteCaseDefinition,
    working_root: Path,
) -> LargeTreeModelFittingObservation:
    tree_path = write_balanced_benchmark_tree(
        working_root / "benchmark-tree.nwk",
        benchmark_taxa(case_definition.taxon_count),
    )
    simulation = simulate_discrete_traits(
        tree_path,
        states=list(case_definition.states),
        transition_rate=case_definition.transition_rate,
        root_state=case_definition.root_state,
        seed=case_definition.seed,
    )
    traits_path = write_discrete_trait_table(working_root / "traits.tsv", simulation)
    notes = [
        "simulates a binary equal-rates discrete trait on a rooted balanced tree and refits the ER likelihood surface on the same 100-taxon tree",
        "records the owned discrete optimizer iteration and function-evaluation counts directly from the fit report",
    ]
    try:
        fit_report, runtime_seconds, peak_memory_bytes = measure(
            lambda: fit_discrete_mk_model(
                tree_path,
                traits_path,
                trait="state",
                model="equal-rates",
            )
        )
    except Exception as error:
        return failed_observation(
            case_id=case_definition.case_id,
            tier=case_definition.tier,
            trait_kind="discrete",
            fit_surface=case_definition.fit_surface,
            taxon_count=case_definition.taxon_count,
            threshold=case_definition.threshold,
            geiger_match_tolerance=case_definition.geiger_match_tolerance,
            notes=[*notes, f"owned fit raised {type(error).__name__}: {error}"],
        )

    optimizer_diagnostics = fit_report.optimizer_diagnostics
    optimizer_step_count = optimizer_diagnostics.iteration_count
    geiger_reference = GEIGER_LARGE_TREE_MODEL_FITTING_REFERENCE_PAYLOADS.get(
        case_definition.case_id
    )
    (
        geiger_runtime_seconds,
        geiger_optimizer_step_count,
        geiger_parameter_name,
        geiger_parameter_value,
        geiger_rate,
        geiger_log_likelihood,
        geiger_aic,
        geiger_aicc,
    ) = discrete_geiger_fields(geiger_reference)
    representative_rate = representative_discrete_rate(fit_report)
    parameter_delta = optional_delta(None, geiger_parameter_value)
    rate_delta = optional_delta(representative_rate, geiger_rate)
    log_likelihood_delta = optional_delta(
        fit_report.log_likelihood,
        geiger_log_likelihood,
    )
    aic_delta = optional_delta(fit_report.aic, geiger_aic)
    matches_geiger_reference = match_summary(
        [
            rate_delta,
            log_likelihood_delta,
            aic_delta,
        ],
        tolerance=case_definition.geiger_match_tolerance,
        reference_available=geiger_reference is not None,
    )
    (
        runtime_within_threshold,
        peak_memory_within_threshold,
        optimizer_step_within_threshold,
        performance_threshold_passed,
    ) = evaluate_threshold(
        threshold=case_definition.threshold,
        runtime_seconds=runtime_seconds,
        peak_memory_bytes=peak_memory_bytes,
        optimizer_step_count=optimizer_step_count,
    )
    unstable_review = (
        fit_report.overparameterized
        or not optimizer_diagnostics.converged
        or optimizer_diagnostics.hit_lower_parameter_bound
        or optimizer_diagnostics.hit_upper_parameter_bound
    )
    too_slow_review = performance_threshold_passed is False and (
        runtime_within_threshold is False
    )
    if unstable_review:
        notes.append(
            "owned fit triggered optimizer-boundary or overparameterization review on the large discrete surface"
        )
    if too_slow_review:
        notes.append(
            "owned fit completed but exceeded the governed runtime budget for this benchmark case"
        )
    return LargeTreeModelFittingObservation(
        case_id=case_definition.case_id,
        tier=case_definition.tier,
        trait_kind="discrete",
        fit_surface=case_definition.fit_surface,
        taxon_count=case_definition.taxon_count,
        status="ok",
        runtime_seconds=runtime_seconds,
        peak_memory_bytes=peak_memory_bytes,
        memory_observation_kind="python-tracemalloc",
        optimizer_name=optimizer_diagnostics.optimizer_name,
        optimizer_iteration_count=optimizer_diagnostics.iteration_count,
        optimizer_function_evaluation_count=(
            optimizer_diagnostics.function_evaluation_count
        ),
        converged=optimizer_diagnostics.converged,
        hit_lower_parameter_boundary=optimizer_diagnostics.hit_lower_parameter_bound,
        hit_upper_parameter_boundary=optimizer_diagnostics.hit_upper_parameter_bound,
        unstable_review=unstable_review,
        too_slow_review=too_slow_review,
        stable_conclusion_supported=not unstable_review,
        threshold=case_definition.threshold,
        runtime_within_threshold=runtime_within_threshold,
        peak_memory_within_threshold=peak_memory_within_threshold,
        optimizer_step_within_threshold=optimizer_step_within_threshold,
        performance_threshold_passed=performance_threshold_passed,
        geiger_reference_available=geiger_reference is not None,
        geiger_runtime_seconds=geiger_runtime_seconds,
        geiger_optimizer_step_count=geiger_optimizer_step_count,
        geiger_parameter_name=geiger_parameter_name,
        geiger_parameter_value=geiger_parameter_value,
        geiger_rate=geiger_rate,
        geiger_log_likelihood=geiger_log_likelihood,
        geiger_aic=geiger_aic,
        geiger_aicc=geiger_aicc,
        parameter_delta=parameter_delta,
        rate_delta=rate_delta,
        log_likelihood_delta=log_likelihood_delta,
        aic_delta=aic_delta,
        geiger_match_tolerance=case_definition.geiger_match_tolerance,
        matches_geiger_reference=matches_geiger_reference,
        notes=notes,
    )


def failed_observation(
    *,
    case_id: str,
    tier: str,
    trait_kind: str,
    fit_surface: str,
    taxon_count: int,
    threshold: LargeTreeModelFittingThreshold,
    geiger_match_tolerance: float,
    notes: list[str],
) -> LargeTreeModelFittingObservation:
    return LargeTreeModelFittingObservation(
        case_id=case_id,
        tier=tier,
        trait_kind=trait_kind,
        fit_surface=fit_surface,
        taxon_count=taxon_count,
        status="failed",
        runtime_seconds=None,
        peak_memory_bytes=None,
        memory_observation_kind=None,
        optimizer_name=None,
        optimizer_iteration_count=None,
        optimizer_function_evaluation_count=None,
        converged=None,
        hit_lower_parameter_boundary=None,
        hit_upper_parameter_boundary=None,
        unstable_review=True,
        too_slow_review=False,
        stable_conclusion_supported=None,
        threshold=threshold,
        runtime_within_threshold=None,
        peak_memory_within_threshold=None,
        optimizer_step_within_threshold=None,
        performance_threshold_passed=None,
        geiger_reference_available=False,
        geiger_runtime_seconds=None,
        geiger_optimizer_step_count=None,
        geiger_parameter_name=None,
        geiger_parameter_value=None,
        geiger_rate=None,
        geiger_log_likelihood=None,
        geiger_aic=None,
        geiger_aicc=None,
        parameter_delta=None,
        rate_delta=None,
        log_likelihood_delta=None,
        aic_delta=None,
        geiger_match_tolerance=geiger_match_tolerance,
        matches_geiger_reference=None,
        notes=notes,
    )


def evaluate_threshold(
    *,
    threshold: LargeTreeModelFittingThreshold,
    runtime_seconds: float | None,
    peak_memory_bytes: int | None,
    optimizer_step_count: int | None,
) -> tuple[bool | None, bool | None, bool | None, bool | None]:
    runtime_within_threshold = (
        None
        if runtime_seconds is None
        else runtime_seconds <= threshold.max_runtime_seconds
    )
    peak_memory_within_threshold = (
        None
        if peak_memory_bytes is None
        else peak_memory_bytes <= threshold.max_peak_memory_bytes
    )
    optimizer_step_within_threshold = (
        True
        if optimizer_step_count is None and threshold.max_optimizer_step_count == 0
        else (
            None
            if optimizer_step_count is None
            else optimizer_step_count <= threshold.max_optimizer_step_count
        )
    )
    values = (
        runtime_within_threshold,
        peak_memory_within_threshold,
        optimizer_step_within_threshold,
    )
    if any(value is None for value in values):
        return (
            runtime_within_threshold,
            peak_memory_within_threshold,
            optimizer_step_within_threshold,
            None,
        )
    return (
        runtime_within_threshold,
        peak_memory_within_threshold,
        optimizer_step_within_threshold,
        all(value is True for value in values),
    )


def continuous_geiger_fields(
    payload: dict[str, object] | None,
) -> tuple[
    float | None,
    int | None,
    str | None,
    float | None,
    float | None,
    float | None,
    float | None,
    float | None,
]:
    if payload is None:
        return (None, None, None, None, None, None, None, None)
    summary = payload.get("fit_summary")
    if not isinstance(summary, dict):
        return (None, None, None, None, None, None, None, None)
    optimizer_result = summary.get("optimizer_result")
    optimizer_step_count = None
    if isinstance(optimizer_result, dict):
        optimizer_step_count = optional_int(
            optimizer_result.get("attempt_count")
            or optimizer_result.get("converged_attempt_count")
        )
    return (
        optional_float(summary.get("runtime_seconds")),
        optimizer_step_count,
        optional_string(summary.get("parameter_name")),
        optional_float(summary.get("parameter_value")),
        optional_float(summary.get("rate")),
        optional_float(summary.get("log_likelihood")),
        optional_float(summary.get("aic")),
        optional_float(summary.get("aicc")),
    )


def discrete_geiger_fields(
    payload: dict[str, object] | None,
) -> tuple[
    float | None,
    int | None,
    str | None,
    float | None,
    float | None,
    float | None,
    float | None,
    float | None,
]:
    if payload is None:
        return (None, None, None, None, None, None, None, None)
    summary = payload.get("fit_summary")
    if not isinstance(summary, dict):
        return (None, None, None, None, None, None, None, None)
    optimizer_result = summary.get("optimizer_result")
    optimizer_step_count = None
    if isinstance(optimizer_result, dict):
        optimizer_step_count = optional_int(
            optimizer_result.get("attempt_count")
            or optimizer_result.get("converged_attempt_count")
        )
    return (
        optional_float(summary.get("runtime_seconds")),
        optimizer_step_count,
        optional_string(summary.get("parameter_name")),
        optional_float(summary.get("parameter_value")),
        optional_float(summary.get("representative_rate")),
        optional_float(summary.get("log_likelihood")),
        optional_float(summary.get("aic")),
        optional_float(summary.get("aicc")),
    )


def representative_discrete_rate(fit_report) -> float | None:
    allowed_rows = [
        row.rate for row in fit_report.transition_rate_rows if row.transition_allowed
    ]
    if not allowed_rows:
        return None
    return allowed_rows[0]


def match_summary(
    deltas: list[float | None],
    *,
    tolerance: float,
    reference_available: bool,
) -> bool | None:
    if not reference_available:
        return None
    comparable = [delta for delta in deltas if delta is not None]
    if not comparable:
        return False
    return all(delta <= tolerance for delta in comparable)


def optional_delta(left: float | None, right: float | None) -> float | None:
    if left is None or right is None:
        return None
    return abs(left - right)


def benchmark_taxa(count: int) -> list[str]:
    return [f"taxon_{index:04d}" for index in range(1, count + 1)]


def write_balanced_benchmark_tree(path: Path, taxa: list[str]) -> Path:
    leaves = [TreeNode(name=taxon, branch_length=0.1) for taxon in taxa]
    while len(leaves) > 1:
        next_level: list[TreeNode] = []
        for index in range(0, len(leaves), 2):
            left = leaves[index]
            right = leaves[index + 1] if index + 1 < len(leaves) else None
            if right is None:
                left.branch_length = round((left.branch_length or 0.0) + 0.1, 15)
                next_level.append(left)
                continue
            next_level.append(TreeNode(children=[left, right], branch_length=0.1))
        leaves = next_level
    root = leaves[0]
    root.branch_length = None
    return write_newick(path, PhyloTree(root=root, source_format="newick"))
