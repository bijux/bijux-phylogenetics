from __future__ import annotations

import tempfile
from pathlib import Path

from .large_tree_model_fitting.artifact_outputs import (
    write_large_tree_model_fitting_observation_table,
    write_large_tree_model_fitting_summary_table,
)
from .large_tree_model_fitting.case_definitions import (
    ContinuousCaseDefinition as _ContinuousCaseDefinition,
    DiscreteCaseDefinition as _DiscreteCaseDefinition,
    case_definitions_for_tier as _select_case_definitions_for_tier,
)
from .large_tree_model_fitting.contracts import (
    LargeTreeModelFittingBenchmarkBundle,
    LargeTreeModelFittingBenchmarkReport,
    LargeTreeModelFittingObservation,
    LargeTreeModelFittingThreshold,
)
from .large_tree_model_fitting.observation_runner import (
    evaluate_threshold as _evaluate_threshold_impl,
    run_case as _run_case,
)


def benchmark_large_tree_model_fitting(
    *,
    tier: str = "small",
) -> LargeTreeModelFittingBenchmarkReport:
    """Benchmark large-tree continuous and discrete fits beyond toy governed cases."""
    case_definitions = _case_definitions_for_tier(tier)
    observations: list[LargeTreeModelFittingObservation] = []
    with tempfile.TemporaryDirectory(prefix=f"bijux-large-fit-{tier}-") as tmpdir:
        root = Path(tmpdir)
        for case_definition in case_definitions:
            case_root = root / case_definition.case_id
            case_root.mkdir(parents=True, exist_ok=True)
            observations.append(_run_case(case_definition, case_root))
    return LargeTreeModelFittingBenchmarkReport(
        tier=tier,
        observations=observations,
        case_count=len(observations),
        geiger_match_case_count=sum(
            1 for row in observations if row.matches_geiger_reference is True
        ),
        threshold_pass_case_count=sum(
            1 for row in observations if row.performance_threshold_passed is True
        ),
        too_slow_case_count=sum(1 for row in observations if row.too_slow_review),
        unstable_case_count=sum(1 for row in observations if row.unstable_review),
        limitations=[
            "runtime and peak memory record the owned bijux fit in the active Python process; stored geiger references provide fit-surface comparison rather than same-process memory parity",
            "continuous optimizer_iteration_count records governed profile-evaluation steps because the owned fitcontinuous surface uses bounded one-parameter profile search rather than a separate iterative optimizer counter",
            "the heavy tier adds one 512-taxon Brownian fit and records threshold review explicitly when that case exceeds the governed runtime or memory budget instead of implying broad large-tree optimizer parity",
        ],
    )


def write_large_tree_model_fitting_bundle(
    output_root: Path,
    *,
    tier: str = "small",
) -> LargeTreeModelFittingBenchmarkBundle:
    """Write a governed large-tree benchmark bundle for one tier."""
    if output_root.exists():
        for path in sorted(output_root.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                path.rmdir()
    output_root.mkdir(parents=True, exist_ok=True)
    report = benchmark_large_tree_model_fitting(tier=tier)
    summary_path = write_large_tree_model_fitting_summary_table(
        output_root / "summary.tsv",
        report,
    )
    observation_table_path = write_large_tree_model_fitting_observation_table(
        output_root / "observations.tsv",
        report,
    )
    return LargeTreeModelFittingBenchmarkBundle(
        output_root=output_root,
        summary_path=summary_path,
        observation_table_path=observation_table_path,
    )


def _case_definitions_for_tier(
    tier: str,
) -> tuple[_ContinuousCaseDefinition | _DiscreteCaseDefinition, ...]:
    return _select_case_definitions_for_tier(tier)


def _evaluate_threshold(
    *,
    threshold: LargeTreeModelFittingThreshold,
    runtime_seconds: float | None,
    peak_memory_bytes: int | None,
    optimizer_step_count: int | None,
) -> tuple[bool | None, bool | None, bool | None, bool | None]:
    return _evaluate_threshold_impl(
        threshold=threshold,
        runtime_seconds=runtime_seconds,
        peak_memory_bytes=peak_memory_bytes,
        optimizer_step_count=optimizer_step_count,
    )


__all__ = [
    "LargeTreeModelFittingBenchmarkBundle",
    "LargeTreeModelFittingBenchmarkReport",
    "LargeTreeModelFittingObservation",
    "LargeTreeModelFittingThreshold",
    "_case_definitions_for_tier",
    "_evaluate_threshold",
    "benchmark_large_tree_model_fitting",
    "write_large_tree_model_fitting_bundle",
    "write_large_tree_model_fitting_observation_table",
    "write_large_tree_model_fitting_summary_table",
]
