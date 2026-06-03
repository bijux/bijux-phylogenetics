from __future__ import annotations

from pathlib import Path
import tempfile

from .artifact_outputs import (
    write_large_tree_model_fitting_observation_table,
    write_large_tree_model_fitting_summary_table,
)
from .case_definitions import (
    ContinuousCaseDefinition,
    DiscreteCaseDefinition,
    case_definitions_for_tier,
)
from .contracts import (
    LargeTreeModelFittingBenchmarkBundle,
    LargeTreeModelFittingBenchmarkReport,
    LargeTreeModelFittingObservation,
)
from .observation_runner import run_case


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
            observations.append(run_case(case_definition, case_root))
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
) -> tuple[ContinuousCaseDefinition | DiscreteCaseDefinition, ...]:
    return case_definitions_for_tier(tier)
