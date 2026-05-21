from __future__ import annotations

from pathlib import Path
import tempfile

from .._fixtures import (
    comparative_stress_payload,
    large_alignment_stress_payload,
    resolve_stress_tier_config,
    supermatrix_stress_payload,
    table_generation_stress_payload,
    tree_set_stress_payload,
)
from .._measurement import measure_stress_workload
from ..contracts import LargeDatasetStressObservation, LargeDatasetStressSuiteReport


def benchmark_large_dataset_stress_suite(
    *,
    tier: str = "small",
) -> LargeDatasetStressSuiteReport:
    """Benchmark large owned workloads across one governed stress tier."""
    config = resolve_stress_tier_config(tier)
    observations: list[LargeDatasetStressObservation] = []
    limitations = [
        "resource peaks are measured with python tracemalloc where possible and reuse stage-level engine memory observations when an owned workflow already records them",
        "timeout_seconds is a workload budget recorded for review; only engine-backed workflows enforce it internally during execution",
    ]
    with tempfile.TemporaryDirectory(prefix=f"bijux-stress-{config.tier}-") as tmpdir:
        root = Path(tmpdir)
        workloads = [
            lambda: large_alignment_stress_payload(
                root=root / "alignment", config=config
            ),
            lambda: supermatrix_stress_payload(
                root=root / "supermatrix", config=config
            ),
            lambda: tree_set_stress_payload(root=root / "tree-set", config=config),
            lambda: comparative_stress_payload(
                root=root / "comparative", config=config
            ),
            lambda: table_generation_stress_payload(
                root=root / "tables", config=config
            ),
        ]
        for workload in workloads:
            payload, runtime_seconds, peak_memory_bytes, memory_observation_kind = (
                measure_stress_workload(workload)
            )
            observations.append(
                LargeDatasetStressObservation(
                    workload=payload.workload,
                    tier=config.tier,
                    timeout_seconds=config.timeout_seconds,
                    input_size_bytes=payload.input_size_bytes,
                    sequence_count=payload.sequence_count,
                    alignment_length=payload.alignment_length,
                    tree_count=payload.tree_count,
                    taxon_count=payload.taxon_count,
                    locus_count=payload.locus_count,
                    runtime_seconds=round(runtime_seconds, 15),
                    peak_memory_bytes=peak_memory_bytes,
                    memory_observation_kind=memory_observation_kind,
                    output_row_count=payload.output_row_count,
                    notes=payload.notes,
                )
            )
    return LargeDatasetStressSuiteReport(
        tier=config.tier,
        observations=observations,
        limitations=limitations,
    )
