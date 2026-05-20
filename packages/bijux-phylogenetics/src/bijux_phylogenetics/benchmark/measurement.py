from __future__ import annotations

import time
import tracemalloc

from .contracts import (
    BenchmarkObservation,
    LargeAlignmentScalingObservation,
    LargeTreeSetScalingObservation,
    _StressObservationPayload,
)


def max_runtime_seconds(observations) -> float:
    return round(max(row.runtime_seconds for row in observations), 15)


def max_peak_memory_bytes(observations) -> int:
    return max(row.peak_memory_bytes for row in observations)


def measure(
    label: str, item_count: int, *, replicates: int, callback
) -> BenchmarkObservation:
    runtimes: list[float] = []
    peak_memory = 0
    for _ in range(replicates):
        tracemalloc.start()
        started = time.perf_counter()
        callback()
        elapsed = time.perf_counter() - started
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        runtimes.append(elapsed)
        peak_memory = max(peak_memory, peak)
    return BenchmarkObservation(
        label=label,
        item_count=item_count,
        runtime_seconds=round(sum(runtimes) / len(runtimes), 15),
        peak_memory_bytes=peak_memory,
    )


def measure_large_alignment_observation(
    label: str,
    *,
    sequence_count: int,
    alignment_length: int,
    replicates: int,
    callback,
) -> LargeAlignmentScalingObservation:
    runtimes: list[float] = []
    peak_memory = 0
    aligned_site_count = sequence_count * alignment_length
    for _ in range(replicates):
        tracemalloc.start()
        started = time.perf_counter()
        callback()
        elapsed = time.perf_counter() - started
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        runtimes.append(elapsed)
        peak_memory = max(peak_memory, peak)
    return LargeAlignmentScalingObservation(
        label=label,
        sequence_count=sequence_count,
        alignment_length=alignment_length,
        aligned_site_count=aligned_site_count,
        runtime_seconds=round(sum(runtimes) / len(runtimes), 15),
        peak_memory_bytes=peak_memory,
    )


def measure_large_tree_set_observation(
    label: str,
    *,
    tree_count: int,
    tip_count: int,
    replicates: int,
    callback,
) -> LargeTreeSetScalingObservation:
    runtimes: list[float] = []
    peak_memory = 0
    pair_count = tree_count * max(tree_count - 1, 0) // 2
    for _ in range(replicates):
        tracemalloc.start()
        started = time.perf_counter()
        callback()
        elapsed = time.perf_counter() - started
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        runtimes.append(elapsed)
        peak_memory = max(peak_memory, peak)
    return LargeTreeSetScalingObservation(
        label=label,
        tree_count=tree_count,
        tip_count=tip_count,
        pair_count=pair_count,
        runtime_seconds=round(sum(runtimes) / len(runtimes), 15),
        peak_memory_bytes=peak_memory,
    )


def summarize_memory_observation_kinds(kinds: list[str]) -> str:
    distinct = list(dict.fromkeys(kind for kind in kinds if kind))
    if not distinct:
        return "python-tracemalloc"
    if len(distinct) == 1:
        return distinct[0]
    return "mixed"


def measure_stress_workload(
    callback,
) -> tuple[_StressObservationPayload, float, int, str]:
    tracemalloc.start()
    started = time.perf_counter()
    payload, observed_peak_memory_bytes, observed_memory_kind = callback()
    elapsed_seconds = time.perf_counter() - started
    _, tracemalloc_peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_memory_bytes = max(
        tracemalloc_peak_bytes,
        0 if observed_peak_memory_bytes is None else observed_peak_memory_bytes,
    )
    return (
        payload,
        elapsed_seconds,
        peak_memory_bytes,
        observed_memory_kind or "python-tracemalloc",
    )
