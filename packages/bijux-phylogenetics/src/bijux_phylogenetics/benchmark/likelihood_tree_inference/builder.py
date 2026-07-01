from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import time
import tracemalloc

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.likelihood import (
    bootstrap_nucleotide_likelihood_tree_inference_from_alignment,
    infer_nucleotide_likelihood_tree_from_alignment,
)
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_bipartition,
    canonical_clade_id,
    informative_unrooted_splits,
    robinson_foulds_metrics,
    split_sort_key,
)

from .contracts import (
    MaximumLikelihoodWrapperBenchmarkBundle,
    MaximumLikelihoodWrapperBenchmarkObservation,
    MaximumLikelihoodWrapperBenchmarkReport,
    MaximumLikelihoodWrapperBenchmarkSummaryRow,
    MaximumLikelihoodWrapperSupportComparison,
)
from .wrapper_references import MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS
from .artifact_outputs import (
    write_maximum_likelihood_wrapper_comparison_observation_table,
    write_maximum_likelihood_wrapper_comparison_run_json,
    write_maximum_likelihood_wrapper_comparison_summary_table,
)

_MAXIMUM_LIKELIHOOD_WRAPPER_BENCHMARK_STATUSES = (
    "tolerance-match",
    "expected-model-assumption-difference",
    "unsupported-case",
    "native-bug",
)
_LIKELIHOOD_TOLERANCE = 1e-3
_SUPPORT_TOLERANCE = 1e-9


@dataclass(frozen=True, slots=True)
class _MeasuredNativeTreeInference:
    runtime_seconds: float
    peak_memory_bytes: int
    memory_observation_kind: str
    log_likelihood: float
    tree_newick: str
    topology_signature: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _MeasuredNativeBootstrapSupport:
    runtime_seconds: float
    peak_memory_bytes: int
    memory_observation_kind: str
    log_likelihood: float
    tree_newick: str
    topology_signature: tuple[str, ...]
    support_by_split_id: dict[str, float]


def benchmark_maximum_likelihood_wrapper_comparison() -> (
    MaximumLikelihoodWrapperBenchmarkReport
):
    """Benchmark native ML inference against cached IQ-TREE/FastTree/RAxML-style evidence."""
    native_tree = _measure_native_tree_inference()
    native_support = _measure_native_bootstrap_support()
    observations = [
        _build_iqtree_maximum_likelihood_observation(native_tree),
        _build_iqtree_bootstrap_support_observation(native_support),
        _build_fasttree_maximum_likelihood_observation(native_tree),
        _build_raxml_unsupported_observation(),
    ]
    summary_rows = _build_summary_rows(observations)
    tolerance_match_case_count = sum(
        1 for row in observations if row.status == "tolerance-match"
    )
    expected_model_assumption_difference_case_count = sum(
        1
        for row in observations
        if row.status == "expected-model-assumption-difference"
    )
    unsupported_case_count = sum(
        1 for row in observations if row.status == "unsupported-case"
    )
    native_bug_case_count = sum(1 for row in observations if row.status == "native-bug")
    blocking_case_count = sum(1 for row in observations if row.blocking)
    supported_case_count = sum(1 for row in observations if row.supported)
    return MaximumLikelihoodWrapperBenchmarkReport(
        observations=observations,
        summary_rows=summary_rows,
        case_count=len(observations),
        supported_case_count=supported_case_count,
        tolerance_match_case_count=tolerance_match_case_count,
        expected_model_assumption_difference_case_count=(
            expected_model_assumption_difference_case_count
        ),
        unsupported_case_count=unsupported_case_count,
        native_bug_case_count=native_bug_case_count,
        blocking_case_count=blocking_case_count,
        all_supported_cases_clear=blocking_case_count == 0,
        limitations=[
            "native runtime and peak memory are measured in-process with python tracemalloc, while cached wrapper peak memory comes from macos /usr/bin/time maximum resident set size; absolute memory values therefore remain lane-specific rather than same-instrument parity",
            "the iqtree support lane compares one canonical unrooted split against rooted native bootstrap clades so equivalent internal edges match even when rooting conventions differ",
            "the fasttree lane remains an expected model-assumption difference because FastTree uses an approximate JC+CAT likelihood surface and the cached tiny-alignment support output does not expose a separate comparable split-support table",
            "the raxml lane remains explicitly unsupported until one governed small-reference corpus with runtime, memory, likelihood, topology, and support evidence is checked in",
        ],
    )


def write_maximum_likelihood_wrapper_comparison_bundle(
    output_root: Path,
) -> MaximumLikelihoodWrapperBenchmarkBundle:
    """Write one governed artifact bundle for the maximum-likelihood wrapper benchmark."""
    if output_root.exists():
        for path in sorted(output_root.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                path.rmdir()
    output_root.mkdir(parents=True, exist_ok=True)
    report = benchmark_maximum_likelihood_wrapper_comparison()
    summary_path = write_maximum_likelihood_wrapper_comparison_summary_table(
        output_root / "summary.tsv",
        report,
    )
    observation_table_path = (
        write_maximum_likelihood_wrapper_comparison_observation_table(
            output_root / "observations.tsv",
            report,
        )
    )
    run_json_path = write_maximum_likelihood_wrapper_comparison_run_json(
        output_root / "run.json",
        report,
    )
    return MaximumLikelihoodWrapperBenchmarkBundle(
        output_root=output_root,
        summary_path=summary_path,
        observation_table_path=observation_table_path,
        run_json_path=run_json_path,
    )


def _measure_native_tree_inference() -> _MeasuredNativeTreeInference:
    alignment_path = _alignment_path()

    def run():
        return infer_nucleotide_likelihood_tree_from_alignment(
            alignment_path,
            model_name="jc69",
            search_method="nni",
            start_tree_count=3,
            start_tree_seed=17,
            upper_branch_length_bound=10.0,
        )

    report, runtime_seconds, peak_memory_bytes = _measure_callback(run)
    return _MeasuredNativeTreeInference(
        runtime_seconds=runtime_seconds,
        peak_memory_bytes=peak_memory_bytes,
        memory_observation_kind="python-tracemalloc",
        log_likelihood=report.best_final_log_likelihood,
        tree_newick=report.best_final_tree_newick,
        topology_signature=_unrooted_topology_signature(report.best_final_tree_newick),
    )


def _measure_native_bootstrap_support() -> _MeasuredNativeBootstrapSupport:
    alignment_path = _alignment_path()

    def run():
        return bootstrap_nucleotide_likelihood_tree_inference_from_alignment(
            alignment_path,
            model_name="jc69",
            search_method="nni",
            start_tree_count=2,
            start_tree_seed=5,
            replicate_count=1,
            bootstrap_seed=9,
            upper_branch_length_bound=10.0,
        )

    report, runtime_seconds, peak_memory_bytes = _measure_callback(run)
    return _MeasuredNativeBootstrapSupport(
        runtime_seconds=runtime_seconds,
        peak_memory_bytes=peak_memory_bytes,
        memory_observation_kind="python-tracemalloc",
        log_likelihood=report.reference_log_likelihood,
        tree_newick=report.reference_tree_newick,
        topology_signature=_unrooted_topology_signature(report.reference_tree_newick),
        support_by_split_id=_native_support_by_split_id(report),
    )


def _measure_callback(callback):
    tracemalloc.start()
    started = time.perf_counter()
    try:
        report = callback()
    finally:
        elapsed_seconds = time.perf_counter() - started
        _, peak_memory_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    return report, round(elapsed_seconds, 15), peak_memory_bytes


def _build_iqtree_maximum_likelihood_observation(
    native_tree: _MeasuredNativeTreeInference,
) -> MaximumLikelihoodWrapperBenchmarkObservation:
    reference = MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS[
        "iqtree-jc69-maximum-likelihood"
    ]
    wrapper_tree_newick = str(reference["tree_newick"])
    topology_equal, rf_distance, normalized_rf = _compare_topology(
        native_tree.tree_newick,
        wrapper_tree_newick,
    )
    wrapper_log_likelihood = float(reference["log_likelihood"])
    log_likelihood_delta = native_tree.log_likelihood - wrapper_log_likelihood
    status = (
        "tolerance-match"
        if topology_equal and abs(log_likelihood_delta) <= _LIKELIHOOD_TOLERANCE
        else "native-bug"
    )
    return MaximumLikelihoodWrapperBenchmarkObservation(
        case_id="iqtree-jc69-maximum-likelihood",
        wrapper_engine="IQ-TREE-style",
        native_surface="infer_nucleotide_likelihood_tree_from_alignment",
        wrapper_surface="cached IQ-TREE JC maximum-likelihood reference",
        comparison_policy="topology-and-likelihood-tolerance",
        status=status,
        supported=True,
        blocking=status == "native-bug",
        native_runtime_seconds=native_tree.runtime_seconds,
        native_peak_memory_bytes=native_tree.peak_memory_bytes,
        native_memory_observation_kind=native_tree.memory_observation_kind,
        wrapper_runtime_seconds=float(reference["runtime_seconds"]),
        wrapper_peak_memory_bytes=int(reference["peak_memory_bytes"]),
        wrapper_memory_observation_kind=str(reference["memory_observation_kind"]),
        runtime_delta_seconds=(
            native_tree.runtime_seconds - float(reference["runtime_seconds"])
        ),
        peak_memory_delta_bytes=(
            native_tree.peak_memory_bytes - int(reference["peak_memory_bytes"])
        ),
        native_log_likelihood=native_tree.log_likelihood,
        wrapper_log_likelihood=wrapper_log_likelihood,
        log_likelihood_delta=log_likelihood_delta,
        native_tree_newick=native_tree.tree_newick,
        wrapper_tree_newick=wrapper_tree_newick,
        native_topology_signature=native_tree.topology_signature,
        wrapper_topology_signature=_unrooted_topology_signature(wrapper_tree_newick),
        robinson_foulds_distance=rf_distance,
        normalized_robinson_foulds=normalized_rf,
        topology_equal=topology_equal,
        support_comparisons=[],
        notes=list(reference["notes"]),
    )


def _build_iqtree_bootstrap_support_observation(
    native_support: _MeasuredNativeBootstrapSupport,
) -> MaximumLikelihoodWrapperBenchmarkObservation:
    reference = MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS[
        "iqtree-jc69-ultrafast-bootstrap-support"
    ]
    wrapper_tree_newick = str(reference["tree_newick"])
    support_comparisons = _support_comparisons(
        native_support.support_by_split_id,
        reference["support_rows"],
    )
    max_support_delta = max(
        (abs(row.support_percent_delta) for row in support_comparisons),
        default=0.0,
    )
    status = (
        "tolerance-match"
        if max_support_delta <= _SUPPORT_TOLERANCE and support_comparisons
        else "native-bug"
    )
    return MaximumLikelihoodWrapperBenchmarkObservation(
        case_id="iqtree-jc69-ultrafast-bootstrap-support",
        wrapper_engine="IQ-TREE-style",
        native_surface="bootstrap_nucleotide_likelihood_tree_inference_from_alignment",
        wrapper_surface="cached IQ-TREE JC ultrafast-bootstrap reference",
        comparison_policy="canonical-split-support-tolerance",
        status=status,
        supported=True,
        blocking=status == "native-bug",
        native_runtime_seconds=native_support.runtime_seconds,
        native_peak_memory_bytes=native_support.peak_memory_bytes,
        native_memory_observation_kind=native_support.memory_observation_kind,
        wrapper_runtime_seconds=float(reference["runtime_seconds"]),
        wrapper_peak_memory_bytes=int(reference["peak_memory_bytes"]),
        wrapper_memory_observation_kind=str(reference["memory_observation_kind"]),
        runtime_delta_seconds=(
            native_support.runtime_seconds - float(reference["runtime_seconds"])
        ),
        peak_memory_delta_bytes=(
            native_support.peak_memory_bytes - int(reference["peak_memory_bytes"])
        ),
        native_log_likelihood=native_support.log_likelihood,
        wrapper_log_likelihood=float(reference["log_likelihood"]),
        log_likelihood_delta=(
            native_support.log_likelihood - float(reference["log_likelihood"])
        ),
        native_tree_newick=native_support.tree_newick,
        wrapper_tree_newick=wrapper_tree_newick,
        native_topology_signature=native_support.topology_signature,
        wrapper_topology_signature=_unrooted_topology_signature(wrapper_tree_newick),
        robinson_foulds_distance=0,
        normalized_robinson_foulds=0.0,
        topology_equal=True,
        support_comparisons=support_comparisons,
        notes=list(reference["notes"]),
    )


def _build_fasttree_maximum_likelihood_observation(
    native_tree: _MeasuredNativeTreeInference,
) -> MaximumLikelihoodWrapperBenchmarkObservation:
    reference = MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS[
        "fasttree-jc69-approximate-maximum-likelihood"
    ]
    wrapper_tree_newick = str(reference["tree_newick"])
    topology_equal, rf_distance, normalized_rf = _compare_topology(
        native_tree.tree_newick,
        wrapper_tree_newick,
    )
    wrapper_log_likelihood = float(reference["log_likelihood"])
    return MaximumLikelihoodWrapperBenchmarkObservation(
        case_id="fasttree-jc69-approximate-maximum-likelihood",
        wrapper_engine="FastTree-style",
        native_surface="infer_nucleotide_likelihood_tree_from_alignment",
        wrapper_surface="cached FastTree JC+CAT approximate ML reference",
        comparison_policy="expected-assumption-difference",
        status="expected-model-assumption-difference",
        supported=True,
        blocking=False,
        native_runtime_seconds=native_tree.runtime_seconds,
        native_peak_memory_bytes=native_tree.peak_memory_bytes,
        native_memory_observation_kind=native_tree.memory_observation_kind,
        wrapper_runtime_seconds=float(reference["runtime_seconds"]),
        wrapper_peak_memory_bytes=int(reference["peak_memory_bytes"]),
        wrapper_memory_observation_kind=str(reference["memory_observation_kind"]),
        runtime_delta_seconds=(
            native_tree.runtime_seconds - float(reference["runtime_seconds"])
        ),
        peak_memory_delta_bytes=(
            native_tree.peak_memory_bytes - int(reference["peak_memory_bytes"])
        ),
        native_log_likelihood=native_tree.log_likelihood,
        wrapper_log_likelihood=wrapper_log_likelihood,
        log_likelihood_delta=native_tree.log_likelihood - wrapper_log_likelihood,
        native_tree_newick=native_tree.tree_newick,
        wrapper_tree_newick=wrapper_tree_newick,
        native_topology_signature=native_tree.topology_signature,
        wrapper_topology_signature=_unrooted_topology_signature(wrapper_tree_newick),
        robinson_foulds_distance=rf_distance,
        normalized_robinson_foulds=normalized_rf,
        topology_equal=topology_equal,
        support_comparisons=[],
        notes=list(reference["notes"]),
    )


def _build_raxml_unsupported_observation() -> (
    MaximumLikelihoodWrapperBenchmarkObservation
):
    return MaximumLikelihoodWrapperBenchmarkObservation(
        case_id="raxml-small-maximum-likelihood-reference",
        wrapper_engine="RAxML-style",
        native_surface="none",
        wrapper_surface="governed RAxML small ML benchmark reference",
        comparison_policy="unsupported",
        status="unsupported-case",
        supported=False,
        blocking=False,
        native_runtime_seconds=None,
        native_peak_memory_bytes=None,
        native_memory_observation_kind=None,
        wrapper_runtime_seconds=None,
        wrapper_peak_memory_bytes=None,
        wrapper_memory_observation_kind=None,
        runtime_delta_seconds=None,
        peak_memory_delta_bytes=None,
        native_log_likelihood=None,
        wrapper_log_likelihood=None,
        log_likelihood_delta=None,
        native_tree_newick=None,
        wrapper_tree_newick=None,
        native_topology_signature=(),
        wrapper_topology_signature=(),
        robinson_foulds_distance=None,
        normalized_robinson_foulds=None,
        topology_equal=None,
        support_comparisons=[],
        notes=[
            "the repository does not yet ship one governed RAxML-style small-reference corpus with runtime, memory, likelihood, topology, and support evidence",
            "this lane remains explicitly unsupported until that cached reference corpus exists",
        ],
    )


def _support_comparisons(
    native_support_by_split_id: dict[str, float],
    wrapper_support_rows: object,
) -> list[MaximumLikelihoodWrapperSupportComparison]:
    wrapper_support_by_split_id = {
        canonical_clade_id(frozenset(row["split_taxa"])): float(row["support_percent"])
        for row in list(wrapper_support_rows)
    }
    split_ids = sorted(
        set(native_support_by_split_id) | set(wrapper_support_by_split_id),
        key=lambda value: split_sort_key(frozenset(value.split("|"))),
    )
    comparisons: list[MaximumLikelihoodWrapperSupportComparison] = []
    for split_id in split_ids:
        native_support_percent = native_support_by_split_id.get(split_id)
        wrapper_support_percent = wrapper_support_by_split_id.get(split_id)
        support_percent_delta = (
            None
            if native_support_percent is None or wrapper_support_percent is None
            else native_support_percent - wrapper_support_percent
        )
        comparisons.append(
            MaximumLikelihoodWrapperSupportComparison(
                split_taxa=split_id.split("|"),
                native_support_percent=native_support_percent,
                wrapper_support_percent=wrapper_support_percent,
                support_percent_delta=support_percent_delta,
            )
        )
    return comparisons


def _native_support_by_split_id(report) -> dict[str, float]:
    universe = set(loads_newick(report.reference_tree_newick).tip_names)
    support_by_split_id: dict[str, float] = {}
    for row in report.clade_support_rows:
        split = canonical_bipartition(set(row.descendant_taxa), universe)
        split_id = canonical_clade_id(split)
        existing = support_by_split_id.get(split_id)
        if existing is None:
            support_by_split_id[split_id] = row.support_percent
            continue
        if not math.isclose(existing, row.support_percent, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError(
                f"native bootstrap reported inconsistent support for split {split_id}"
            )
    return support_by_split_id


def _compare_topology(
    native_tree_newick: str,
    wrapper_tree_newick: str,
) -> tuple[bool, int, float]:
    metrics = robinson_foulds_metrics(
        loads_newick(native_tree_newick),
        loads_newick(wrapper_tree_newick),
        shared_taxa=set(loads_newick(native_tree_newick).tip_names),
        rf_mode="unrooted",
    )
    return metrics.distance == 0, metrics.distance, metrics.normalized_distance


def _unrooted_topology_signature(tree_newick: str) -> tuple[str, ...]:
    splits = informative_unrooted_splits(loads_newick(tree_newick))
    return tuple(
        canonical_clade_id(signature)
        for signature in sorted(splits, key=split_sort_key)
    )


def _build_summary_rows(
    observations: list[MaximumLikelihoodWrapperBenchmarkObservation],
) -> list[MaximumLikelihoodWrapperBenchmarkSummaryRow]:
    rows: list[MaximumLikelihoodWrapperBenchmarkSummaryRow] = []
    for status in _MAXIMUM_LIKELIHOOD_WRAPPER_BENCHMARK_STATUSES:
        selected = [row for row in observations if row.status == status]
        rows.append(
            MaximumLikelihoodWrapperBenchmarkSummaryRow(
                status=status,
                case_count=len(selected),
                blocking_case_count=sum(1 for row in selected if row.blocking),
                case_ids=[row.case_id for row in selected],
            )
        )
    return rows


def _alignment_path() -> Path:
    return Path(__file__).resolve().parents[6] / str(
        MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS[
            "iqtree-jc69-maximum-likelihood"
        ]["alignment_path"]
    )
