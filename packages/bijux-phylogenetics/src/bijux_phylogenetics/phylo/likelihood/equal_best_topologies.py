from __future__ import annotations

from dataclasses import dataclass, field
import math

from bijux_phylogenetics.io.newick import dumps_newick, loads_newick
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodEqualBestTreeReport,
    NucleotideLikelihoodEqualBestTreeRow,
)
from bijux_phylogenetics.phylo.topology import rooted_topology_fingerprint
from bijux_phylogenetics.trees.tree_sets.consensus import (
    _build_consensus_tree_with_threshold_from_trees,
)


@dataclass(slots=True)
class _EqualBestTopologyAccumulator:
    best_log_likelihood: float | None = None
    rows_by_topology_fingerprint: dict[str, NucleotideLikelihoodEqualBestTreeRow] = (
        field(default_factory=dict)
    )


def validate_nucleotide_likelihood_equal_best_likelihood_tolerance(
    tolerance: float,
) -> float:
    """Validate one equal-best likelihood retention tolerance."""
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError(
            "equal-best likelihood tolerance must be finite and nonnegative"
        )
    return float(tolerance)


def validate_nucleotide_likelihood_equal_best_tree_cap(cap: int) -> int:
    """Validate one equal-best tree retention cap."""
    if cap < 1:
        raise ValueError("equal-best tree retention cap must be at least one")
    return int(cap)


def initialize_nucleotide_likelihood_equal_best_topology_accumulator() -> (
    _EqualBestTopologyAccumulator
):
    """Start one equal-best topology accumulator for a native likelihood search."""
    return _EqualBestTopologyAccumulator()


def record_nucleotide_likelihood_equal_best_topology(
    accumulator: _EqualBestTopologyAccumulator,
    *,
    tree_newick: str,
    log_likelihood: float,
    likelihood_tolerance: float,
) -> None:
    """Record one evaluated tree into the current equal-best topology set."""
    validated_tolerance = (
        validate_nucleotide_likelihood_equal_best_likelihood_tolerance(
            likelihood_tolerance
        )
    )
    topology_fingerprint = rooted_topology_fingerprint(loads_newick(tree_newick))
    candidate_row = NucleotideLikelihoodEqualBestTreeRow(
        retained_rank=0,
        topology_fingerprint=topology_fingerprint,
        tree_newick=tree_newick,
        log_likelihood=float(log_likelihood),
    )
    if (
        accumulator.best_log_likelihood is None
        or log_likelihood > accumulator.best_log_likelihood + validated_tolerance
    ):
        accumulator.best_log_likelihood = float(log_likelihood)
        accumulator.rows_by_topology_fingerprint = {
            topology_fingerprint: candidate_row,
        }
        return
    if accumulator.best_log_likelihood - log_likelihood > validated_tolerance:
        return
    existing_row = accumulator.rows_by_topology_fingerprint.get(topology_fingerprint)
    if existing_row is None or log_likelihood > existing_row.log_likelihood:
        accumulator.rows_by_topology_fingerprint[topology_fingerprint] = candidate_row


def build_nucleotide_likelihood_equal_best_tree_report(
    accumulator: _EqualBestTopologyAccumulator,
    *,
    likelihood_tolerance: float,
    retention_cap: int,
) -> NucleotideLikelihoodEqualBestTreeReport:
    """Build one retained equal-best tree report from one finished accumulator."""
    validated_tolerance = (
        validate_nucleotide_likelihood_equal_best_likelihood_tolerance(
            likelihood_tolerance
        )
    )
    validated_cap = validate_nucleotide_likelihood_equal_best_tree_cap(retention_cap)
    if (
        accumulator.best_log_likelihood is None
        or not accumulator.rows_by_topology_fingerprint
    ):
        raise ValueError(
            "equal-best topology report requires at least one recorded tree"
        )
    sorted_rows = sorted(
        accumulator.rows_by_topology_fingerprint.values(),
        key=lambda row: (
            -row.log_likelihood,
            row.topology_fingerprint,
            row.tree_newick,
        ),
    )
    retained_rows = [
        NucleotideLikelihoodEqualBestTreeRow(
            retained_rank=index,
            topology_fingerprint=row.topology_fingerprint,
            tree_newick=row.tree_newick,
            log_likelihood=row.log_likelihood,
        )
        for index, row in enumerate(sorted_rows[:validated_cap], start=1)
    ]
    consensus_tree, _ = _build_consensus_tree_with_threshold_from_trees(
        [loads_newick(row.tree_newick) for row in retained_rows],
        threshold=1.0,
    )
    return NucleotideLikelihoodEqualBestTreeReport(
        likelihood_tolerance=validated_tolerance,
        retention_cap=validated_cap,
        retained_tree_count=len(retained_rows),
        omitted_tree_count=max(0, len(sorted_rows) - len(retained_rows)),
        best_log_likelihood=float(accumulator.best_log_likelihood),
        consensus_method="strict",
        consensus_newick=dumps_newick(consensus_tree),
        rows=retained_rows,
    )
