from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.phylo.topology.tree import TreeNode
from bijux_phylogenetics.runtime.errors import InvalidDistanceMatrixError

from .shared import _iter_ultrametric_violations, _pair_key

_TIE_TOLERANCE = 1e-12
_ZERO_TOLERANCE = 1e-15

ClusterKey = tuple[str, ...]


@dataclass(frozen=True, slots=True)
class LinkageMergeRow:
    """One deterministic merge step during hierarchical linkage clustering."""

    merge_index: int
    left_cluster: str
    right_cluster: str
    left_cluster_size: int
    right_cluster_size: int
    pair_distance: float
    merge_height: float
    resulting_cluster: str
    resulting_cluster_size: int


@dataclass(frozen=True, slots=True)
class LinkageClusterHeightRow:
    """One internal-cluster height recovered during linkage clustering."""

    merge_index: int
    cluster: str
    height: float


@dataclass(slots=True)
class ActiveCluster:
    """One active cluster tracked during deterministic linkage clustering."""

    key: ClusterKey
    node: TreeNode
    height: float
    size: int


def cluster_label(cluster: ClusterKey) -> str:
    """Return the durable human-readable label for one active cluster."""
    return "|".join(cluster)


def cluster_pair_key(
    left: ClusterKey, right: ClusterKey
) -> tuple[ClusterKey, ClusterKey]:
    """Return one stable unordered key for a cluster-distance lookup."""
    return (left, right) if left < right else (right, left)


def normalized_branch_length(value: float) -> float:
    """Normalize tiny floating-point residue to exact zero branch length."""
    if abs(value) <= _ZERO_TOLERANCE:
        return 0.0
    return value


def require_distance(
    distance_lookup: dict[tuple[str, str], float],
    left_identifier: str,
    right_identifier: str,
) -> float:
    """Require one finite directional distance from a full lookup."""
    try:
        value = distance_lookup[(left_identifier, right_identifier)]
    except KeyError as error:
        raise InvalidDistanceMatrixError(
            f"distance matrix is missing pair {left_identifier}/{right_identifier}"
        ) from error
    if not math.isfinite(value):
        raise InvalidDistanceMatrixError(
            f"distance matrix contains non-finite distance for {left_identifier}/{right_identifier}"
        )
    return value


def validate_linkage_distance_lookup(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> None:
    """Validate one full symmetric distance lookup before linkage clustering."""
    if len(identifiers) < 2:
        raise InvalidDistanceMatrixError(
            "distance matrix must contain at least two taxa"
        )
    if len(set(identifiers)) != len(identifiers):
        raise InvalidDistanceMatrixError(
            "distance matrix contains duplicated taxon labels"
        )
    for left_identifier in identifiers:
        diagonal = require_distance(distance_lookup, left_identifier, left_identifier)
        if not math.isclose(diagonal, 0.0, abs_tol=_TIE_TOLERANCE, rel_tol=0.0):
            raise InvalidDistanceMatrixError(
                "distance matrix has nonzero diagonal entries"
            )
        if diagonal < 0.0:
            raise InvalidDistanceMatrixError(
                "distance matrix contains negative distances"
            )
        for right_identifier in identifiers:
            if left_identifier >= right_identifier:
                continue
            left_to_right = require_distance(
                distance_lookup,
                left_identifier,
                right_identifier,
            )
            right_to_left = require_distance(
                distance_lookup,
                right_identifier,
                left_identifier,
            )
            if left_to_right < 0.0 or right_to_left < 0.0:
                raise InvalidDistanceMatrixError(
                    "distance matrix contains negative distances"
                )
            if not math.isclose(
                left_to_right,
                right_to_left,
                rel_tol=_TIE_TOLERANCE,
                abs_tol=_TIE_TOLERANCE,
            ):
                raise InvalidDistanceMatrixError(
                    "distance matrix contains asymmetric directional entries"
                )


def clustering_assumption_warnings(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
    *,
    method_name: str,
) -> tuple[bool, list[str]]:
    """Return one explicit ultrametric-compatibility warning set for a dendrogram method."""
    pair_distances = {
        _pair_key(left_identifier, right_identifier): distance
        for (left_identifier, right_identifier), distance in distance_lookup.items()
        if left_identifier != right_identifier
    }
    violations = _iter_ultrametric_violations(
        identifiers,
        pair_distances,
        tolerance=_TIE_TOLERANCE,
    )
    if not violations:
        return True, []
    return (
        False,
        [
            f"pairwise distances are not ultrametric, so {method_name}'s strict clock-like assumption is violated"
        ],
    )


def choose_linkage_merge_pair(
    active_distances: dict[tuple[ClusterKey, ClusterKey], float],
    active_keys: list[ClusterKey],
    *,
    method_name: str,
) -> tuple[ClusterKey, ClusterKey]:
    """Choose one deterministic minimum-distance cluster pair."""
    best_pair: tuple[ClusterKey, ClusterKey] | None = None
    best_distance: float | None = None
    for left_index, left_key in enumerate(active_keys):
        for right_key in active_keys[left_index + 1 :]:
            pair_key = cluster_pair_key(left_key, right_key)
            pair_distance = active_distances[pair_key]
            if best_distance is None or pair_distance < best_distance - _TIE_TOLERANCE:
                best_pair = pair_key
                best_distance = pair_distance
                continue
            if (
                best_distance is not None
                and math.isclose(
                    pair_distance,
                    best_distance,
                    rel_tol=_TIE_TOLERANCE,
                    abs_tol=_TIE_TOLERANCE,
                )
                and pair_key > best_pair
            ):
                best_pair = pair_key
    if best_pair is None:
        raise InvalidDistanceMatrixError(
            f"distance matrix did not yield a valid {method_name.upper()} merge pair"
        )
    return best_pair
