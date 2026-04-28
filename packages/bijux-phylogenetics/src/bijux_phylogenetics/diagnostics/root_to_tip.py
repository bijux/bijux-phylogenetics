from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import _load_tree


@dataclass(slots=True)
class RootToTipDistance:
    tip: str
    distance: float | None


@dataclass(slots=True)
class RootToTipDistanceReport:
    path: Path
    source_format: str
    distances: list[RootToTipDistance]


def compute_root_to_tip_distances(path: Path, *, source_format: str | None = None) -> RootToTipDistanceReport:
    """Compute one root-to-tip distance per leaf."""
    tree = _load_tree(path, source_format=source_format)
    return RootToTipDistanceReport(
        path=path,
        source_format=tree.source_format,
        distances=[
            RootToTipDistance(tip=tip_name, distance=distance)
            for tip_name, distance in tree.root_to_tip_pairs()
            if tip_name is not None
        ],
    )


def write_root_to_tip_tsv(path: Path, report: RootToTipDistanceReport) -> Path:
    """Write root-to-tip distances as a TSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["tip\tdistance"]
    lines.extend(
        f"{row.tip}\t{'' if row.distance is None else format(row.distance, '.15g')}" for row in report.distances
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
