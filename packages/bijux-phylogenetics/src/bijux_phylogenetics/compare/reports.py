from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_clade_overlap,
    compare_support_values,
    compare_tree_paths,
    detect_clade_changes,
)
from bijux_phylogenetics.render.html import write_html_report


@dataclass(slots=True)
class ComparisonReportBuildResult:
    output_path: Path
    topology: object
    clades: object
    changes: object
    support: object
    branch_lengths: object
    input_checksums: dict[str, str]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_tree_comparison_report(
    left_path: Path, right_path: Path, *, out_path: Path
) -> ComparisonReportBuildResult:
    """Build a deterministic HTML report for two-tree comparison."""
    topology = compare_tree_paths(left_path, right_path)
    clades = compare_clade_overlap([left_path, right_path])
    changes = detect_clade_changes(left_path, right_path)
    support = compare_support_values(left_path, right_path)
    branch_lengths = compare_branch_lengths(left_path, right_path)
    input_checksums = {
        str(left_path): _sha256(left_path),
        str(right_path): _sha256(right_path),
    }
    sections = [
        (
            "comparison-metrics",
            json.dumps(asdict(topology), default=str, indent=2, sort_keys=True),
        ),
        (
            "clade-comparison",
            json.dumps(asdict(clades), default=str, indent=2, sort_keys=True),
        ),
        (
            "clade-changes",
            json.dumps(asdict(changes), default=str, indent=2, sort_keys=True),
        ),
        (
            "support-comparison",
            json.dumps(asdict(support), default=str, indent=2, sort_keys=True),
        ),
        (
            "support-conflicts",
            json.dumps(
                [asdict(row) for row in support.conflicting_clades],
                default=str,
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "branch-length-comparison",
            json.dumps(asdict(branch_lengths), default=str, indent=2, sort_keys=True),
        ),
        ("input-checksums", json.dumps(input_checksums, indent=2, sort_keys=True)),
    ]
    write_html_report(
        title="Bijux Tree Comparison Report", sections=sections, out_path=out_path
    )
    return ComparisonReportBuildResult(
        output_path=out_path,
        topology=topology,
        clades=clades,
        changes=changes,
        support=support,
        branch_lengths=branch_lengths,
        input_checksums=input_checksums,
    )
