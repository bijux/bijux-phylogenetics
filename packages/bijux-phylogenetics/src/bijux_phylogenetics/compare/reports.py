from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from bijux_phylogenetics.compare.topology import compare_branch_lengths, compare_support_values, compare_tree_paths
from bijux_phylogenetics.render.html import write_html_report


@dataclass(slots=True)
class ComparisonReportBuildResult:
    output_path: Path
    topology: object
    support: object
    branch_lengths: object
    input_checksums: dict[str, str]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_tree_comparison_report(left_path: Path, right_path: Path, *, out_path: Path) -> ComparisonReportBuildResult:
    """Build a deterministic HTML report for two-tree comparison."""
    topology = compare_tree_paths(left_path, right_path)
    support = compare_support_values(left_path, right_path)
    branch_lengths = compare_branch_lengths(left_path, right_path)
    input_checksums = {str(left_path): _sha256(left_path), str(right_path): _sha256(right_path)}
    sections = [
        ("comparison-metrics", json.dumps(asdict(topology), default=str, indent=2, sort_keys=True)),
        ("support-comparison", json.dumps(asdict(support), default=str, indent=2, sort_keys=True)),
        ("branch-length-comparison", json.dumps(asdict(branch_lengths), default=str, indent=2, sort_keys=True)),
        ("input-checksums", json.dumps(input_checksums, indent=2, sort_keys=True)),
    ]
    write_html_report(title="Bijux Tree Comparison Report", sections=sections, out_path=out_path)
    return ComparisonReportBuildResult(
        output_path=out_path,
        topology=topology,
        support=support,
        branch_lengths=branch_lengths,
        input_checksums=input_checksums,
    )
