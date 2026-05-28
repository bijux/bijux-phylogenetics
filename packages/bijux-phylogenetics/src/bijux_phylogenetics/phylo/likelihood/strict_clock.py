from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .models import StrictClockLikelihoodReport


def fit_strict_clock_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model: str = "jc69",
    lower_clock_rate_bound: float = 1e-6,
    upper_clock_rate_bound: float = 5.0,
) -> StrictClockLikelihoodReport:
    raise NotImplementedError


def fit_strict_clock_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model: str = "jc69",
    lower_clock_rate_bound: float = 1e-6,
    upper_clock_rate_bound: float = 5.0,
) -> StrictClockLikelihoodReport:
    raise NotImplementedError


def write_strict_clock_branch_table(
    path: Path,
    report: StrictClockLikelihoodReport,
) -> Path:
    raise NotImplementedError


def write_strict_clock_run_json(
    path: Path,
    report: StrictClockLikelihoodReport,
) -> Path:
    raise NotImplementedError


def write_strict_clock_likelihood_artifacts(
    out_dir: Path,
    report: StrictClockLikelihoodReport,
) -> dict[str, Path]:
    raise NotImplementedError
