from __future__ import annotations

import csv
from pathlib import Path
from statistics import median

from .contracts import (
    PosteriorBranchLengthSummaryReport,
    PosteriorBranchLengthSummaryRow,
)
from .inventory import _analyze_tree_set, _require_exact_taxa
from .posterior_statistics import (
    effective_sample_size,
    highest_posterior_density_interval,
)
from .topology import _format_clade

_POSTERIOR_BRANCH_LENGTH_HPD_MASS = 0.95


def summarize_posterior_branch_lengths(
    path: Path,
) -> PosteriorBranchLengthSummaryReport:
    """Summarize posterior branch lengths by rooted clade identity across one tree set."""
    analysis = _analyze_tree_set(path)
    exact_taxa = _require_exact_taxa(analysis)
    rows = []
    for clade, branch_lengths in sorted(
        analysis.clade_branch_lengths.items(),
        key=lambda item: (-len(item[1]), _format_clade(item[0])),
    ):
        lower_hpd, upper_hpd = highest_posterior_density_interval(
            branch_lengths,
            mass=_POSTERIOR_BRANCH_LENGTH_HPD_MASS,
        )
        matched_tree_count = len(branch_lengths)
        rows.append(
            PosteriorBranchLengthSummaryRow(
                clade=_format_clade(clade),
                matched_tree_count=matched_tree_count,
                posterior_tree_count=len(analysis.trees),
                clade_frequency=round(matched_tree_count / len(analysis.trees), 15),
                mean_branch_length=round(sum(branch_lengths) / matched_tree_count, 15),
                median_branch_length=round(float(median(branch_lengths)), 15),
                hpd_95_lower=round(lower_hpd, 15),
                hpd_95_upper=round(upper_hpd, 15),
                effective_sample_size=effective_sample_size(branch_lengths),
            )
        )
    return PosteriorBranchLengthSummaryReport(
        path=analysis.path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        hpd_mass=_POSTERIOR_BRANCH_LENGTH_HPD_MASS,
        rows=rows,
    )


def write_posterior_branch_length_summary_table(
    path: Path,
    report: PosteriorBranchLengthSummaryReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "clade",
                "matched_tree_count",
                "posterior_tree_count",
                "clade_frequency",
                "mean_branch_length",
                "median_branch_length",
                "hpd_95_lower",
                "hpd_95_upper",
                "effective_sample_size",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.rows:
            writer.writerow(
                {
                    "clade": row.clade,
                    "matched_tree_count": row.matched_tree_count,
                    "posterior_tree_count": row.posterior_tree_count,
                    "clade_frequency": format(row.clade_frequency, ".15g"),
                    "mean_branch_length": format(row.mean_branch_length, ".15g"),
                    "median_branch_length": format(row.median_branch_length, ".15g"),
                    "hpd_95_lower": format(row.hpd_95_lower, ".15g"),
                    "hpd_95_upper": format(row.hpd_95_upper, ".15g"),
                    "effective_sample_size": (
                        ""
                        if row.effective_sample_size is None
                        else format(row.effective_sample_size, ".15g")
                    ),
                }
            )
    return path
