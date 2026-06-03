from __future__ import annotations

import csv
from pathlib import Path

from .contracts import TreeSetCredibleCladeRow, TreeSetCredibleCladeSetReport
from .inventory import _analyze_tree_set, _require_exact_taxa
from .topology import _format_clade


def compute_credible_clade_set(
    path: Path,
    *,
    credible_threshold: float = 0.95,
) -> TreeSetCredibleCladeSetReport:
    """Select clades by descending posterior frequency until the credible threshold is reached."""
    if not 0.0 < credible_threshold <= 1.0:
        raise ValueError(
            "credible_threshold must be greater than 0 and at most 1, "
            f"got {credible_threshold}"
        )
    analysis = _analyze_tree_set(path)
    exact_taxa = _require_exact_taxa(analysis)
    clade_counts = analysis.clade_counts or {}
    tree_count = len(analysis.trees)
    ranked_clades = sorted(
        clade_counts.items(),
        key=lambda item: (-item[1], _format_clade(item[0])),
    )

    included_rows: list[TreeSetCredibleCladeRow] = []
    excluded_rows: list[TreeSetCredibleCladeRow] = []
    cumulative_frequency = 0.0
    threshold_reached = False
    for index, (clade, count) in enumerate(ranked_clades, start=1):
        frequency = round(count / tree_count, 15)
        cumulative_frequency = round(cumulative_frequency + frequency, 15)
        row = TreeSetCredibleCladeRow(
            inclusion_rank=index,
            clade=_format_clade(clade),
            tree_count=count,
            frequency=frequency,
            cumulative_frequency=cumulative_frequency,
        )
        if threshold_reached:
            excluded_rows.append(row)
            continue
        included_rows.append(row)
        if cumulative_frequency >= credible_threshold:
            threshold_reached = True

    included_cumulative_frequency = (
        0.0 if not included_rows else included_rows[-1].cumulative_frequency
    )
    return TreeSetCredibleCladeSetReport(
        path=analysis.path,
        tree_count=tree_count,
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        credible_threshold=credible_threshold,
        included_clade_count=len(included_rows),
        excluded_clade_count=len(excluded_rows),
        included_cumulative_frequency=included_cumulative_frequency,
        included_clades=included_rows,
        excluded_clades=excluded_rows,
    )


def write_credible_clade_set_included_table(
    path: Path,
    report: TreeSetCredibleCladeSetReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "inclusion_rank",
                "clade",
                "tree_count",
                "frequency",
                "cumulative_frequency",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.included_clades:
            writer.writerow(
                {
                    "inclusion_rank": row.inclusion_rank,
                    "clade": row.clade,
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                    "cumulative_frequency": format(row.cumulative_frequency, ".15g"),
                }
            )
    return path


def write_credible_clade_set_excluded_table(
    path: Path,
    report: TreeSetCredibleCladeSetReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "inclusion_rank",
                "clade",
                "tree_count",
                "frequency",
                "cumulative_frequency",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.excluded_clades:
            writer.writerow(
                {
                    "inclusion_rank": row.inclusion_rank,
                    "clade": row.clade,
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                    "cumulative_frequency": format(row.cumulative_frequency, ".15g"),
                }
            )
    return path


def write_credible_clade_set_artifacts(
    out_dir: Path,
    report: TreeSetCredibleCladeSetReport,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    included_table_path = write_credible_clade_set_included_table(
        out_dir / "credible-clades.tsv",
        report,
    )
    excluded_table_path = write_credible_clade_set_excluded_table(
        out_dir / "excluded-clades.tsv",
        report,
    )
    return {
        "credible_clades_path": included_table_path,
        "excluded_clades_path": excluded_table_path,
    }
