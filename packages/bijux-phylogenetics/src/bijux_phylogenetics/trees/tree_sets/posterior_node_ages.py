from __future__ import annotations

import csv
from pathlib import Path
from statistics import median

from bijux_phylogenetics.phylo.branch_lengths.node_ages import (
    RootedCladeNodeAgeObservation,
    rooted_clade_node_age_map,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .contracts import PosteriorNodeAgeSummaryReport, PosteriorNodeAgeSummaryRow
from .inventory import _analyze_tree_set, _require_exact_taxa
from .posterior_statistics import (
    effective_sample_size,
    highest_posterior_density_interval,
)
from .topology import _format_clade

_POSTERIOR_NODE_AGE_HPD_MASS = 0.95


def summarize_posterior_node_ages(
    path: Path,
) -> PosteriorNodeAgeSummaryReport:
    """Summarize posterior node ages by rooted clade identity across one dated tree set."""
    analysis = _analyze_tree_set(path)
    exact_taxa = _require_exact_taxa(analysis)
    ages_by_clade: dict[frozenset[str], list[float]] = {}
    observation_by_clade: dict[frozenset[str], RootedCladeNodeAgeObservation] = {}
    for record, tree in zip(analysis.records, analysis.trees, strict=True):
        try:
            _root_age, clade_ages = rooted_clade_node_age_map(
                tree,
                tree_path=analysis.path,
            )
        except PhylogeneticsError as error:
            raise type(error)(
                error.message,
                code=error.code,
                details={**error.details, "source_tree_index": record.index},
            ) from error
        for clade, observation in clade_ages.items():
            observation_by_clade.setdefault(clade, observation)
            ages_by_clade.setdefault(clade, []).append(observation.age)
    rows: list[PosteriorNodeAgeSummaryRow] = []
    for clade, node_ages in sorted(
        ages_by_clade.items(),
        key=lambda item: (-len(item[1]), _format_clade(item[0])),
    ):
        lower_hpd, upper_hpd = highest_posterior_density_interval(
            node_ages,
            mass=_POSTERIOR_NODE_AGE_HPD_MASS,
        )
        matched_tree_count = len(node_ages)
        observation = observation_by_clade[clade]
        rows.append(
            PosteriorNodeAgeSummaryRow(
                clade=_format_clade(clade),
                node_kind=observation.node_kind,
                matched_tree_count=matched_tree_count,
                posterior_tree_count=len(analysis.trees),
                clade_frequency=round(matched_tree_count / len(analysis.trees), 15),
                mean_node_age=round(sum(node_ages) / matched_tree_count, 15),
                median_node_age=round(float(median(node_ages)), 15),
                hpd_95_lower=round(lower_hpd, 15),
                hpd_95_upper=round(upper_hpd, 15),
                effective_sample_size=effective_sample_size(node_ages),
            )
        )
    return PosteriorNodeAgeSummaryReport(
        path=analysis.path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        hpd_mass=_POSTERIOR_NODE_AGE_HPD_MASS,
        rows=rows,
    )


def write_posterior_node_age_summary_table(
    path: Path,
    report: PosteriorNodeAgeSummaryReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "clade",
                "node_kind",
                "matched_tree_count",
                "posterior_tree_count",
                "clade_frequency",
                "mean_node_age",
                "median_node_age",
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
                    "node_kind": row.node_kind,
                    "matched_tree_count": row.matched_tree_count,
                    "posterior_tree_count": row.posterior_tree_count,
                    "clade_frequency": format(row.clade_frequency, ".15g"),
                    "mean_node_age": format(row.mean_node_age, ".15g"),
                    "median_node_age": format(row.median_node_age, ".15g"),
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
