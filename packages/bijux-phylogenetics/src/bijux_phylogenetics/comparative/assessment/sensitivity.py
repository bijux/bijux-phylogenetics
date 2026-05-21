from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative.common import load_comparative_dataset
from bijux_phylogenetics.comparative.continuous import (
    fit_brownian_motion_model,
    fit_ornstein_uhlenbeck_model,
)


@dataclass(slots=True)
class LeaveOneTaxonOutRow:
    """Influence of dropping one taxon from a continuous-trait fit."""

    dropped_taxon: str
    taxon_count: int
    primary_parameter: float
    delta_primary_parameter: float
    log_likelihood: float
    delta_log_likelihood: float


@dataclass(slots=True)
class ComparativeSensitivityReport:
    """Leave-one-taxon-out sensitivity over one comparative trait model."""

    tree_path: Path
    traits_path: Path
    trait: str
    model: str
    baseline_primary_parameter: float
    baseline_log_likelihood: float
    rows: list[LeaveOneTaxonOutRow]
    most_influential_taxa: list[str]


def run_comparative_sensitivity_analysis(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    model: str,
    taxon_column: str | None = None,
) -> ComparativeSensitivityReport:
    """Quantify leave-one-taxon-out sensitivity for a standalone BM or OU fit."""
    if model not in {"brownian", "ou"}:
        raise ValueError(f"unsupported comparative sensitivity model: {model}")
    baseline_bm = fit_brownian_motion_model(
        tree_path, traits_path, trait=trait, taxon_column=taxon_column
    )
    baseline_ou = (
        None
        if model == "brownian"
        else fit_ornstein_uhlenbeck_model(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
        )
    )
    baseline_primary = baseline_bm.rate if model == "brownian" else baseline_ou.alpha
    baseline_log_likelihood = (
        baseline_bm.log_likelihood
        if model == "brownian"
        else baseline_ou.log_likelihood
    )
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=4,
        require_rooted=True,
        require_binary=False,
    )
    rows: list[LeaveOneTaxonOutRow] = []
    root = Path(tree_path)
    table = Path(traits_path)
    for dropped_taxon in dataset.taxa:
        reduced_tree, _ = _write_reduced_comparative_inputs(
            root, dataset.taxa, dropped_taxon
        )
        reduced_traits = _write_reduced_trait_table(
            table, dataset.taxa, dropped_taxon, taxon_column=dataset.taxon_column
        )
        try:
            if model == "brownian":
                reduced = fit_brownian_motion_model(
                    reduced_tree,
                    reduced_traits,
                    trait=trait,
                    taxon_column=dataset.taxon_column,
                )
                primary = reduced.rate
            else:
                reduced = fit_ornstein_uhlenbeck_model(
                    reduced_tree,
                    reduced_traits,
                    trait=trait,
                    taxon_column=dataset.taxon_column,
                )
                primary = reduced.alpha
            rows.append(
                LeaveOneTaxonOutRow(
                    dropped_taxon=dropped_taxon,
                    taxon_count=reduced.taxon_count,
                    primary_parameter=primary,
                    delta_primary_parameter=primary - baseline_primary,
                    log_likelihood=reduced.log_likelihood,
                    delta_log_likelihood=reduced.log_likelihood
                    - baseline_log_likelihood,
                )
            )
        finally:
            reduced_tree.unlink(missing_ok=True)
            reduced_traits.unlink(missing_ok=True)
    most_influential_taxa = [
        row.dropped_taxon
        for row in sorted(
            rows,
            key=lambda row: (
                abs(row.delta_log_likelihood),
                abs(row.delta_primary_parameter),
                row.dropped_taxon,
            ),
            reverse=True,
        )[:3]
    ]
    return ComparativeSensitivityReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        model=model,
        baseline_primary_parameter=baseline_primary,
        baseline_log_likelihood=baseline_log_likelihood,
        rows=sorted(rows, key=lambda row: row.dropped_taxon),
        most_influential_taxa=most_influential_taxa,
    )


def _write_reduced_comparative_inputs(
    tree_path: Path, taxa: list[str], dropped_taxon: str
) -> tuple[Path, list[str]]:
    from bijux_phylogenetics.io.newick import write_newick
    from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa

    kept_taxa = [taxon for taxon in taxa if taxon != dropped_taxon]
    reduced_tree, _ = prune_tree_to_requested_taxa(tree_path, kept_taxa)
    out_path = Path(
        tempfile.mkstemp(prefix=f"bijux-comparative-{dropped_taxon}-", suffix=".nwk")[1]
    )
    write_newick(out_path, reduced_tree)
    return out_path, kept_taxa


def _write_reduced_trait_table(
    traits_path: Path,
    taxa: list[str],
    dropped_taxon: str,
    *,
    taxon_column: str,
) -> Path:
    from bijux_phylogenetics.datasets.study_inputs import (
        load_taxon_table,
        write_taxon_rows,
    )

    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    kept_taxa = {taxon for taxon in taxa if taxon != dropped_taxon}
    rows = [row for row in table.rows if row[table.taxon_column] in kept_taxa]
    out_path = Path(
        tempfile.mkstemp(
            prefix=f"bijux-comparative-{dropped_taxon}-", suffix=traits_path.suffix
        )[1]
    )
    write_taxon_rows(out_path, columns=table.columns, rows=rows)
    return out_path
