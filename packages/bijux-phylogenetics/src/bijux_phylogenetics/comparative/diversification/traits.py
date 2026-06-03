from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import DiversificationAnalysisError

from .models import (
    TraitDependentDiversificationReport,
    TraitDependentDiversificationState,
)
from .trees import (
    find_smallest_covering_node,
    node_age,
    node_depths,
    validate_time_tree_for_diversification,
)


def run_trait_dependent_diversification_analysis(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> TraitDependentDiversificationReport:
    """Summarize simple state-linked diversification rates when trait states form interpretable clades."""
    validate_time_tree_for_diversification(tree_path)
    tree = load_tree(tree_path)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    if trait not in table.columns:
        raise DiversificationAnalysisError(
            f"trait table does not contain column '{trait}'"
        )
    tree_taxa = set(tree.tip_names)
    rows_by_taxon = {
        row[table.taxon_column]: row
        for row in table.rows
        if row[table.taxon_column] in tree_taxa and row[trait].strip()
    }
    observed_states = sorted({row[trait].strip() for row in rows_by_taxon.values()})
    depths = node_depths(tree)
    states: list[TraitDependentDiversificationState] = []
    warnings: list[str] = []
    for state in observed_states:
        taxa = sorted(
            taxon for taxon, row in rows_by_taxon.items() if row[trait].strip() == state
        )
        state_warnings: list[str] = []
        if len(taxa) < 2:
            state_warnings.append("state is represented by fewer than two taxa")
            states.append(
                TraitDependentDiversificationState(
                    state=state,
                    taxon_count=len(taxa),
                    taxa=taxa,
                    monophyletic=False,
                    crown_age=None,
                    diversification_rate=None,
                    warnings=state_warnings,
                )
            )
            continue
        covering_node, covering_taxa = find_smallest_covering_node(tree, set(taxa))
        monophyletic = covering_taxa == taxa
        crown_age = node_age(tree, depths, covering_node)
        diversification_rate = (
            float(format(math.log(len(taxa)) / crown_age, ".15g"))
            if monophyletic and crown_age > 0.0
            else None
        )
        if not monophyletic:
            state_warnings.append("state taxa are not monophyletic in the input tree")
        states.append(
            TraitDependentDiversificationState(
                state=state,
                taxon_count=len(taxa),
                taxa=taxa,
                monophyletic=monophyletic,
                crown_age=crown_age if monophyletic else None,
                diversification_rate=diversification_rate,
                warnings=state_warnings,
            )
        )
        warnings.extend(state_warnings)
    return TraitDependentDiversificationReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        observed_states=observed_states,
        states=states,
        warnings=sorted(set(warnings)),
    )


def write_trait_dependent_diversification_table(
    path: Path, report: TraitDependentDiversificationReport
) -> Path:
    """Export state-linked diversification summaries as a deterministic TSV."""
    rows = [
        {
            "state": row.state,
            "taxon_count": str(row.taxon_count),
            "taxa": ",".join(row.taxa),
            "monophyletic": str(row.monophyletic).lower(),
            "crown_age": "" if row.crown_age is None else format(row.crown_age, ".15g"),
            "diversification_rate": ""
            if row.diversification_rate is None
            else format(row.diversification_rate, ".15g"),
            "warnings": "; ".join(row.warnings),
        }
        for row in report.states
    ]
    return write_taxon_rows(
        path,
        columns=[
            "state",
            "taxon_count",
            "taxa",
            "monophyletic",
            "crown_age",
            "diversification_rate",
            "warnings",
        ],
        rows=rows,
    )
