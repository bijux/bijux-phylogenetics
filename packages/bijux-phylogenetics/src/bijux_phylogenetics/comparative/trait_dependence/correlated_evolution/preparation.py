from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .contracts import CorrelatedTraitExclusion


@dataclass(slots=True)
class _PreparedTraitRows:
    analysis_kind: str
    analyzed_taxa: list[str]
    analyzed_rows: list[dict[str, str]]
    excluded_taxa: list[CorrelatedTraitExclusion]
    left_state_order: list[str]
    right_state_order: list[str]
    warnings: list[str]


def _prepare_shared_trait_rows(
    *,
    tree,
    table,
    left_trait: str,
    right_trait: str,
    analysis_kind: str,
) -> _PreparedTraitRows:
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    tree_taxa = set(tree.tip_names)
    table_taxa = set(table.taxa)
    excluded_taxa: list[CorrelatedTraitExclusion] = []
    excluded_taxa.extend(
        CorrelatedTraitExclusion(
            taxon=taxon,
            reason="missing_from_trait_table",
            missing_traits=[left_trait, right_trait],
        )
        for taxon in sorted(tree_taxa - table_taxa)
    )
    excluded_taxa.extend(
        CorrelatedTraitExclusion(
            taxon=taxon,
            reason="missing_from_tree",
            missing_traits=[],
        )
        for taxon in sorted(table_taxa - tree_taxa)
    )
    candidate_rows: list[dict[str, str]] = []
    for taxon in sorted(tree_taxa & table_taxa):
        row = rows_by_taxon[taxon]
        missing_traits = [
            trait for trait in (left_trait, right_trait) if not row[trait].strip()
        ]
        if missing_traits:
            excluded_taxa.append(
                CorrelatedTraitExclusion(
                    taxon=taxon,
                    reason="missing_trait_value",
                    missing_traits=missing_traits,
                )
            )
            continue
        candidate_rows.append(row)
    resolved_kind = _resolve_analysis_kind(
        candidate_rows=candidate_rows,
        left_trait=left_trait,
        right_trait=right_trait,
        requested=analysis_kind,
    )
    warnings: list[str] = []
    if resolved_kind == "continuous-brownian-contrasts":
        analyzed_rows: list[dict[str, str]] = []
        for row in candidate_rows:
            invalid_traits = [
                trait
                for trait in (left_trait, right_trait)
                if _parse_float_or_none(row[trait]) is None
            ]
            if invalid_traits:
                excluded_taxa.append(
                    CorrelatedTraitExclusion(
                        taxon=row[table.taxon_column],
                        reason="non_numeric_trait_value",
                        missing_traits=invalid_traits,
                    )
                )
                continue
            analyzed_rows.append(row)
        if len(analyzed_rows) < 4:
            raise ComparativeMethodError(
                "continuous correlated-trait evolution requires at least four analyzable taxa"
            )
        return _PreparedTraitRows(
            analysis_kind=resolved_kind,
            analyzed_taxa=[row[table.taxon_column] for row in analyzed_rows],
            analyzed_rows=analyzed_rows,
            excluded_taxa=sorted(excluded_taxa, key=lambda row: row.taxon),
            left_state_order=[],
            right_state_order=[],
            warnings=warnings,
        )
    analyzed_rows = candidate_rows
    if len(analyzed_rows) < 4:
        raise ComparativeMethodError(
            "binary correlated-trait evolution requires at least four analyzable taxa"
        )
    left_state_order = sorted({row[left_trait] for row in analyzed_rows})
    right_state_order = sorted({row[right_trait] for row in analyzed_rows})
    if len(left_state_order) != 2:
        raise ComparativeMethodError(
            f"binary correlated-trait evolution requires exactly two observed states for '{left_trait}'"
        )
    if len(right_state_order) != 2:
        raise ComparativeMethodError(
            f"binary correlated-trait evolution requires exactly two observed states for '{right_trait}'"
        )
    if len({f"{row[left_trait]}|{row[right_trait]}" for row in analyzed_rows}) < 4:
        warnings.append(
            "one or more binary joint states are absent, so binary coupling inference may be weakly identified"
        )
    warnings.append(
        "binary correlated-trait evolution uses a joint-state discrete transition pseudo-likelihood rather than a full Pagel maximum-likelihood fit"
    )
    return _PreparedTraitRows(
        analysis_kind=resolved_kind,
        analyzed_taxa=[row[table.taxon_column] for row in analyzed_rows],
        analyzed_rows=analyzed_rows,
        excluded_taxa=sorted(excluded_taxa, key=lambda row: row.taxon),
        left_state_order=left_state_order,
        right_state_order=right_state_order,
        warnings=warnings,
    )


def _resolve_analysis_kind(
    *,
    candidate_rows: list[dict[str, str]],
    left_trait: str,
    right_trait: str,
    requested: str,
) -> str:
    if not candidate_rows:
        raise ComparativeMethodError(
            "correlated trait evolution does not retain any shared non-missing taxa"
        )
    left_values = [row[left_trait] for row in candidate_rows]
    right_values = [row[right_trait] for row in candidate_rows]
    if requested == "continuous":
        return "continuous-brownian-contrasts"
    if requested == "binary":
        return "binary-joint-state"
    if len(set(left_values)) == 2 and len(set(right_values)) == 2:
        return "binary-joint-state"
    if all(_parse_float_or_none(value) is not None for value in left_values) and all(
        _parse_float_or_none(value) is not None for value in right_values
    ):
        return "continuous-brownian-contrasts"
    raise ComparativeMethodError(
        "auto correlated-trait analysis requires either two numeric traits or two binary traits"
    )


def _parse_float_or_none(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None
