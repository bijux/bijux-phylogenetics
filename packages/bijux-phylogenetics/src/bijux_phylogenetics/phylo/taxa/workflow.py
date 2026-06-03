from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import (
    detect_missing_trait_values,
    load_taxon_table,
)
from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.trees import detect_tree_format, load_tree


@dataclass(frozen=True, slots=True)
class TaxonWorkflowLossEvent:
    """One stage where a taxon was lost or made unsafe."""

    stage: str
    reason: str


@dataclass(slots=True)
class TaxonWorkflowLossRow:
    """Full workflow loss trace for one taxon."""

    taxon: str
    present_in_tree: bool
    loss_events: list[TaxonWorkflowLossEvent]
    first_loss_stage: str | None
    retained_for_reporting: bool


@dataclass(slots=True)
class TaxonWorkflowLossReport:
    """Trace where taxa disappear across one dataset workflow."""

    tree_path: Path
    metadata_path: Path
    traits_path: Path
    alignment_path: Path | None
    filtered_alignment_path: Path | None
    inference_tree_path: Path | None
    reported_taxa_path: Path | None
    rows: list[TaxonWorkflowLossRow]
    loss_stage_counts: dict[str, int]


@dataclass(slots=True)
class TaxonRunSource:
    """One named run or artifact contributing a taxon set."""

    label: str
    path: Path
    kind: str
    taxa: list[str]


@dataclass(slots=True)
class TaxonStabilityRow:
    """Retention summary for one taxon across repeated runs."""

    taxon: str
    present_in_runs: list[str]
    absent_from_runs: list[str]
    retention_fraction: float
    stable: bool


@dataclass(slots=True)
class TaxonStabilityReport:
    """Taxon retention stability across multiple workflows or artifacts."""

    sources: list[TaxonRunSource]
    shared_taxa: list[str]
    union_taxa: list[str]
    stable_taxa: list[str]
    unstable_taxa: list[str]
    rows: list[TaxonStabilityRow]


def _try_tree(path: Path) -> list[str] | None:
    try:
        tree_format = detect_tree_format(path)
    except Exception:
        return None
    if tree_format not in {"newick", "nexus", "phyloxml"}:
        return None
    try:
        return sorted(load_tree(path, source_format=tree_format).tip_names)
    except Exception:
        return None


def _try_alignment(path: Path) -> list[str] | None:
    try:
        return sorted(record.identifier for record in load_fasta_alignment(path))
    except Exception:
        return None


def _load_taxa_from_any_path(path: Path) -> tuple[str, list[str]]:
    tree_taxa = _try_tree(path)
    if tree_taxa is not None:
        return "tree", tree_taxa
    alignment_taxa = _try_alignment(path)
    if alignment_taxa is not None:
        return "alignment", alignment_taxa
    table = load_taxon_table(path)
    return "table", sorted(table.taxa)


def _missing_trait_taxa(traits_path: Path) -> set[str]:
    report = detect_missing_trait_values(traits_path)
    return {row.taxon for row in report.missing_values}


def build_taxon_workflow_loss_report(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    alignment_path: Path | None = None,
    filtered_alignment_path: Path | None = None,
    inference_tree_path: Path | None = None,
    reported_taxa_path: Path | None = None,
) -> TaxonWorkflowLossReport:
    """Trace where each taxon is lost across a pragmatic phylogenetics workflow."""
    tree_taxa = set(load_tree(tree_path).tip_names)
    metadata_taxa = set(load_taxon_table(metadata_path).taxa)
    trait_taxa = set(load_taxon_table(traits_path).taxa)
    trait_missingness_taxa = _missing_trait_taxa(traits_path)
    alignment_taxa = (
        set()
        if alignment_path is None
        else {record.identifier for record in load_fasta_alignment(alignment_path)}
    )
    filtered_alignment_taxa = (
        set()
        if filtered_alignment_path is None
        else {
            record.identifier
            for record in load_fasta_alignment(filtered_alignment_path)
        }
    )
    inference_tree_taxa = (
        set()
        if inference_tree_path is None
        else set(load_tree(inference_tree_path).tip_names)
    )
    reported_taxa = (
        set()
        if reported_taxa_path is None
        else set(_load_taxa_from_any_path(reported_taxa_path)[1])
    )
    union_taxa = sorted(
        tree_taxa
        | metadata_taxa
        | trait_taxa
        | alignment_taxa
        | filtered_alignment_taxa
        | inference_tree_taxa
        | reported_taxa
    )
    rows: list[TaxonWorkflowLossRow] = []
    stage_counts: dict[str, int] = {}
    for taxon in union_taxa:
        events: list[TaxonWorkflowLossEvent] = []
        if taxon not in tree_taxa:
            events.append(
                TaxonWorkflowLossEvent(stage="tree", reason="absent_from_tree")
            )
        if (
            alignment_path is not None
            and taxon in tree_taxa
            and taxon not in alignment_taxa
        ):
            events.append(
                TaxonWorkflowLossEvent(
                    stage="alignment", reason="absent_from_alignment"
                )
            )
        if filtered_alignment_path is not None:
            if (
                alignment_path is not None
                and taxon in alignment_taxa
                and taxon not in filtered_alignment_taxa
            ):
                events.append(
                    TaxonWorkflowLossEvent(
                        stage="alignment_filtering",
                        reason="removed_during_alignment_filtering",
                    )
                )
            elif (
                alignment_path is None
                and taxon in tree_taxa
                and taxon not in filtered_alignment_taxa
            ):
                events.append(
                    TaxonWorkflowLossEvent(
                        stage="alignment_filtering",
                        reason="absent_from_filtered_alignment",
                    )
                )
        if taxon not in metadata_taxa:
            events.append(
                TaxonWorkflowLossEvent(stage="metadata", reason="absent_from_metadata")
            )
        if taxon not in trait_taxa:
            events.append(
                TaxonWorkflowLossEvent(stage="traits", reason="absent_from_trait_table")
            )
        elif taxon in trait_missingness_taxa:
            events.append(
                TaxonWorkflowLossEvent(
                    stage="trait_missingness", reason="one_or_more_trait_values_missing"
                )
            )
        if (
            inference_tree_path is not None
            and taxon in tree_taxa
            and taxon not in inference_tree_taxa
        ):
            events.append(
                TaxonWorkflowLossEvent(
                    stage="inference", reason="absent_from_inference_tree"
                )
            )
        if (
            reported_taxa_path is not None
            and taxon in tree_taxa
            and taxon not in reported_taxa
        ):
            events.append(
                TaxonWorkflowLossEvent(
                    stage="reporting", reason="absent_from_report_output"
                )
            )
        first_loss_stage = None if not events else events[0].stage
        if first_loss_stage is not None:
            stage_counts[first_loss_stage] = stage_counts.get(first_loss_stage, 0) + 1
        retained_for_reporting = bool(
            taxon in tree_taxa
            and (reported_taxa_path is None or taxon in reported_taxa)
        )
        rows.append(
            TaxonWorkflowLossRow(
                taxon=taxon,
                present_in_tree=taxon in tree_taxa,
                loss_events=events,
                first_loss_stage=first_loss_stage,
                retained_for_reporting=retained_for_reporting
                and not any(event.stage == "trait_missingness" for event in events),
            )
        )
    return TaxonWorkflowLossReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        alignment_path=alignment_path,
        filtered_alignment_path=filtered_alignment_path,
        inference_tree_path=inference_tree_path,
        reported_taxa_path=reported_taxa_path,
        rows=rows,
        loss_stage_counts=dict(sorted(stage_counts.items())),
    )


def load_taxon_run_source(*, label: str, path: Path) -> TaxonRunSource:
    """Load one run source from a tree, alignment, or taxon-keyed table."""
    kind, taxa = _load_taxa_from_any_path(path)
    return TaxonRunSource(label=label, path=path, kind=kind, taxa=taxa)


def build_taxon_stability_report(sources: list[TaxonRunSource]) -> TaxonStabilityReport:
    """Compare taxon retention across repeated workflows or outputs."""
    if len(sources) < 2:
        raise ValueError(
            "taxon stability reporting requires at least two named sources"
        )
    union_taxa = sorted({taxon for source in sources for taxon in source.taxa})
    shared_taxa = sorted(set.intersection(*(set(source.taxa) for source in sources)))
    rows: list[TaxonStabilityRow] = []
    for taxon in union_taxa:
        present_in_runs = [source.label for source in sources if taxon in source.taxa]
        absent_from_runs = [
            source.label for source in sources if taxon not in source.taxa
        ]
        retention_fraction = len(present_in_runs) / len(sources)
        stable = len(present_in_runs) in {0, len(sources)}
        rows.append(
            TaxonStabilityRow(
                taxon=taxon,
                present_in_runs=present_in_runs,
                absent_from_runs=absent_from_runs,
                retention_fraction=retention_fraction,
                stable=stable,
            )
        )
    return TaxonStabilityReport(
        sources=sources,
        shared_taxa=shared_taxa,
        union_taxa=union_taxa,
        stable_taxa=[row.taxon for row in rows if row.stable],
        unstable_taxa=[row.taxon for row in rows if not row.stable],
        rows=rows,
    )
