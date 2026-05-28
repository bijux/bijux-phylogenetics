from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..contracts import MultispeciesCoalescentReport


def write_multispecies_coalescent_event_table(
    path: Path,
    report: MultispeciesCoalescentReport,
) -> Path:
    """Write one row per multispecies-coalescent event."""
    return write_taxon_rows(
        path,
        columns=[
            "event_index",
            "species_branch",
            "branch_role",
            "descendant_species",
            "population_size",
            "branch_start_age",
            "branch_end_age",
            "event_age",
            "waiting_time",
            "input_lineage_count",
            "output_lineage_count",
            "left_gene_clade",
            "right_gene_clade",
            "resulting_gene_clade",
        ],
        rows=[
            {
                "event_index": row.event_index,
                "species_branch": row.species_branch,
                "branch_role": row.branch_role,
                "descendant_species": "|".join(row.descendant_species),
                "population_size": format(row.population_size, ".15g"),
                "branch_start_age": format(row.branch_start_age, ".15g"),
                "branch_end_age": (
                    ""
                    if row.branch_end_age is None
                    else format(row.branch_end_age, ".15g")
                ),
                "event_age": format(row.event_age, ".15g"),
                "waiting_time": format(row.waiting_time, ".15g"),
                "input_lineage_count": row.input_lineage_count,
                "output_lineage_count": row.output_lineage_count,
                "left_gene_clade": "|".join(row.left_gene_clade),
                "right_gene_clade": "|".join(row.right_gene_clade),
                "resulting_gene_clade": "|".join(row.resulting_gene_clade),
            }
            for row in report.event_rows
        ],
    )


def write_multispecies_coalescent_branch_table(
    path: Path,
    report: MultispeciesCoalescentReport,
) -> Path:
    """Write one row per species-tree branch traversal in one gene-tree simulation."""
    return write_taxon_rows(
        path,
        columns=[
            "species_branch",
            "branch_role",
            "descendant_species",
            "branch_duration",
            "population_size",
            "lineage_count_entering",
            "coalescent_event_count",
            "lineage_count_exiting",
            "extra_lineage_count",
            "included_in_deep_coalescence_total",
        ],
        rows=[
            {
                "species_branch": row.species_branch,
                "branch_role": row.branch_role,
                "descendant_species": "|".join(row.descendant_species),
                "branch_duration": (
                    ""
                    if row.branch_duration is None
                    else format(row.branch_duration, ".15g")
                ),
                "population_size": format(row.population_size, ".15g"),
                "lineage_count_entering": row.lineage_count_entering,
                "coalescent_event_count": row.coalescent_event_count,
                "lineage_count_exiting": row.lineage_count_exiting,
                "extra_lineage_count": row.extra_lineage_count,
                "included_in_deep_coalescence_total": row.included_in_deep_coalescence_total,
            }
            for row in report.branch_rows
        ],
    )
