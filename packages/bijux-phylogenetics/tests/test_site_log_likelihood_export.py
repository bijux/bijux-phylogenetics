from __future__ import annotations

from collections import Counter
import csv
import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_f81_tree_likelihood,
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
    evaluate_nucleotide_site_log_likelihoods_from_alignment,
    write_site_log_likelihood_table,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_nucleotide_site_log_likelihood_rows_match_selected_model_totals() -> None:
    jc69_tree = load_tree(fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"))
    jc69_records = load_fasta_alignment(
        fixture("alignments", "jc69_site_pattern_alignment.fasta")
    )
    jc69_total = evaluate_jc69_tree_likelihood(jc69_tree, jc69_records)
    jc69_report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="jc69",
    )

    assert jc69_report.model_name == "JC69"
    assert jc69_report.site_count == 10
    assert jc69_report.pattern_count == 6
    assert len(jc69_report.site_log_likelihoods) == 10
    assert math.isclose(
        jc69_report.log_likelihood,
        jc69_total.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        sum(row.log_likelihood for row in jc69_report.site_log_likelihoods),
        jc69_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    k80_tree = load_tree(fixture("trees", "k80_likelihood_tree_2_taxa.nwk"))
    k80_records = load_fasta_alignment(
        fixture("alignments", "k80_likelihood_alignment_2_taxa.fasta")
    )
    k80_total = evaluate_k80_tree_likelihood(k80_tree, k80_records, kappa=4.0)
    k80_report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "k80_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "k80_likelihood_alignment_2_taxa.fasta"),
        model_name="k80",
        kappa=4.0,
    )

    assert k80_report.parameter_values == {"kappa": 4.0}
    assert len(k80_report.site_log_likelihoods) == 4
    assert math.isclose(
        k80_report.log_likelihood,
        k80_total.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    f81_tree = load_tree(fixture("trees", "f81_likelihood_tree_2_taxa.nwk"))
    f81_records = load_fasta_alignment(
        fixture("alignments", "f81_likelihood_alignment_2_taxa.fasta")
    )
    f81_base_frequencies = {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}
    f81_total = evaluate_f81_tree_likelihood(
        f81_tree,
        f81_records,
        base_frequencies=f81_base_frequencies,
    )
    f81_report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "f81_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "f81_likelihood_alignment_2_taxa.fasta"),
        model_name="f81",
        base_frequencies=f81_base_frequencies,
    )

    assert math.isclose(
        f81_report.parameter_values["base_frequency_a"],
        0.4,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        f81_report.log_likelihood,
        f81_total.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    hky85_tree = load_tree(fixture("trees", "hky85_likelihood_tree_2_taxa.nwk"))
    hky85_records = load_fasta_alignment(
        fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta")
    )
    hky85_base_frequencies = {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}
    hky85_total = evaluate_hky85_tree_likelihood(
        hky85_tree,
        hky85_records,
        kappa=4.0,
        base_frequencies=hky85_base_frequencies,
    )
    hky85_report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "hky85_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta"),
        model_name="hky85",
        kappa=4.0,
        base_frequencies=hky85_base_frequencies,
    )

    assert math.isclose(
        hky85_report.parameter_values["kappa"],
        4.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        hky85_report.log_likelihood,
        hky85_total.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    gtr_tree = load_tree(fixture("trees", "gtr_likelihood_tree_2_taxa.nwk"))
    gtr_records = load_fasta_alignment(
        fixture("alignments", "gtr_likelihood_alignment_2_taxa.fasta")
    )
    gtr_base_frequencies = {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}
    gtr_exchangeabilities = {
        "AC": 1.0,
        "AG": 4.5,
        "AT": 0.8,
        "CG": 1.6,
        "CT": 2.4,
        "GT": 3.1,
    }
    gtr_total = evaluate_gtr_tree_likelihood(
        gtr_tree,
        gtr_records,
        exchangeabilities=gtr_exchangeabilities,
        base_frequencies=gtr_base_frequencies,
    )
    gtr_report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "gtr_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "gtr_likelihood_alignment_2_taxa.fasta"),
        model_name="gtr",
        exchangeabilities=gtr_exchangeabilities,
        base_frequencies=gtr_base_frequencies,
    )

    assert math.isclose(
        gtr_report.parameter_values["exchangeability_ac"],
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        gtr_report.parameter_values["exchangeability_ag"],
        4.5,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        gtr_report.log_likelihood,
        gtr_total.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_jc69_site_log_likelihood_export_expands_repeated_patterns() -> None:
    report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="jc69",
    )
    counts_by_pattern = Counter(row.pattern_id for row in report.site_log_likelihoods)
    weights_by_pattern = {
        row.pattern_id: row.pattern_weight for row in report.site_log_likelihoods
    }
    positions_by_pattern = {
        pattern_id: sorted(
            row.site_position
            for row in report.site_log_likelihoods
            if row.pattern_id == pattern_id
        )
        for pattern_id in counts_by_pattern
    }

    assert report.expansion_policy == "expanded-site-rows"
    assert counts_by_pattern == {
        "pattern-1": 2,
        "pattern-2": 2,
        "pattern-3": 1,
        "pattern-4": 1,
        "pattern-5": 1,
        "pattern-6": 3,
    }
    assert weights_by_pattern == {
        "pattern-1": 2,
        "pattern-2": 2,
        "pattern-3": 1,
        "pattern-4": 1,
        "pattern-5": 1,
        "pattern-6": 3,
    }
    assert positions_by_pattern == {
        "pattern-1": [1, 2],
        "pattern-2": [3, 4],
        "pattern-3": [5],
        "pattern-4": [6],
        "pattern-5": [7],
        "pattern-6": [8, 9, 10],
    }


def test_write_site_log_likelihood_table_writes_expanded_tsv(tmp_path: Path) -> None:
    report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="jc69",
    )
    path = write_site_log_likelihood_table(
        tmp_path / "site_log_likelihoods.tsv", report
    )

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert path.name == "site_log_likelihoods.tsv"
    assert len(rows) == 10
    assert rows[0]["model_name"] == "JC69"
    assert rows[0]["taxon_order"] == "A|B|C|D"
    assert rows[0]["pattern_id"] == "pattern-1"
    assert rows[0]["pattern_weight"] == "2"
    assert rows[0]["site_position"] == "1"
    assert rows[0]["site_states"] == "A|A|G|T"
    assert [
        row["site_position"] for row in rows if row["pattern_id"] == "pattern-6"
    ] == [
        "8",
        "9",
        "10",
    ]
    assert math.isclose(
        sum(float(row["log_likelihood"]) for row in rows),
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
