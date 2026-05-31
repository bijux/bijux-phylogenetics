from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.likelihood import (
    rooted_topology_fingerprint_from_newick,
    search_nucleotide_likelihood_nni_from_alignment,
    search_nucleotide_likelihood_spr_from_alignment,
    search_nucleotide_likelihood_tbr_from_alignment,
    write_nucleotide_likelihood_nni_artifacts,
)
from bijux_phylogenetics.phylo.topology.clades import informative_rooted_clades

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_nni_retains_equal_best_topologies_and_strict_consensus() -> None:
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_equal_best_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
        equal_best_likelihood_tolerance=1e-15,
        equal_best_tree_cap=4,
    )

    equal_best = report.equal_best_tree_report
    assert equal_best.likelihood_tolerance == 1e-15
    assert equal_best.retention_cap == 4
    assert equal_best.retained_tree_count == 3
    assert equal_best.omitted_tree_count == 0
    assert equal_best.consensus_method == "strict"
    assert len(equal_best.rows) == 3
    assert len({row.topology_fingerprint for row in equal_best.rows}) == 3
    assert all(
        math.isclose(
            row.log_likelihood,
            equal_best.best_log_likelihood,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for row in equal_best.rows
    )
    consensus_tree = loads_newick(equal_best.consensus_newick)
    consensus_clades = informative_rooted_clades(consensus_tree, {"A", "B", "C", "D"})
    assert consensus_clades == {frozenset({"A", "C"})}


def test_likelihood_nni_equal_best_tree_cap_omits_extra_ties() -> None:
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_equal_best_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
        equal_best_likelihood_tolerance=1e-15,
        equal_best_tree_cap=1,
    )

    equal_best = report.equal_best_tree_report
    assert equal_best.retained_tree_count == 1
    assert equal_best.omitted_tree_count == 2
    assert rooted_topology_fingerprint_from_newick(
        equal_best.consensus_newick
    ) == rooted_topology_fingerprint_from_newick(equal_best.rows[0].tree_newick)


def test_write_nucleotide_likelihood_nni_run_json_records_equal_best_tree_report(
    tmp_path: Path,
) -> None:
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_equal_best_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
        equal_best_likelihood_tolerance=1e-15,
        equal_best_tree_cap=4,
    )

    outputs = write_nucleotide_likelihood_nni_artifacts(
        tmp_path / "likelihood-nni-equal-best-run",
        report,
    )

    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["equal_best_tree_report"]["likelihood_tolerance"] == 1e-15
    assert payload["equal_best_tree_report"]["retained_tree_count"] == 3
    assert payload["equal_best_tree_report"]["omitted_tree_count"] == 0
    assert len(payload["equal_best_tree_report"]["rows"]) == 3


def test_likelihood_spr_and_tbr_reports_surface_one_retained_best_tree_by_default() -> (
    None
):
    spr_report = search_nucleotide_likelihood_spr_from_alignment(
        fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
        evaluation_budget=2,
    )
    tbr_report = search_nucleotide_likelihood_tbr_from_alignment(
        fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )

    assert spr_report.equal_best_tree_report.retained_tree_count == 1
    assert spr_report.equal_best_tree_report.omitted_tree_count == 0
    assert (
        spr_report.equal_best_tree_report.rows[0].tree_newick
        == spr_report.final_tree_newick
    )
    assert rooted_topology_fingerprint_from_newick(
        spr_report.equal_best_tree_report.consensus_newick
    ) == rooted_topology_fingerprint_from_newick(spr_report.final_tree_newick)
    assert tbr_report.equal_best_tree_report.retained_tree_count == 1
    assert tbr_report.equal_best_tree_report.omitted_tree_count == 2
    assert (
        tbr_report.equal_best_tree_report.rows[0].tree_newick
        == tbr_report.final_tree_newick
    )
    assert rooted_topology_fingerprint_from_newick(
        tbr_report.equal_best_tree_report.consensus_newick
    ) == rooted_topology_fingerprint_from_newick(tbr_report.final_tree_newick)
