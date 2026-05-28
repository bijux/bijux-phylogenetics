from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import detect_rogue_taxa, write_rogue_taxon_table


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_write_rogue_taxon_table_writes_ranked_consensus_deltas(
    tmp_path: Path,
) -> None:
    output = tmp_path / "rogue-taxon-ranking.tsv"
    report = detect_rogue_taxa(fixture("rogue_taxon_tree_set.nwk"))

    write_rogue_taxon_table(output, report)

    assert output.read_text(encoding="utf-8") == (
        "rank\ttaxon\tmean_terminal_branch_length\tbaseline_consensus_resolution\tpruned_consensus_resolution\tconsensus_resolution_delta\tbaseline_mean_support_percent\tpruned_mean_support_percent\tmean_support_percent_delta\tbaseline_mean_normalized_robinson_foulds\tpruned_mean_normalized_robinson_foulds\tnormalized_robinson_foulds_stability_delta\tbaseline_rooted_topology_count\tpruned_rooted_topology_count\trooted_topology_count_delta\tbaseline_dominant_topology_frequency\tpruned_dominant_topology_frequency\tdominant_topology_frequency_delta\tpruned_consensus_newick\n"
        "1\tE\t0.1\t0.666666666666667\t1\t0.333333333333333\t80\t100\t20\t0.533333333333333\t0\t0.533333333333333\t5\t1\t4\t0.2\t1\t0.8\t((A:1.02,B:0.1)100:0.12,(C:0.1,D:0.12)100:0.12);\n"
        "2\tA\t1\t0.666666666666667\t0.5\t-0.166666666666667\t80\t80\t0\t0.533333333333333\t0.6\t-0.066666666666667\t5\t4\t1\t0.2\t0.4\t0.2\t(B:0.18,(C:0.1,D:0.1)80:0.1,E:0.12);\n"
        "3\tB\t0.1\t0.666666666666667\t0.5\t-0.166666666666667\t80\t80\t0\t0.533333333333333\t0.6\t-0.066666666666667\t5\t4\t1\t0.2\t0.4\t0.2\t(A:1.08,(C:0.1,D:0.1)80:0.1,E:0.1);\n"
        "4\tC\t0.1\t0.666666666666667\t0.5\t-0.166666666666667\t80\t80\t0\t0.533333333333333\t0.6\t-0.066666666666667\t5\t4\t1\t0.2\t0.4\t0.2\t((A:1,B:0.1)80:0.1,D:0.18,E:0.1);\n"
        "5\tD\t0.1\t0.666666666666667\t0.5\t-0.166666666666667\t80\t80\t0\t0.533333333333333\t0.6\t-0.066666666666667\t5\t4\t1\t0.2\t0.4\t0.2\t((A:1,B:0.1)80:0.1,C:0.18,E:0.12);\n"
    )
