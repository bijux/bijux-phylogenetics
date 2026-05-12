from __future__ import annotations

from pathlib import Path
import pytest
import shutil
from statistics import stdev

from bijux_phylogenetics.bayesian import (
    assess_mrbayes_convergence,
    compute_mrbayes_effective_sample_sizes,
    parse_mrbayes_consensus_tree,
    parse_mrbayes_mcmc_diagnostics,
    parse_mrbayes_parameter_traces,
    parse_mrbayes_posterior_tree_samples,
    prepare_mrbayes_analysis,
    render_bayesian_posterior_report,
    run_mrbayes_posterior_inference,
    summarize_mrbayes_parameter_diagnostics,
    summarize_mrbayes_posterior_trees,
    write_mrbayes_parameter_summary_table,
)
from bijux_phylogenetics.errors import EngineWorkflowError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_mrbayes(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

if "--version" in sys.argv[1:]:
    print("MrBayes v3.2.7a fixture")
    raise SystemExit(0)

nexus_path = Path(sys.argv[1])
trace_path = Path(f"{nexus_path}.run1.p")
tree_path = Path(f"{nexus_path}.run1.t")
mcmc_path = Path(f"{nexus_path}.mcmc")
consensus_path = Path(f"{nexus_path}.con.tre")
trace_path.write_text(
    "Gen\\tLnL\\tTL\\talpha\\n"
    "0\\t-110.0\\t0.40\\t0.90\\n"
    "100\\t-108.0\\t0.41\\t0.95\\n"
    "200\\t-107.0\\t0.42\\t1.00\\n"
    "300\\t-106.5\\t0.43\\t1.05\\n",
    encoding="utf-8",
)
tree_path.write_text(
    "#NEXUS\\n"
    "begin trees;\\n"
    "tree gen1 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "tree gen2 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "tree gen3 = [&R] ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);\\n"
    "tree gen4 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "end;\\n",
    encoding="utf-8",
)
mcmc_path.write_text(
    "[ID: 1]\\n"
    "[   Gen -- Generation]\\n"
    "Gen\\tMove$acc_run1\\tSwap(1<>2)$acc(1)\\tAvgStdDev(s)\\n"
    "100\\t0.5\\t0.75\\t0.20\\n"
    "200\\tNA\\t1.0\\t0.10\\n",
    encoding="utf-8",
)
consensus_path.write_text(
    "#NEXUS\\n"
    "begin trees;\\n"
    "tree con_50_majrule = [&R] ((A[&prob=1.0,prob(percent)=\\\"100\\\"]:0.1,B[&prob=1.0,prob(percent)=\\\"100\\\"]:0.1)"
    "[&prob=0.75,prob(percent)=\\\"75\\\"]:0.2,(C[&prob=1.0,prob(percent)=\\\"100\\\"]:0.1,D[&prob=1.0,prob(percent)=\\\"100\\\"]:0.1)"
    "[&prob=0.5,prob(percent)=\\\"50\\\"]:0.2);\\n"
    "end;\\n",
    encoding="utf-8",
)
print("warning: mrbayes fixture posterior run", file=sys.stderr)
""",
    )


def _real_mrbayes_executable() -> str | None:
    return shutil.which("mb")


def test_prepare_mrbayes_analysis_writes_nexus_with_run_settings(
    tmp_path: Path,
) -> None:
    alignment_path = fixture("alignments/example_alignment.fasta")
    nexus_path = tmp_path / "analysis.nex"

    report = prepare_mrbayes_analysis(
        alignment_path,
        nexus_path,
        model="gtr",
        rates="gamma",
        ngen=2000,
        nchains=4,
        samplefreq=50,
        printfreq=50,
        burnin_fraction=0.25,
    )

    text = nexus_path.read_text(encoding="utf-8")
    assert report.taxon_count == 4
    assert report.character_count == 8
    assert report.partition_count == 1
    assert "dimensions ntax=4 nchar=8;" in text
    assert "lset nst=6 rates=gamma;" in text
    assert "mcmcp ngen=2000 nchains=4 samplefreq=50 printfreq=50;" in text


def test_prepare_mrbayes_analysis_writes_partition_definitions_for_multilocus_matrix(
    tmp_path: Path,
) -> None:
    alignment_path = fixture("alignments/example_multilocus_alignment.fasta")
    partition_path = fixture("alignments/example_multilocus_partitions.txt")
    nexus_path = tmp_path / "partitioned-analysis.nex"

    report = prepare_mrbayes_analysis(
        alignment_path,
        nexus_path,
        partition_path=partition_path,
        model="gtr",
        rates="gamma",
    )

    text = nexus_path.read_text(encoding="utf-8")
    assert report.partition_path == partition_path
    assert report.partition_count == 3
    assert report.partition_names == ["gene_alpha", "gene_beta", "gene_gamma"]
    assert report.partition_data_types == ["DNA"]
    assert report.partition_warnings == []
    assert "charset gene_alpha = 1-4;" in text
    assert "charset gene_beta = 5-9;" in text
    assert "charset gene_gamma = 10-12;" in text
    assert "partition loci = 3: gene_alpha, gene_beta, gene_gamma;" in text
    assert "set partition=loci;" in text
    assert "prset applyto=(all) ratepr=variable;" in text
    assert "lset applyto=(all) nst=6 rates=gamma;" in text


def test_prepare_mrbayes_analysis_rejects_partition_datatype_mismatch(
    tmp_path: Path,
) -> None:
    alignment_path = fixture("alignments/example_multilocus_alignment.fasta")
    partition_path = tmp_path / "bad-partitions.txt"
    partition_path.write_text(
        "PROTEIN,gene_alpha = 1-4\nPROTEIN,gene_beta = 5-9\n",
        encoding="utf-8",
    )

    with pytest.raises(
        EngineWorkflowError,
        match="MrBayes preparation requires partition datatypes that match the alignment alphabet DNA",
    ):
        prepare_mrbayes_analysis(
            alignment_path,
            tmp_path / "bad-analysis.nex",
            partition_path=partition_path,
        )


def test_run_mrbayes_and_summarize_posterior_outputs(tmp_path: Path) -> None:
    executable = _fake_mrbayes(tmp_path / "mb-fixture")
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)

    run_report = run_mrbayes_posterior_inference(
        nexus_path, executable=executable, resume=False
    )
    resumed = run_mrbayes_posterior_inference(
        nexus_path, executable=executable, resume=True
    )
    posterior_tree_path = run_report.output_paths["posterior_trees"]

    consensus_tree, summary = summarize_mrbayes_posterior_trees(
        posterior_tree_path, burnin_fraction=0.25
    )

    assert run_report.run.warning_lines == ["warning: mrbayes fixture posterior run"]
    assert resumed.resumed is True
    assert run_report.output_paths["mcmc_diagnostics"].exists()
    assert run_report.output_paths["consensus_tree"].exists()
    assert summary.total_tree_count == 4
    assert summary.burnin_tree_count == 1
    assert summary.kept_tree_count == 3
    assert summary.rooted_topology_count == 2
    assert summary.filtered_tree_set_path.exists()
    assert consensus_tree.tip_count == 4


def test_prepare_mrbayes_analysis_is_accepted_by_real_mrbayes_on_partitioned_input(
    tmp_path: Path,
) -> None:
    executable = _real_mrbayes_executable()
    if executable is None:
        pytest.skip("real MrBayes executable is not available")
    nexus_path = tmp_path / "partitioned-analysis.nex"
    prepare_mrbayes_analysis(
        fixture("alignments/example_multilocus_alignment.fasta"),
        nexus_path,
        partition_path=fixture("alignments/example_multilocus_partitions.txt"),
        model="gtr",
        rates="gamma",
        ngen=20,
        samplefreq=10,
        printfreq=10,
        burnin_fraction=0.25,
    )

    report = run_mrbayes_posterior_inference(
        nexus_path, executable=executable, resume=False
    )

    assert report.output_paths["posterior_trees"].exists()
    assert report.output_paths["parameter_traces"].exists()
    assert report.output_paths["mcmc_diagnostics"].exists()
    assert report.output_paths["consensus_tree"].exists()
    assert (
        parse_mrbayes_parameter_traces(
            report.output_paths["parameter_traces"]
        ).row_count
        > 0
    )
    assert (
        parse_mrbayes_posterior_tree_samples(
            report.output_paths["posterior_trees"]
        ).tree_count
        > 0
    )
    assert (
        parse_mrbayes_mcmc_diagnostics(report.output_paths["mcmc_diagnostics"]).row_count
        > 0
    )
    assert (
        parse_mrbayes_consensus_tree(report.output_paths["consensus_tree"])[1]
        .annotated_node_count
        > 0
    )


def test_parse_mrbayes_traces_and_compute_effective_sample_sizes(
    tmp_path: Path,
) -> None:
    executable = _fake_mrbayes(tmp_path / "mb-fixture")
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)
    run_report = run_mrbayes_posterior_inference(
        nexus_path, executable=executable, resume=False
    )
    trace_path = run_report.output_paths["parameter_traces"]

    trace_report = parse_mrbayes_parameter_traces(trace_path)
    ess_report = compute_mrbayes_effective_sample_sizes(trace_path)

    assert trace_report.row_count == 4
    assert trace_report.columns == ["LnL", "TL", "alpha"]
    assert trace_report.rows[1].generation == 100
    assert ess_report.sample_count == 4
    assert [row.parameter for row in ess_report.effective_sample_sizes] == [
        "LnL",
        "TL",
        "alpha",
    ]
    assert all(
        row.effective_sample_size > 0 for row in ess_report.effective_sample_sizes
    )


def test_parse_mrbayes_posterior_tree_samples_and_consensus_tree(
    tmp_path: Path,
) -> None:
    executable = _fake_mrbayes(tmp_path / "mb-fixture")
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)
    run_report = run_mrbayes_posterior_inference(
        nexus_path, executable=executable, resume=False
    )

    tree_report = parse_mrbayes_posterior_tree_samples(
        run_report.output_paths["posterior_trees"]
    )
    consensus_tree, consensus_report = parse_mrbayes_consensus_tree(
        run_report.output_paths["consensus_tree"]
    )

    assert tree_report.tree_count == 4
    assert tree_report.rooted_tree_count == 4
    assert tree_report.sampled_generations == [1, 2, 3, 4]
    assert tree_report.tip_names == ["A", "B", "C", "D"]
    assert tree_report.trees[0].newick == "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);"
    assert consensus_tree.tip_count == 4
    assert consensus_report.tree_name == "con_50_majrule"
    assert consensus_report.rooted is True
    assert consensus_report.annotated_node_count == 6
    assert consensus_report.minimum_posterior_probability == 0.5
    assert consensus_report.maximum_posterior_probability == 1.0
    assert consensus_report.minimum_posterior_probability_percent == 50.0
    assert consensus_report.maximum_posterior_probability_percent == 100.0


def test_parse_mrbayes_mcmc_diagnostics(
    tmp_path: Path,
) -> None:
    executable = _fake_mrbayes(tmp_path / "mb-fixture")
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)
    run_report = run_mrbayes_posterior_inference(
        nexus_path, executable=executable, resume=False
    )

    report = parse_mrbayes_mcmc_diagnostics(run_report.output_paths["mcmc_diagnostics"])

    assert report.row_count == 2
    assert report.columns == ["Move$acc_run1", "Swap(1<>2)$acc(1)", "AvgStdDev(s)"]
    assert len(report.comment_lines) == 2
    assert report.rows[0].generation == 100
    assert report.rows[0].values["AvgStdDev(s)"] == 0.2
    assert report.rows[1].values["Move$acc_run1"] is None


def test_parse_real_mrbayes_output_fixture() -> None:
    trace_report = parse_mrbayes_parameter_traces(
        fixture("mrbayes/partitioned-analysis.run1.p")
    )
    tree_report = parse_mrbayes_posterior_tree_samples(
        fixture("mrbayes/partitioned-analysis.run1.t")
    )
    mcmc_report = parse_mrbayes_mcmc_diagnostics(
        fixture("mrbayes/partitioned-analysis.mcmc")
    )
    consensus_tree, consensus_report = parse_mrbayes_consensus_tree(
        fixture("mrbayes/partitioned-analysis.con.tre")
    )

    assert trace_report.row_count == 3
    assert trace_report.columns[:3] == ["LnL", "LnPr", "TL{all}"]
    assert trace_report.rows[-1].generation == 20
    assert tree_report.tree_count == 3
    assert tree_report.rooted_tree_count == 0
    assert tree_report.sampled_generations == [0, 10, 20]
    assert sorted(tree_report.tip_names) == [
        "TaxonA",
        "TaxonB",
        "TaxonC",
        "TaxonD",
        "TaxonE",
    ]
    assert mcmc_report.row_count == 1
    assert mcmc_report.rows[0].generation == 20
    assert mcmc_report.rows[0].values["AvgStdDev(s)"] == pytest.approx(0.235702)
    assert consensus_tree.tip_count == 5
    assert consensus_report.tree_name == "con_50_majrule"
    assert consensus_report.rooted is False
    assert consensus_report.annotated_node_count == 6
    assert consensus_report.minimum_posterior_probability == pytest.approx(0.5)
    assert consensus_report.maximum_posterior_probability_percent == pytest.approx(
        100.0
    )


def test_assess_mrbayes_convergence_flags_low_ess_and_mean_drift(
    tmp_path: Path,
) -> None:
    executable = _fake_mrbayes(tmp_path / "mb-fixture")
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)
    run_report = run_mrbayes_posterior_inference(
        nexus_path, executable=executable, resume=False
    )
    trace_path = run_report.output_paths["parameter_traces"]

    report = assess_mrbayes_convergence(
        trace_path, ess_threshold=5.0, mean_shift_threshold=0.1
    )

    assert report.sample_count == 4
    assert report.converged is False
    assert {warning["code"] for warning in report.warnings} == {"low-ess", "mean-drift"}
    assert [summary["parameter"] for summary in report.parameter_summaries] == [
        "LnL",
        "TL",
        "alpha",
    ]


def test_summarize_mrbayes_parameter_diagnostics_supports_burnin_and_hpd(
    tmp_path: Path,
) -> None:
    trace_path = tmp_path / "diagnostics.run1.p"
    summary_path = tmp_path / "diagnostics-summary.tsv"
    posterior_values = [*range(24), 100.0]
    rows = ["Gen\tLnL\tTL"]
    for index, value in enumerate(posterior_values):
        rows.append(f"{index * 10}\t{value}\t0.5")
    trace_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    report = summarize_mrbayes_parameter_diagnostics(trace_path, burnin_fraction=0.2)
    output_path = write_mrbayes_parameter_summary_table(summary_path, report)
    lnl = next(
        summary for summary in report.parameter_summaries if summary.parameter == "LnL"
    )
    kept_values = posterior_values[5:]

    assert report.burnin_row_count == 5
    assert report.kept_row_count == 20
    assert report.first_kept_generation == 50
    assert report.last_kept_generation == 240
    assert lnl.mean == pytest.approx(18.3)
    assert lnl.median == pytest.approx(14.5)
    assert lnl.standard_deviation == pytest.approx(round(stdev(kept_values), 6))
    assert lnl.hpd_95_lower == pytest.approx(5.0)
    assert lnl.hpd_95_upper == pytest.approx(23.0)
    assert output_path == summary_path
    text = summary_path.read_text(encoding="utf-8")
    assert "median\tstandard_deviation" in text
    assert "first_kept_generation\tlast_kept_generation" in text


def test_render_bayesian_posterior_report_writes_consensus_and_convergence_sections(
    tmp_path: Path,
) -> None:
    executable = _fake_mrbayes(tmp_path / "mb-fixture")
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)
    run_report = run_mrbayes_posterior_inference(
        nexus_path, executable=executable, resume=False
    )
    output_path = tmp_path / "posterior-report.html"

    report = render_bayesian_posterior_report(
        posterior_tree_path=run_report.output_paths["posterior_trees"],
        trace_path=run_report.output_paths["parameter_traces"],
        out_path=output_path,
        burnin_fraction=0.25,
        ess_threshold=5.0,
        mean_shift_threshold=0.1,
    )

    html = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert report.kept_tree_count == 3
    assert report.warning_count >= 1
    assert "posterior-summary" in html
    assert "convergence" in html
    assert "clade-frequencies" in html
