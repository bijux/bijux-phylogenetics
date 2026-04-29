from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.bayesian import (
    compute_mrbayes_effective_sample_sizes,
    parse_mrbayes_parameter_traces,
    prepare_mrbayes_analysis,
    run_mrbayes_posterior_inference,
    summarize_mrbayes_posterior_trees,
)


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
prefix = nexus_path.with_suffix("")
prefix.with_suffix(".run1.p").write_text(
    "Gen\\tLnL\\tTL\\talpha\\n"
    "0\\t-110.0\\t0.40\\t0.90\\n"
    "100\\t-108.0\\t0.41\\t0.95\\n"
    "200\\t-107.0\\t0.42\\t1.00\\n"
    "300\\t-106.5\\t0.43\\t1.05\\n",
    encoding="utf-8",
)
prefix.with_suffix(".run1.t").write_text(
    "#NEXUS\\n"
    "begin trees;\\n"
    "tree gen1 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "tree gen2 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "tree gen3 = [&R] ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);\\n"
    "tree gen4 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "end;\\n",
    encoding="utf-8",
)
print("warning: mrbayes fixture posterior run", file=sys.stderr)
""",
    )


def test_prepare_mrbayes_analysis_writes_nexus_with_run_settings(tmp_path: Path) -> None:
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
    assert "dimensions ntax=4 nchar=8;" in text
    assert "lset nst=6 rates=gamma;" in text
    assert "mcmcp ngen=2000 nchains=4 samplefreq=50 printfreq=50;" in text


def test_run_mrbayes_and_summarize_posterior_outputs(tmp_path: Path) -> None:
    executable = _fake_mrbayes(tmp_path / "mb-fixture")
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)

    run_report = run_mrbayes_posterior_inference(nexus_path, executable=executable, resume=False)
    resumed = run_mrbayes_posterior_inference(nexus_path, executable=executable, resume=True)
    posterior_tree_path = run_report.output_paths["posterior_trees"]

    consensus_tree, summary = summarize_mrbayes_posterior_trees(posterior_tree_path, burnin_fraction=0.25)

    assert run_report.run.warning_lines == ["warning: mrbayes fixture posterior run"]
    assert resumed.resumed is True
    assert summary.total_tree_count == 4
    assert summary.burnin_tree_count == 1
    assert summary.kept_tree_count == 3
    assert summary.rooted_topology_count == 2
    assert summary.filtered_tree_set_path.exists()
    assert consensus_tree.tip_count == 4


def test_parse_mrbayes_traces_and_compute_effective_sample_sizes(tmp_path: Path) -> None:
    executable = _fake_mrbayes(tmp_path / "mb-fixture")
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)
    run_report = run_mrbayes_posterior_inference(nexus_path, executable=executable, resume=False)
    trace_path = run_report.output_paths["parameter_traces"]

    trace_report = parse_mrbayes_parameter_traces(trace_path)
    ess_report = compute_mrbayes_effective_sample_sizes(trace_path)

    assert trace_report.row_count == 4
    assert trace_report.columns == ["LnL", "TL", "alpha"]
    assert trace_report.rows[1].generation == 100
    assert ess_report.sample_count == 4
    assert [row.parameter for row in ess_report.effective_sample_sizes] == ["LnL", "TL", "alpha"]
    assert all(row.effective_sample_size > 0 for row in ess_report.effective_sample_sizes)
