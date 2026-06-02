from __future__ import annotations

import json
from pathlib import Path
from statistics import stdev

import pytest

from bijux_phylogenetics.bayesian import (
    assess_mrbayes_burnin_sensitivity,
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
from bijux_phylogenetics.bayesian.posterior_sets.tree_sets import (
    subsample_mrbayes_posterior_tree_set,
    write_posterior_tree_subsample,
    write_posterior_tree_subsample_table,
)
from bijux_phylogenetics.runtime.errors import (
    EngineUnavailableError,
    EngineWorkflowError,
)

pytestmark = pytest.mark.engine_contract

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

if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
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


def _fake_mrbayes_version_variant(
    path: Path,
    *,
    version_text: str,
    log_likelihood: float,
) -> Path:
    return _write_executable(
        path,
        f"""#!/usr/bin/env python3
import sys
from pathlib import Path

if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
    print("{version_text}")
    raise SystemExit(0)

nexus_path = Path(sys.argv[1])
trace_path = Path(f"{{nexus_path}}.run1.p")
tree_path = Path(f"{{nexus_path}}.run1.t")
mcmc_path = Path(f"{{nexus_path}}.mcmc")
consensus_path = Path(f"{{nexus_path}}.con.tre")
trace_path.write_text(
    "Gen\\tLnL\\tTL\\talpha\\n"
    "0\\t{log_likelihood:.1f}\\t0.40\\t0.90\\n"
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
""",
    )


def _fake_mrbayes_timeout(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
import time

if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
    print("MrBayes v3.2.7a fixture")
    raise SystemExit(0)

time.sleep(1.0)
""",
    )


def _fake_mrbayes_malformed_outputs(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
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
    "100\\tbad\\t0.41\\t0.95\\n",
    encoding="utf-8",
)
tree_path.write_text(
    "#NEXUS\\n"
    "begin trees;\\n"
    "tree gen1 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "end;\\n",
    encoding="utf-8",
)
mcmc_path.write_text(
    "Gen\\tMove$acc_run1\\tSwap(1<>2)$acc(1)\\tAvgStdDev(s)\\n"
    "100\\t0.5\\t0.75\\t0.20\\n",
    encoding="utf-8",
)
consensus_path.write_text(
    "#NEXUS\\n"
    "begin trees;\\n"
    "tree con_50_majrule = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "end;\\n",
    encoding="utf-8",
)
""",
    )


def _fake_mrbayes_missing_consensus_output(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
    print("MrBayes v3.2.7a fixture")
    raise SystemExit(0)

nexus_path = Path(sys.argv[1])
trace_path = Path(f"{nexus_path}.run1.p")
tree_path = Path(f"{nexus_path}.run1.t")
mcmc_path = Path(f"{nexus_path}.mcmc")
trace_path.write_text(
    "Gen\\tLnL\\tTL\\talpha\\n"
    "0\\t-110.0\\t0.40\\t0.90\\n"
    "100\\t-108.0\\t0.41\\t0.95\\n",
    encoding="utf-8",
)
tree_path.write_text(
    "#NEXUS\\n"
    "begin trees;\\n"
    "tree gen1 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "tree gen2 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "end;\\n",
    encoding="utf-8",
)
mcmc_path.write_text(
    "Gen\\tMove$acc_run1\\tSwap(1<>2)$acc(1)\\tAvgStdDev(s)\\n"
    "100\\t0.5\\t0.75\\t0.20\\n",
    encoding="utf-8",
)
""",
    )


def _fake_mrbayes_inconsistent_taxa(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
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
    "[&prob=0.75,prob(percent)=\\\"75\\\"]:0.2,(C[&prob=1.0,prob(percent)=\\\"100\\\"]:0.1,E[&prob=1.0,prob(percent)=\\\"100\\\"]:0.1)"
    "[&prob=0.5,prob(percent)=\\\"50\\\"]:0.2);\\n"
    "end;\\n",
    encoding="utf-8",
)
""",
    )


def _fake_mrbayes_killed(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import os
import signal
import sys

if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
    print("MrBayes v3.2.7a fixture")
    raise SystemExit(0)

os.kill(os.getpid(), signal.SIGTERM)
""",
    )


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
    assert summary.filtered_tree_set_path.parent != posterior_tree_path.parent
    assert consensus_tree.tip_count == 4


def test_run_mrbayes_rebuilds_outputs_after_version_change(tmp_path: Path) -> None:
    first_executable = _fake_mrbayes_version_variant(
        tmp_path / "mb-fixture-first",
        version_text="MrBayes v3.2.7a fixture",
        log_likelihood=-110.0,
    )
    second_executable = _fake_mrbayes_version_variant(
        tmp_path / "mb-fixture-second",
        version_text="MrBayes v3.2.8 fixture",
        log_likelihood=-111.0,
    )
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)

    run_mrbayes_posterior_inference(
        nexus_path,
        executable=first_executable,
        resume=False,
    )
    rebuilt = run_mrbayes_posterior_inference(
        nexus_path,
        executable=second_executable,
        resume=True,
    )

    assert rebuilt.resumed is False
    trace_text = rebuilt.output_paths["parameter_traces"].read_text(encoding="utf-8")
    assert "-111.0" in trace_text


def test_run_mrbayes_posterior_inference_reports_missing_executable_without_marker(
    tmp_path: Path,
) -> None:
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)

    with pytest.raises(EngineUnavailableError, match="was not found"):
        run_mrbayes_posterior_inference(
            nexus_path,
            executable=tmp_path / "missing-mrbayes",
        )

    assert list(tmp_path.glob("*.incomplete.json")) == []


def test_run_mrbayes_posterior_inference_reports_missing_analysis_with_structured_error(
    tmp_path: Path,
) -> None:
    missing_path = tmp_path / "missing-analysis.nex"

    with pytest.raises(EngineWorkflowError) as error:
        run_mrbayes_posterior_inference(missing_path)

    assert error.value.code == "mrbayes_analysis_missing_file"
    assert error.value.details["artifact_kind"] == "mrbayes-analysis-nexus"
    assert error.value.details["path"] == str(missing_path)


def test_run_mrbayes_posterior_inference_times_out_and_marks_incomplete_run(
    tmp_path: Path,
) -> None:
    executable = _fake_mrbayes_timeout(tmp_path / "mb-timeout")
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)

    with pytest.raises(EngineWorkflowError, match="timed out"):
        run_mrbayes_posterior_inference(
            nexus_path,
            executable=executable,
            timeout_seconds=0.5,
        )

    marker_candidates = sorted(tmp_path.glob("*.incomplete.json"))
    assert len(marker_candidates) == 1
    marker_text = marker_candidates[0].read_text(encoding="utf-8")
    assert '"timed_out": true' in marker_text


def test_run_mrbayes_posterior_inference_marks_killed_process_incomplete(
    tmp_path: Path,
) -> None:
    executable = _fake_mrbayes_killed(tmp_path / "mb-killed")
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)

    with pytest.raises(EngineWorkflowError, match="failed with exit code"):
        run_mrbayes_posterior_inference(
            nexus_path,
            executable=executable,
        )

    marker_candidates = sorted(tmp_path.glob("*.incomplete.json"))
    assert len(marker_candidates) == 1
    marker_text = marker_candidates[0].read_text(encoding="utf-8")
    assert '"exit_code": -15' in marker_text


def test_run_mrbayes_posterior_inference_rejects_or_cleans_malformed_outputs(
    tmp_path: Path,
) -> None:
    malformed = _fake_mrbayes_malformed_outputs(tmp_path / "mb-malformed")
    valid = _fake_mrbayes(tmp_path / "mb-valid")
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)

    with pytest.raises(
        EngineWorkflowError,
        match="non-numeric value for 'LnL'",
    ) as error:
        run_mrbayes_posterior_inference(
            nexus_path,
            executable=malformed,
        )

    assert error.value.code == "mrbayes_trace_invalid_parameter_value"
    manifest_path = nexus_path.with_suffix("").with_suffix(".manifest.json")
    marker_path = manifest_path.with_suffix(".incomplete.json")
    assert marker_path.exists()
    marker_payload = json.loads(marker_path.read_text(encoding="utf-8"))
    assert marker_payload["failure_reason"] == "mrbayes_trace_invalid_parameter_value"
    assert marker_payload["missing_output_names"] == []
    observed_outputs = {
        item["output_name"]: item for item in marker_payload["observed_outputs"]
    }
    assert observed_outputs["parameter_traces"]["exists"] is True
    assert observed_outputs["parameter_traces"]["path_kind"] == "file"
    assert observed_outputs["posterior_trees"]["exists"] is True
    assert observed_outputs["posterior_trees"]["path_kind"] == "file"
    assert observed_outputs["consensus_tree"]["exists"] is True
    assert observed_outputs["consensus_tree"]["path_kind"] == "file"

    with pytest.raises(EngineWorkflowError, match="incomplete outputs") as rejected:
        run_mrbayes_posterior_inference(
            nexus_path,
            executable=valid,
            resume=True,
            incomplete_run_policy="reject",
        )
    assert rejected.value.code == "engine_incomplete_outputs_present"
    assert (
        rejected.value.details["failure_reason"]
        == "mrbayes_trace_invalid_parameter_value"
    )
    assert (
        rejected.value.details["observed_outputs"] == marker_payload["observed_outputs"]
    )

    report = run_mrbayes_posterior_inference(
        nexus_path,
        executable=valid,
        resume=True,
        incomplete_run_policy="clean",
    )
    assert report.output_paths["parameter_traces"].exists()
    assert report.run.runtime_seconds >= 0.0
    assert report.config == {"timeout_seconds": None}
    assert marker_path.exists() is False


def test_run_mrbayes_posterior_inference_reports_structured_missing_outputs(
    tmp_path: Path,
) -> None:
    executable = _fake_mrbayes_missing_consensus_output(
        tmp_path / "mb-missing-consensus"
    )
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)

    with pytest.raises(EngineWorkflowError) as error:
        run_mrbayes_posterior_inference(
            nexus_path,
            executable=executable,
        )

    assert error.value.code == "engine_required_output_missing"
    assert error.value.details["workflow"] == "posterior-tree-inference"
    assert error.value.details["engine_name"] == "MrBayes"
    assert error.value.details["missing_outputs"] == [
        {
            "output_name": "consensus_tree",
            "path": str(Path(f"{nexus_path}.con.tre")),
        }
    ]


def test_run_mrbayes_posterior_inference_rejects_inconsistent_taxa(
    tmp_path: Path,
) -> None:
    executable = _fake_mrbayes_inconsistent_taxa(tmp_path / "mb-inconsistent")
    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(fixture("alignments/example_alignment.fasta"), nexus_path)

    with pytest.raises(EngineWorkflowError) as error:
        run_mrbayes_posterior_inference(
            nexus_path,
            executable=executable,
        )

    assert error.value.code == "mrbayes_outputs_inconsistent_taxa"
    assert error.value.details["consensus_tree_taxa"] == ["A", "B", "C", "E"]
    assert error.value.details["posterior_tree_taxa"] == ["A", "B", "C", "D"]


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


def test_parse_mrbayes_parameter_traces_accepts_warning_heavy_fixture() -> None:
    report = parse_mrbayes_parameter_traces(
        fixture("engine_outputs/mrbayes/trace-warning-heavy.run1.p")
    )

    assert report.row_count == 2
    assert report.columns == ["LnL", "TL", "alpha"]
    assert report.rows[0].generation == 0
    assert report.rows[1].values["alpha"] == pytest.approx(0.95)


def test_parse_mrbayes_parameter_traces_reports_truncated_fixture() -> None:
    with pytest.raises(EngineWorkflowError) as error:
        parse_mrbayes_parameter_traces(
            fixture("engine_outputs/mrbayes/trace-truncated.run1.p")
        )

    assert error.value.code == "mrbayes_trace_missing_parameter_value"
    assert error.value.details["artifact_kind"] == "mrbayes-trace"
    assert error.value.details["column"] == "alpha"
    assert error.value.details["row_number"] == 3
    assert error.value.details["expected_section"] == "sampled parameter row"


def test_parse_mrbayes_parameter_traces_reports_malformed_fixture() -> None:
    with pytest.raises(EngineWorkflowError) as error:
        parse_mrbayes_parameter_traces(
            fixture("engine_outputs/mrbayes/trace-malformed.run1.p")
        )

    assert error.value.code == "mrbayes_trace_invalid_parameter_value"
    assert error.value.details["artifact_kind"] == "mrbayes-trace"
    assert error.value.details["column"] == "LnL"
    assert error.value.details["row_number"] == 3
    assert error.value.details["expected_section"] == "sampled parameter row"


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


def test_parse_mrbayes_mcmc_diagnostics_accepts_warning_heavy_fixture() -> None:
    report = parse_mrbayes_mcmc_diagnostics(
        fixture("engine_outputs/mrbayes/mcmc-warning-heavy.mcmc")
    )

    assert report.row_count == 2
    assert len(report.comment_lines) == 2
    assert report.columns == ["Move$acc_run1", "Swap(1<>2)$acc(1)", "AvgStdDev(s)"]
    assert report.rows[1].values["Move$acc_run1"] is None


def test_parse_mrbayes_mcmc_diagnostics_reports_truncated_fixture() -> None:
    with pytest.raises(EngineWorkflowError) as error:
        parse_mrbayes_mcmc_diagnostics(
            fixture("engine_outputs/mrbayes/mcmc-truncated.mcmc")
        )

    assert error.value.code == "mrbayes_mcmc_missing_parameter_value"
    assert error.value.details["artifact_kind"] == "mrbayes-mcmc"
    assert error.value.details["column"] == "AvgStdDev(s)"
    assert error.value.details["row_number"] == 3
    assert error.value.details["expected_section"] == "sampled diagnostics row"


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


def test_parse_mrbayes_tree_and_consensus_accept_version_variant_fixtures() -> None:
    tree_report = parse_mrbayes_posterior_tree_samples(
        fixture("engine_outputs/mrbayes/posterior-tree-version-variant.run1.t")
    )
    consensus_tree, consensus_report = parse_mrbayes_consensus_tree(
        fixture("engine_outputs/mrbayes/consensus-version-variant.con.tre")
    )

    assert tree_report.tree_count == 1
    assert tree_report.rooted_tree_count == 0
    assert tree_report.sampled_generations == [10]
    assert tree_report.tip_names == ["Taxon A", "Taxon B's sample", "Taxon_C"]
    assert "Taxon B''s sample" in tree_report.trees[0].newick
    assert consensus_tree.tip_names == ["Taxon A", "Taxon B's sample", "Taxon_C"]
    assert consensus_report.annotated_node_count == 2
    assert consensus_report.minimum_posterior_probability == pytest.approx(0.75)
    assert consensus_report.maximum_posterior_probability_percent == pytest.approx(
        100.0
    )


def test_parse_mrbayes_posterior_tree_samples_reports_invalid_translate_fixture() -> (
    None
):
    with pytest.raises(EngineWorkflowError) as error:
        parse_mrbayes_posterior_tree_samples(
            fixture("engine_outputs/mrbayes/posterior-tree-invalid-translate.run1.t")
        )

    assert error.value.code == "mrbayes_tree_invalid_translate_block"
    assert error.value.details["artifact_kind"] == "mrbayes-posterior-trees"
    assert error.value.details["expected_section"] == "translate block"


def test_subsample_mrbayes_posterior_tree_set_preserves_generation_metadata(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "posterior.run1.t"
    retained_path = tmp_path / "posterior-subsample.nwk"
    table_path = tmp_path / "posterior-subsample.tsv"
    tree_path.write_text(
        "#NEXUS\n"
        "begin trees;\n"
        "tree gen1 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n"
        "tree gen2 = [&R] ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);\n"
        "tree gen3 = [&R] ((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);\n"
        "tree gen4 = [&R] ((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.4);\n"
        "end;\n",
        encoding="utf-8",
    )

    report = subsample_mrbayes_posterior_tree_set(
        tree_path,
        method="random",
        sample_count=2,
        burnin_fraction=0.25,
        random_seed=11,
    )
    write_posterior_tree_subsample(retained_path, report)
    write_posterior_tree_subsample_table(table_path, report)

    assert report.burnin_tree_count == 1
    assert report.pre_subsampling_tree_count == 3
    assert report.retained_tree_count == 2
    assert report.retained_source_indices == [3, 4]
    assert [tree.generation for tree in report.trees] == [3, 4]
    assert [tree.tree_name for tree in report.trees] == ["gen3", "gen4"]
    assert retained_path.read_text(encoding="utf-8").count("\n") == 2
    assert "gen3\t\t3\ttrue" in table_path.read_text(encoding="utf-8")


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


def test_assess_mrbayes_burnin_sensitivity_reports_parameter_and_clade_instability(
    tmp_path: Path,
) -> None:
    posterior_path = tmp_path / "burnin-sensitive.run1.t"
    trace_path = tmp_path / "burnin-sensitive.run1.p"
    posterior_path.write_text(
        "#NEXUS\n"
        "begin trees;\n"
        "tree gen1 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n"
        "tree gen2 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n"
        "tree gen3 = [&R] ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);\n"
        "tree gen4 = [&R] ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);\n"
        "end;\n",
        encoding="utf-8",
    )
    rows = ["Gen\tLnL\tTL"]
    for index in range(19):
        rows.append(f"{index * 10}\t0.0\t0.5")
    rows.append("190\t100.0\t0.5")
    trace_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    report = assess_mrbayes_burnin_sensitivity(
        posterior_path,
        trace_path=trace_path,
        burnin_fractions=(0.0, 0.95),
    )

    ab_shift = next(shift for shift in report.clade_shifts if shift.clade == "A|B")
    lnl_shift = next(
        shift for shift in report.parameter_shifts if shift.parameter == "LnL"
    )

    assert report.changed_consensus_count == 1
    assert report.unstable_parameter_count >= 1
    assert report.unstable_clade_count >= 1
    assert report.slices[1].first_kept_generation == 190
    assert report.slices[1].lnl_mean == 100.0
    assert ab_shift.crosses_majority_threshold is True
    assert lnl_shift.unstable is True
    assert lnl_shift.common_hpd_95_lower is None


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
    assert report.method_tier.tier == "parser-only"
    assert report.method_tier.inference_mode == "parser-only"
    assert "method-tier" in html
    assert "parser-only" in html
    assert "posterior-summary" in html
    assert "convergence" in html
    assert "clade-frequencies" in html
    assert "limitations" in html
    assert "limitations" in report.machine_manifest["sections"]
