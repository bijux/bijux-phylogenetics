from __future__ import annotations

import importlib
import json
from pathlib import Path
from statistics import stdev
import xml.etree.ElementTree as ET

import pytest

from bijux_phylogenetics.bayesian.beast import (
    assess_beast_burnin_sensitivity,
    assess_beast_chain_mixing,
    assess_beast_convergence,
    detect_impossible_calibration_constraints,
    parse_beast_log,
    parse_beast_posterior_tree_samples,
    prepare_beast_time_tree_analysis,
    run_beast_posterior_inference,
    summarize_beast_analysis_xml,
    summarize_beast_log,
    summarize_beast_posterior_topology_diversity,
    summarize_beast_posterior_trees,
    validate_beast_analysis_xml,
    validate_beast_posterior_log,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
    write_beast_log_summary_table,
    write_beast_posterior_tree_set,
)
from bijux_phylogenetics.bayesian.evidence import build_bayesian_evidence_package
from bijux_phylogenetics.bayesian.posterior_sets.tree_sets import (
    subsample_beast_posterior_tree_set,
    summarize_maximum_clade_credibility_tree,
    write_posterior_tree_subsample,
    write_posterior_tree_subsample_table,
)
from bijux_phylogenetics.bayesian.presentation.html_reports import (
    render_bayesian_diagnostics_report,
    render_calibration_audit_report,
)
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_beast_posterior_fixture,
)
from bijux_phylogenetics.runtime.errors import (
    EngineUnavailableError,
    EngineWorkflowError,
)
from bijux_phylogenetics.trees import compute_consensus_tree

pytestmark = pytest.mark.engine_contract

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")
SHARED_BEAST_POSTERIOR = get_shared_beast_posterior_fixture(
    "strict_yule_real_posterior"
)


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_core_dataset_imports_after_beast_package_split() -> None:
    dataset_module = importlib.import_module("bijux_phylogenetics.core.dataset")
    assert hasattr(dataset_module, "summarize_dataset_readiness")


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_beast(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "-version" in args:
    print("BEAST v2.7.7 fixture")
    raise SystemExit(0)

xml_path = Path(args[-1])
seed = args[args.index("-seed") + 1]
log_path = xml_path.with_name(f"{xml_path.stem}.{seed}.log")
tree_path = xml_path.with_name(f"{xml_path.stem}.{seed}.trees")
log_path.write_text(
    "Sample\\tposterior\\tlikelihood\\tprior\\ttreeHeight\\tclockRate\\tbirthRate\\n"
    "0\\t-120.0\\t-80.0\\t-40.0\\t1.1\\t0.01\\t0.2\\n"
    "20\\t-118.0\\t-79.0\\t-39.0\\t1.0\\t0.011\\t0.21\\n",
    encoding="utf-8",
)
tree_path.write_text(
    "#NEXUS\\n"
    "Begin trees;\\n"
    "tree STATE_0 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "tree STATE_20 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "End;\\n",
    encoding="utf-8",
)
print("warning: beast fixture posterior run", file=sys.stderr)
""",
    )


def _fake_beast_version_variant(
    path: Path,
    *,
    version_text: str,
    posterior_value: float,
) -> Path:
    return _write_executable(
        path,
        f"""#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "-version" in args:
    print("{version_text}")
    raise SystemExit(0)

xml_path = Path(args[-1])
seed = args[args.index("-seed") + 1]
log_path = xml_path.with_name(f"{{xml_path.stem}}.{{seed}}.log")
tree_path = xml_path.with_name(f"{{xml_path.stem}}.{{seed}}.trees")
log_path.write_text(
    "Sample\\tposterior\\tlikelihood\\tprior\\ttreeHeight\\tclockRate\\tbirthRate\\n"
    "0\\t{posterior_value:.1f}\\t-80.0\\t-40.0\\t1.1\\t0.01\\t0.2\\n"
    "20\\t-118.0\\t-79.0\\t-39.0\\t1.0\\t0.011\\t0.21\\n",
    encoding="utf-8",
)
tree_path.write_text(
    "#NEXUS\\n"
    "Begin trees;\\n"
    "tree STATE_0 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "tree STATE_20 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "End;\\n",
    encoding="utf-8",
)
""",
    )


def _fake_beast_timeout(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
import time

if "-version" in sys.argv[1:]:
    print("BEAST v2.7.7 fixture")
    raise SystemExit(0)

time.sleep(1.0)
""",
    )


def _fake_beast_malformed_outputs(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "-version" in args:
    print("BEAST v2.7.7 fixture")
    raise SystemExit(0)

xml_path = Path(args[-1])
seed = args[args.index("-seed") + 1]
log_path = xml_path.with_name(f"{xml_path.stem}.{seed}.log")
tree_path = xml_path.with_name(f"{xml_path.stem}.{seed}.trees")
log_path.write_text(
    "Sample\\tposterior\\tlikelihood\\tprior\\ttreeHeight\\tclockRate\\tbirthRate\\n"
    "0\\t-120.0\\t-80.0\\t-40.0\\t1.1\\t0.01\\t0.2\\n"
    "20\\tbad\\t-79.0\\t-39.0\\t1.0\\t0.011\\t0.21\\n",
    encoding="utf-8",
)
tree_path.write_text(
    "#NEXUS\\n"
    "Begin trees;\\n"
    "tree STATE_0 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "End;\\n",
    encoding="utf-8",
)
""",
    )


def _fake_beast_missing_tree_output(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "-version" in args:
    print("BEAST v2.7.7 fixture")
    raise SystemExit(0)

xml_path = Path(args[-1])
seed = args[args.index("-seed") + 1]
log_path = xml_path.with_name(f"{xml_path.stem}.{seed}.log")
log_path.write_text(
    "Sample\\tposterior\\tlikelihood\\tprior\\ttreeHeight\\tclockRate\\tbirthRate\\n"
    "0\\t-120.0\\t-80.0\\t-40.0\\t1.1\\t0.01\\t0.2\\n"
    "20\\t-118.0\\t-79.0\\t-39.0\\t1.0\\t0.011\\t0.21\\n",
    encoding="utf-8",
)
""",
    )


def _fake_beast_inconsistent_states(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "-version" in args:
    print("BEAST v2.7.7 fixture")
    raise SystemExit(0)

xml_path = Path(args[-1])
seed = args[args.index("-seed") + 1]
log_path = xml_path.with_name(f"{xml_path.stem}.{seed}.log")
tree_path = xml_path.with_name(f"{xml_path.stem}.{seed}.trees")
log_path.write_text(
    "Sample\\tposterior\\tlikelihood\\tprior\\ttreeHeight\\tclockRate\\tbirthRate\\n"
    "0\\t-120.0\\t-80.0\\t-40.0\\t1.1\\t0.01\\t0.2\\n"
    "20\\t-118.0\\t-79.0\\t-39.0\\t1.0\\t0.011\\t0.21\\n",
    encoding="utf-8",
)
tree_path.write_text(
    "#NEXUS\\n"
    "Begin trees;\\n"
    "tree STATE_0 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "tree STATE_40 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "End;\\n",
    encoding="utf-8",
)
""",
    )


def _fake_beast_killed(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import os
import signal
import sys

if "-version" in sys.argv[1:]:
    print("BEAST v2.7.7 fixture")
    raise SystemExit(0)

os.kill(os.getpid(), signal.SIGTERM)
""",
    )


def test_validate_fossil_calibration_table_accepts_named_and_taxon_targets() -> None:
    report = validate_fossil_calibration_table(
        fixture("example_tree_named_clades.nwk"),
        fixture("example_calibrations.tsv"),
    )

    assert report.calibration_count == 2
    assert report.valid_calibration_count == 2
    assert report.invalid_calibration_count == 0
    assert [calibration.target_kind for calibration in report.calibrations] == [
        "named-clade",
        "taxa",
    ]
    assert report.calibrations[0].taxa == ["A", "B"]
    assert report.calibrations[1].taxa == ["C", "D"]


def test_detect_impossible_calibration_constraints_reports_unknown_and_invalid_targets() -> (
    None
):
    report = detect_impossible_calibration_constraints(
        fixture("example_tree_named_clades.nwk"),
        fixture("example_calibrations_invalid.tsv"),
    )

    assert report.impossible_calibration_ids == ["bad-clade", "bad-empty", "bad-order"]
    assert {issue.code for issue in report.issues} >= {
        "unknown-clade-name",
        "non-monophyletic-target",
        "minimum-exceeds-maximum",
        "missing-target",
        "missing-age-bounds",
    }


def test_validate_tip_dating_metadata_checks_tree_and_alignment_membership() -> None:
    report = validate_tip_dating_metadata(
        fixture("example_tree_named_clades.nwk"),
        fixture("example_tip_dates_invalid.tsv"),
        alignment_path=fixture("example_alignment.fasta"),
    )

    assert report.valid_tip_count == 1
    assert report.invalid_tip_count == 2
    assert report.extra_tip_taxa == ["E"]
    assert report.extra_alignment_taxa == ["E"]
    assert {issue.code for issue in report.issues} >= {
        "taxon-missing-from-tree",
        "taxon-missing-from-alignment",
        "invalid-date",
    }


def test_prepare_beast_time_tree_analysis_writes_clock_prior_calibrations_and_tip_dates(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "analysis.xml"

    report = prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        output_path,
        tree_path=fixture("example_tree_named_clades.nwk"),
        calibration_path=fixture("example_calibrations.tsv"),
        tip_dates_path=fixture("example_tip_dates.tsv"),
        clock_model="relaxed-lognormal",
        tree_prior="birth-death",
        chain_length=200000,
        log_every=500,
    )

    text = output_path.read_text(encoding="utf-8")
    xml = ET.parse(output_path)
    root = xml.getroot()
    assert report.taxon_count == 4
    assert report.character_count == 8
    assert report.calibration_count == 2
    assert report.tip_date_count == 4
    assert report.beast_data_type == "nucleotide"
    assert report.substitution_model == "HKY"
    assert report.starting_tree_source == "provided-tree"
    assert report.warning_count == 3
    assert report.log_path.name == "analysis.$(seed).log"
    assert report.tree_log_path.name == "analysis.$(seed).trees"
    assert root.tag == "beast"
    assert root.find("./data[@id='alignment']") is not None
    assert (
        root.find(
            "./input[@id='branchRates'][@spec='beast.base.evolution.branchratemodel.UCRelaxedClockModel']"
        )
        is not None
    )
    assert (
        root.find(
            "./input[@id='treePrior'][@spec='beast.base.evolution.speciation.BirthDeathGernhard08Model']"
        )
        is not None
    )
    assert root.find("./tree[@id='tree']") is not None
    assert root.find("./tree[@id='tree']/trait[@traitname='date-forward']") is not None
    assert root.find(".//distribution[@id='cal-mammals']") is not None
    assert (
        root.find(
            ".//distribution[@id='cal-mammals']/distr[@spec='beast.base.inference.distribution.Uniform']"
        )
        is not None
    )
    assert (
        root.find(
            ".//distribution[@id='cal-birds']/distr[@spec='beast.base.inference.distribution.Exponential']"
        )
        is not None
    )
    assert (
        "template generator does not infer parametric lognormal shape parameters automatically"
        in " ".join(report.warnings)
    )
    assert (
        "translated a lower-bound-only uniform calibration into an offset exponential prior"
        in " ".join(report.warnings)
    )
    assert "standard birth-death tree prior are exploratory" in " ".join(
        report.warnings
    )
    assert "analysis.$(seed).trees" in text


def test_prepare_beast_time_tree_analysis_supports_protein_alignments_without_starting_tree(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "protein-analysis.xml"

    report = prepare_beast_time_tree_analysis(
        fixture("example_alignment_protein.fasta"),
        output_path,
        clock_model="strict",
        tree_prior="yule",
    )

    root = ET.parse(output_path).getroot()
    assert report.beast_data_type == "aminoacid"
    assert report.substitution_model == "JTT"
    assert report.starting_tree_source == "upgma"
    assert report.warning_count == 0
    assert root.find("./data[@id='alignment'][@dataType='aminoacid']") is not None
    assert root.find("./input[@id='siteModel']/substModel[@spec='JTT']") is not None
    assert (
        root.find("./input[@id='tree'][@spec='beast.base.evolution.tree.ClusterTree']")
        is not None
    )
    assert (
        root.find(
            "./input[@id='branchRates'][@spec='beast.base.evolution.branchratemodel.StrictClockModel']"
        )
        is not None
    )


def test_summarize_beast_analysis_xml_reads_real_fixture() -> None:
    report = summarize_beast_analysis_xml(fixture("beast2_strict_yule_posterior.xml"))

    assert report.valid is True
    assert report.beast_version == "2.7"
    assert report.taxon_count == 4
    assert report.character_count == 8
    assert report.beast_data_type == "nucleotide"
    assert report.substitution_model == "HKY"
    assert report.clock_model == "strict"
    assert report.tree_prior == "yule"
    assert report.starting_tree_source == "upgma"
    assert report.chain_length == 2000
    assert report.posterior_log_path == Path("beast2_strict_yule_posterior.$(seed).log")
    assert report.posterior_tree_path == Path(
        "beast2_strict_yule_posterior.$(seed).trees"
    )
    assert report.calibration_count == 0
    assert report.tip_date_count == 0


def test_validate_beast_analysis_xml_reports_missing_required_outputs(
    tmp_path: Path,
) -> None:
    xml_path = tmp_path / "broken-analysis.xml"
    xml_path.write_text(
        "<?xml version='1.0' encoding='utf-8'?>\n"
        "<beast version='2.7'>\n"
        "  <data id='alignment' dataType='nucleotide'>\n"
        "    <sequence taxon='A'>ACTG</sequence>\n"
        "    <sequence taxon='B'>ACTG</sequence>\n"
        "  </data>\n"
        "  <input spec='HKY' id='hky' />\n"
        "  <input spec='SiteModel' id='siteModel'>\n"
        "    <substModel idref='hky' />\n"
        "  </input>\n"
        "  <input spec='beast.base.evolution.tree.ClusterTree' id='tree' clusterType='upgma'>\n"
        "    <taxa idref='alignment' />\n"
        "  </input>\n"
        "  <input spec='beast.base.evolution.branchratemodel.StrictClockModel' id='branchRates' />\n"
        "  <input spec='beast.base.evolution.speciation.YuleModel' id='treePrior' />\n"
        "  <run spec='MCMC' id='mcmc' chainLength='1000'>\n"
        "    <state><stateNode idref='tree' /></state>\n"
        "    <logger logEvery='100'><log idref='posterior' /></logger>\n"
        "  </run>\n"
        "</beast>\n",
        encoding="utf-8",
    )

    report = validate_beast_analysis_xml(xml_path)

    assert report.valid is False
    assert {issue.code for issue in report.issues} >= {
        "missing-posterior-log",
        "missing-posterior-trees",
    }


def test_summarize_beast_analysis_xml_reports_missing_file_with_structured_error(
    tmp_path: Path,
) -> None:
    missing_path = tmp_path / "missing-analysis.xml"

    with pytest.raises(EngineWorkflowError) as caught:
        summarize_beast_analysis_xml(missing_path)

    assert caught.value.code == "beast_xml_missing_file"
    assert caught.value.details["artifact_kind"] == "beast-analysis-xml"
    assert caught.value.details["path"] == str(missing_path)


def test_run_beast_posterior_inference_writes_outputs_and_resumes(
    tmp_path: Path,
) -> None:
    executable = _fake_beast(tmp_path / "beast-fixture")
    xml_path = tmp_path / "strict-yule.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        xml_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )

    first = run_beast_posterior_inference(
        xml_path,
        executable=executable,
        seed=1,
    )
    resumed = run_beast_posterior_inference(
        xml_path,
        executable=executable,
        seed=1,
        resume=True,
    )

    assert first.output_paths["posterior_log"].exists()
    assert first.output_paths["posterior_trees"].exists()
    assert resumed.resumed is True
    assert parse_beast_log(first.output_paths["posterior_log"]).row_count == 2
    assert (
        parse_beast_posterior_tree_samples(
            first.output_paths["posterior_trees"],
            burnin_fraction=0.0,
        ).kept_tree_count
        == 2
    )


def test_run_beast_posterior_inference_rebuilds_outputs_after_version_change(
    tmp_path: Path,
) -> None:
    first_executable = _fake_beast_version_variant(
        tmp_path / "beast-fixture-first",
        version_text="BEAST v2.7.7 fixture",
        posterior_value=-120.0,
    )
    second_executable = _fake_beast_version_variant(
        tmp_path / "beast-fixture-second",
        version_text="BEAST v2.7.8 fixture",
        posterior_value=-121.0,
    )
    xml_path = tmp_path / "strict-yule.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        xml_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )

    run_beast_posterior_inference(
        xml_path,
        executable=first_executable,
        seed=1,
    )
    rebuilt = run_beast_posterior_inference(
        xml_path,
        executable=second_executable,
        seed=1,
        resume=True,
    )

    assert rebuilt.resumed is False
    log_text = rebuilt.output_paths["posterior_log"].read_text(encoding="utf-8")
    assert "-121.0" in log_text


def test_run_beast_posterior_inference_reports_missing_executable_without_marker(
    tmp_path: Path,
) -> None:
    xml_path = tmp_path / "strict-yule.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        xml_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )

    with pytest.raises(EngineUnavailableError, match="was not found"):
        run_beast_posterior_inference(
            xml_path,
            executable=tmp_path / "missing-beast",
            seed=1,
        )

    assert list(tmp_path.glob("*.incomplete.json")) == []


def test_run_beast_posterior_inference_times_out_and_marks_incomplete_run(
    tmp_path: Path,
) -> None:
    executable = _fake_beast_timeout(tmp_path / "beast-timeout")
    xml_path = tmp_path / "strict-yule.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        xml_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )

    with pytest.raises(EngineWorkflowError, match="timed out"):
        run_beast_posterior_inference(
            xml_path,
            executable=executable,
            seed=1,
            timeout_seconds=0.5,
        )

    marker_candidates = sorted(tmp_path.glob("*.incomplete.json"))
    assert len(marker_candidates) == 1
    marker_text = marker_candidates[0].read_text(encoding="utf-8")
    assert '"timed_out": true' in marker_text


def test_run_beast_posterior_inference_marks_killed_process_incomplete(
    tmp_path: Path,
) -> None:
    executable = _fake_beast_killed(tmp_path / "beast-killed")
    xml_path = tmp_path / "strict-yule.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        xml_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )

    with pytest.raises(EngineWorkflowError, match="failed with exit code"):
        run_beast_posterior_inference(
            xml_path,
            executable=executable,
            seed=1,
        )

    marker_candidates = sorted(tmp_path.glob("*.incomplete.json"))
    assert len(marker_candidates) == 1
    marker_text = marker_candidates[0].read_text(encoding="utf-8")
    assert '"exit_code": -15' in marker_text


def test_run_beast_posterior_inference_rejects_or_cleans_malformed_outputs(
    tmp_path: Path,
) -> None:
    malformed = _fake_beast_malformed_outputs(tmp_path / "beast-malformed")
    valid = _fake_beast(tmp_path / "beast-valid")
    xml_path = tmp_path / "strict-yule.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        xml_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )

    with pytest.raises(EngineWorkflowError, match="non-numeric value"):
        run_beast_posterior_inference(
            xml_path,
            executable=malformed,
            seed=1,
        )

    manifest_path = xml_path.with_suffix(".manifest.json")
    marker_path = manifest_path.with_suffix(".incomplete.json")
    assert marker_path.exists()
    marker_payload = json.loads(marker_path.read_text(encoding="utf-8"))
    assert marker_payload["failure_reason"] == "beast_log_invalid_parameter_value"
    assert marker_payload["missing_output_names"] == []
    observed_outputs = {
        item["output_name"]: item for item in marker_payload["observed_outputs"]
    }
    assert observed_outputs["posterior_log"]["exists"] is True
    assert observed_outputs["posterior_log"]["path_kind"] == "file"
    assert observed_outputs["posterior_trees"]["exists"] is True
    assert observed_outputs["posterior_trees"]["path_kind"] == "file"

    with pytest.raises(EngineWorkflowError, match="incomplete outputs") as rejected:
        run_beast_posterior_inference(
            xml_path,
            executable=valid,
            seed=1,
            resume=True,
            incomplete_run_policy="reject",
        )
    assert rejected.value.code == "engine_incomplete_outputs_present"
    assert (
        rejected.value.details["failure_reason"] == "beast_log_invalid_parameter_value"
    )
    assert (
        rejected.value.details["observed_outputs"] == marker_payload["observed_outputs"]
    )

    report = run_beast_posterior_inference(
        xml_path,
        executable=valid,
        seed=1,
        resume=True,
        incomplete_run_policy="clean",
    )
    assert report.output_paths["posterior_log"].exists()
    assert report.run.runtime_seconds >= 0.0
    assert report.config == {
        "threads": 1,
        "seed": 1,
        "overwrite": True,
        "timeout_seconds": None,
    }
    assert marker_path.exists() is False


def test_run_beast_posterior_inference_reports_structured_missing_outputs(
    tmp_path: Path,
) -> None:
    executable = _fake_beast_missing_tree_output(tmp_path / "beast-missing-tree")
    xml_path = tmp_path / "strict-yule.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        xml_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )

    with pytest.raises(EngineWorkflowError) as error:
        run_beast_posterior_inference(
            xml_path,
            executable=executable,
            seed=1,
        )

    assert error.value.code == "engine_required_output_missing"
    assert error.value.details["workflow"] == "posterior-tree-inference"
    assert error.value.details["engine_name"] == "BEAST"
    assert error.value.details["missing_outputs"] == [
        {
            "output_name": "posterior_trees",
            "path": str(xml_path.with_name("strict-yule.1.trees")),
        }
    ]


def test_run_beast_posterior_inference_rejects_inconsistent_sampled_states(
    tmp_path: Path,
) -> None:
    executable = _fake_beast_inconsistent_states(tmp_path / "beast-inconsistent")
    xml_path = tmp_path / "strict-yule.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        xml_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )

    with pytest.raises(EngineWorkflowError) as error:
        run_beast_posterior_inference(
            xml_path,
            executable=executable,
            seed=1,
        )

    assert error.value.code == "beast_outputs_inconsistent_states"
    assert error.value.details["log_states"] == [0, 20]
    assert error.value.details["tree_states"] == [0, 40]
    assert error.value.details["missing_log_states"] == [40]


def test_prepare_beast_time_tree_analysis_requires_tree_for_calibrations(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "analysis.xml"

    try:
        prepare_beast_time_tree_analysis(
            fixture("example_alignment.fasta"),
            output_path,
            calibration_path=fixture("example_calibrations.tsv"),
        )
    except ValueError as error:
        assert "requires tree_path when calibration_path is provided" in str(error)
    else:
        raise AssertionError(
            "expected BEAST preparation to reject calibrations without a starting tree"
        )


def test_parse_beast_log_and_assess_convergence_return_parameter_summaries() -> None:
    log_report = parse_beast_log(fixture("example_beast.log"))
    convergence = assess_beast_convergence(
        fixture("example_beast.log"),
        ess_threshold=5.0,
        mean_shift_threshold=0.1,
    )

    assert log_report.row_count == 4
    assert log_report.columns == ["posterior", "likelihood", "clockRate", "treeHeight"]
    assert log_report.rows[1].state == 1000
    assert convergence.converged is False
    assert {warning["code"] for warning in convergence.warnings} == {
        "low-ess",
        "mean-drift",
    }


def test_parse_beast_log_reports_structured_missing_and_malformed_errors(
    tmp_path: Path,
) -> None:
    missing_path = tmp_path / "missing.log"
    invalid_path = fixture("engine_outputs/beast/log-malformed.log")

    with pytest.raises(EngineWorkflowError) as missing:
        parse_beast_log(missing_path)
    assert missing.value.code == "beast_log_missing_file"
    assert missing.value.details["artifact_kind"] == "beast-log"
    assert missing.value.details["expected_section"] == "posterior log file"

    with pytest.raises(EngineWorkflowError) as malformed:
        parse_beast_log(invalid_path)
    assert malformed.value.code == "beast_log_invalid_parameter_value"
    assert malformed.value.details["artifact_kind"] == "beast-log"
    assert malformed.value.details["column"] == "posterior"
    assert malformed.value.details["row_number"] == 3
    assert malformed.value.details["expected_section"] == "sampled parameter row"


def test_summarize_beast_log_classifies_parameters_and_writes_summary_table(
    tmp_path: Path,
) -> None:
    log_path = tmp_path / "classified-beast.log"
    summary_path = tmp_path / "classified-beast-summary.tsv"
    log_path.write_text(
        "# BEAST fixture log\n"
        "state\tposterior\tlikelihood\tprior\tclockRate\ttreeHeight\tbirthRate\talpha\n"
        "0\t-510.0\t-490.0\t-20.0\t0.0010\t12.0\t0.30\t0.90\n"
        "1000\t-505.0\t-486.0\t-19.0\t0.0011\t12.3\t0.35\t0.95\n"
        "2000\t-500.0\t-482.0\t-18.0\t0.0012\t12.8\t0.40\t1.00\n"
        "3000\t-497.0\t-479.0\t-18.0\t0.0013\t13.4\t0.45\t1.05\n",
        encoding="utf-8",
    )

    report = summarize_beast_log(log_path, burnin_fraction=0.25)
    output_path = write_beast_log_summary_table(summary_path, report)

    text = summary_path.read_text(encoding="utf-8")
    assert report.burnin_row_count == 1
    assert report.kept_row_count == 3
    assert report.first_kept_state == 1000
    assert report.last_kept_state == 3000
    assert report.posterior_parameters == ["posterior"]
    assert report.likelihood_parameters == ["likelihood"]
    assert report.prior_parameters == ["prior"]
    assert report.clock_parameters == ["clockRate"]
    assert set(report.tree_parameters) == {"treeHeight", "birthRate"}
    assert report.other_parameters == ["alpha"]
    assert output_path == summary_path
    assert "parameter_category\tparameter\tsample_count" in text
    assert "posterior\tposterior\t3" in text
    assert "tree\tbirthRate\t3" in text


def test_assess_beast_convergence_respects_burnin_fraction() -> None:
    convergence = assess_beast_convergence(
        fixture("example_beast.log"),
        burnin_fraction=0.25,
        ess_threshold=5.0,
        mean_shift_threshold=0.1,
    )

    assert convergence.burnin_fraction == 0.25
    assert convergence.burnin_row_count == 1
    assert convergence.sample_count == 3


def test_parse_beast_log_accepts_native_sample_header(tmp_path: Path) -> None:
    log_path = tmp_path / "native-beast.log"
    log_path.write_text(
        "# BEAST runtime log\n"
        "Sample\tposterior\tprior\tlikelihood\tclockRate\ttree.height\tbirthRate\n"
        "0\t-48.5\t-6.8\t-41.7\t0.0010\t0.11\t1.0\n"
        "20\t-47.1\t-7.1\t-40.0\t0.0012\t0.45\t1.2\n",
        encoding="utf-8",
    )

    report = parse_beast_log(log_path)
    summary = summarize_beast_log(log_path, burnin_fraction=0.5)

    assert report.row_count == 2
    assert report.rows[1].state == 20
    assert summary.kept_row_count == 1
    assert summary.prior_parameters == ["prior"]
    assert set(summary.tree_parameters) == {"tree.height", "birthRate"}


def test_parse_beast_log_accepts_warning_heavy_fixture() -> None:
    report = parse_beast_log(fixture("engine_outputs/beast/log-warning-heavy.log"))

    assert report.row_count == 2
    assert report.columns == ["posterior", "prior", "likelihood", "clockRate"]
    assert report.rows[0].state == 0
    assert report.rows[1].state == 20


def test_parse_beast_log_accepts_sample_header_fixture_with_bom() -> None:
    report = parse_beast_log(fixture("engine_outputs/beast/log-sample-header.log"))

    assert report.row_count == 2
    assert report.columns == ["posterior", "prior", "likelihood", "tree.height"]
    assert report.rows[1].state == 20


def test_parse_beast_log_reports_truncated_fixture_with_structured_error() -> None:
    with pytest.raises(EngineWorkflowError) as error:
        parse_beast_log(fixture("engine_outputs/beast/log-truncated.log"))

    assert error.value.code == "beast_log_missing_parameter_value"
    assert error.value.details["artifact_kind"] == "beast-log"
    assert error.value.details["column"] == "likelihood"
    assert error.value.details["row_number"] == 3
    assert error.value.details["expected_section"] == "sampled parameter row"


def test_summarize_beast_log_reads_real_beast_fixture() -> None:
    report = parse_beast_log(fixture("beast2_strict_yule_posterior.log"))
    summary = summarize_beast_log(
        fixture("beast2_strict_yule_posterior.log"),
        burnin_fraction=0.1,
    )

    assert report.row_count >= 50
    assert report.columns[:3] == ["posterior", "prior", "likelihood"]
    assert "clockRate" in report.columns
    assert summary.kept_row_count >= 40
    assert summary.posterior_parameters == ["posterior"]
    assert {"prior", "treePrior"} <= set(summary.prior_parameters)
    assert {"tree.height", "birthRate"} <= set(summary.tree_parameters)
    assert summary.clock_parameters == ["clockRate"]


def test_summarize_beast_log_reports_known_reference_ess(
    tmp_path: Path,
) -> None:
    log_path = tmp_path / "known-ess.log"
    log_path.write_text(
        "state\tposterior\n0\t0\n1\t0\n2\t1\n3\t1\n",
        encoding="utf-8",
    )

    report = summarize_beast_log(log_path, burnin_fraction=0.0)
    posterior = next(
        summary
        for summary in report.parameter_summaries
        if summary.parameter == "posterior"
    )

    assert posterior.effective_sample_size == pytest.approx(2.666667)


@pytest.mark.parametrize(
    ("burnin_fraction", "burnin_count", "kept_count"),
    [
        (0.0, 0, 101),
        (0.1, 10, 91),
        (0.25, 25, 76),
        (0.5, 50, 51),
    ],
)
def test_real_beast_fixtures_cover_standard_burnin_fractions(
    burnin_fraction: float,
    burnin_count: int,
    kept_count: int,
) -> None:
    log_summary = summarize_beast_log(
        fixture("beast2_strict_yule_posterior.log"),
        burnin_fraction=burnin_fraction,
    )
    tree_report = parse_beast_posterior_tree_samples(
        fixture("beast2_strict_yule_posterior.trees"),
        burnin_fraction=burnin_fraction,
    )

    assert log_summary.burnin_row_count == burnin_count
    assert log_summary.kept_row_count == kept_count
    assert tree_report.burnin_tree_count == burnin_count
    assert tree_report.kept_tree_count == kept_count


def test_summarize_beast_log_reports_median_sd_and_hpd_interval(
    tmp_path: Path,
) -> None:
    log_path = tmp_path / "posterior-diagnostics.log"
    summary_path = tmp_path / "posterior-diagnostics.tsv"
    posterior_values = [*range(19), 100.0]
    log_lines = [
        "# posterior diagnostics fixture",
        "state\tposterior\tclockRate",
    ]
    for state, posterior in enumerate(posterior_values):
        log_lines.append(f"{state}\t{posterior}\t0.001")
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")

    report = summarize_beast_log(log_path, burnin_fraction=0.0)
    output_path = write_beast_log_summary_table(summary_path, report)
    posterior = next(
        summary
        for summary in report.parameter_summaries
        if summary.parameter == "posterior"
    )

    assert posterior.mean == pytest.approx(13.55)
    assert posterior.median == pytest.approx(9.5)
    assert posterior.standard_deviation == pytest.approx(
        round(stdev(posterior_values), 6)
    )
    assert posterior.hpd_95_lower == pytest.approx(0.0)
    assert posterior.hpd_95_upper == pytest.approx(18.0)
    assert output_path == summary_path
    text = summary_path.read_text(encoding="utf-8")
    assert "median\tstandard_deviation" in text
    assert "hpd_95_lower\thpd_95_upper" in text


def test_parse_beast_posterior_tree_samples_handles_translate_and_burnin(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "posterior.trees"
    normalized_path = tmp_path / "posterior.nwk"
    tree_path.write_text(
        "#NEXUS\n"
        "Begin trees;\n"
        "  Translate\n"
        "    1 A,\n"
        "    2 B,\n"
        "    3 C,\n"
        "    4 D\n"
        "  ;\n"
        "tree STATE_0 = [&R] ((1:0.1,2:0.1):0.2,(3:0.1,4:0.1):0.2):0.0;\n"
        "tree STATE_20 = [&R] ((1:0.1,3:0.1):0.2,(2:0.1,4:0.1):0.2):0.0[&lnP=-10.0];\n"
        "tree STATE_40 = [&R] ((1:0.1,2:0.1):0.2,(3:0.1,4:0.1):0.2):0.0;\n"
        "End;\n",
        encoding="utf-8",
    )

    report = parse_beast_posterior_tree_samples(tree_path, burnin_fraction=0.25)
    output_path = write_beast_posterior_tree_set(normalized_path, report)
    _tree, mcc_report = summarize_maximum_clade_credibility_tree(
        output_path, burnin_fraction=0.0
    )
    _consensus_tree, consensus_report = compute_consensus_tree(output_path)

    assert report.total_tree_count == 3
    assert report.burnin_tree_count == 0
    assert report.kept_tree_count == 3
    assert report.rooted_tree_count == 3
    assert report.sampled_states == [0, 20, 40]
    assert report.tip_names == ["A", "B", "C", "D"]
    assert report.trees[1].newick == "((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);"
    assert report.clades[0].clade == "A|B"
    assert output_path == normalized_path
    assert mcc_report.kept_tree_count == 3
    assert consensus_report.tree_count == 3


def test_parse_beast_posterior_tree_samples_reports_metadata_annotations(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "posterior-annotations.trees"
    tree_path.write_text(
        "#NEXUS\n"
        "Begin trees;\n"
        "tree STATE_0 = [&R] ((A[&height=0.1]:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2):0.0[&lnP=-10.0,height_95%_HPD={0.2,0.4}];\n"
        "End;\n",
        encoding="utf-8",
    )

    report = parse_beast_posterior_tree_samples(tree_path, burnin_fraction=0.0)
    sample = report.trees[0]

    assert sample.annotation_key_count == 3
    assert sample.annotation_record_count == 3
    assert sample.annotation_keys == ["height", "height_95%_HPD", "lnP"]
    assert sample.annotation_values["height"] == "0.1"
    assert sample.annotation_values["lnP"] == "-10.0"
    assert sample.annotation_values["height_95%_HPD"] == "{0.2,0.4}"


def test_parse_beast_posterior_tree_samples_reads_annotation_rich_beast_fixture() -> (
    None
):
    report = parse_beast_posterior_tree_samples(
        fixture("beast2_annotated_posterior.trees"),
        burnin_fraction=0.0,
    )

    sample = report.trees[0]
    assert sample.annotation_key_count == 6
    assert sample.annotation_keys == [
        "height",
        "height_95%_HPD",
        "lnP",
        "posterior",
        "rate",
        "rate_95%_HPD",
    ]
    assert sample.annotation_values["posterior"] == "0.98"
    assert sample.annotation_values["rate_95%_HPD"] == "{0.5,0.9}"


def test_parse_beast_posterior_tree_samples_accepts_quoted_translate_labels_fixture() -> (
    None
):
    report = parse_beast_posterior_tree_samples(
        fixture("engine_outputs/beast/posterior-tree-version-variant.trees"),
        burnin_fraction=0.0,
    )

    sample = report.trees[0]
    assert report.sampled_states == [10]
    assert sample.tip_names == ["Taxon A", "Taxon B's sample", "Taxon_C"]
    assert "Taxon A" in sample.newick
    assert "Taxon B''s sample" in sample.newick
    assert sample.annotation_values["region"] == '"North, Sea"'
    assert sample.annotation_values["note"] == "'founder''s clade'"


def test_parse_beast_posterior_tree_samples_reports_invalid_translate_fixture() -> None:
    with pytest.raises(EngineWorkflowError) as error:
        parse_beast_posterior_tree_samples(
            fixture("engine_outputs/beast/posterior-tree-invalid-translate.trees"),
            burnin_fraction=0.0,
        )

    assert error.value.code == "beast_tree_invalid_translate_block"
    assert error.value.details["artifact_kind"] == "beast-posterior-trees"
    assert error.value.details["expected_section"] == "translate block"


def test_parse_beast_posterior_tree_samples_reports_structured_parse_errors(
    tmp_path: Path,
) -> None:
    missing_path = tmp_path / "missing.trees"
    malformed_path = tmp_path / "malformed.trees"
    malformed_path.write_text(
        "#NEXUS\n"
        "Begin trees;\n"
        "tree STATE_0 = ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2;\n"
        "End;\n",
        encoding="utf-8",
    )

    with pytest.raises(EngineWorkflowError) as missing:
        parse_beast_posterior_tree_samples(missing_path, burnin_fraction=0.0)
    assert missing.value.code == "beast_tree_missing_file"
    assert missing.value.details["artifact_kind"] == "beast-posterior-trees"

    with pytest.raises(EngineWorkflowError) as malformed:
        parse_beast_posterior_tree_samples(malformed_path, burnin_fraction=0.0)
    assert malformed.value.code == "beast_tree_parse_error"
    assert malformed.value.details["artifact_kind"] == "beast-posterior-trees"
    assert malformed.value.details["tree_name"] == "STATE_0"


def test_parse_beast_posterior_tree_samples_applies_burnin_fraction(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "posterior.trees"
    tree_path.write_text(
        "#NEXUS\n"
        "Begin trees;\n"
        "tree STATE_0 = (A:0.1,(B:0.1,(C:0.1,D:0.1):0.1):0.1):0.0;\n"
        "tree STATE_10 = ((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1):0.0;\n"
        "tree STATE_20 = ((A:0.1,C:0.1):0.1,(B:0.1,D:0.1):0.1):0.0;\n"
        "tree STATE_30 = ((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1):0.0;\n"
        "End;\n",
        encoding="utf-8",
    )

    report = parse_beast_posterior_tree_samples(tree_path, burnin_fraction=0.25)

    assert report.total_tree_count == 4
    assert report.burnin_tree_count == 1
    assert report.kept_tree_count == 3
    assert report.sampled_states == [10, 20, 30]
    assert report.clades[0].clade == "A|B"
    assert report.clades[0].tree_count == 2


def test_subsample_beast_posterior_tree_set_preserves_state_metadata(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "posterior.trees"
    retained_path = tmp_path / "posterior-subsample.nwk"
    table_path = tmp_path / "posterior-subsample.tsv"
    tree_path.write_text(
        "#NEXUS\n"
        "Begin trees;\n"
        "tree STATE_0 = ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2):0.0;\n"
        "tree STATE_10 = ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2):0.0;\n"
        "tree STATE_20 = ((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3):0.0;\n"
        "tree STATE_30 = ((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.4):0.0;\n"
        "End;\n",
        encoding="utf-8",
    )

    report = subsample_beast_posterior_tree_set(
        tree_path,
        method="evenly-spaced",
        thinning_interval=2,
        burnin_fraction=0.25,
    )
    write_posterior_tree_subsample(retained_path, report)
    write_posterior_tree_subsample_table(table_path, report)

    assert report.burnin_tree_count == 1
    assert report.pre_subsampling_tree_count == 3
    assert report.retained_tree_count == 2
    assert report.retained_source_indices == [2, 4]
    assert [tree.state for tree in report.trees] == [10, 30]
    assert [tree.tree_name for tree in report.trees] == ["STATE_10", "STATE_30"]
    assert retained_path.read_text(encoding="utf-8").count("\n") == 2
    assert "STATE_10\t10" in table_path.read_text(encoding="utf-8")


def test_parse_beast_posterior_tree_samples_reads_real_beast_fixture() -> None:
    report = parse_beast_posterior_tree_samples(
        fixture("beast2_strict_yule_posterior.trees"),
        burnin_fraction=0.1,
    )

    assert report.total_tree_count == 101
    assert report.burnin_tree_count == 10
    assert report.kept_tree_count == 91
    assert report.rooted_tree_count == 91
    assert report.sampled_states[0] == 200
    assert report.sampled_states[-1] == 2000
    assert report.tip_names == ["A", "B", "C", "D"]
    assert report.clades


def test_summarize_beast_posterior_trees_builds_majority_rule_consensus(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "posterior.trees"
    tree_path.write_text(
        "#NEXUS\n"
        "Begin trees;\n"
        "tree STATE_0 = ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2):0.0;\n"
        "tree STATE_10 = ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2):0.0;\n"
        "tree STATE_20 = ((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.4):0.0;\n"
        "tree STATE_30 = ((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3):0.0;\n"
        "End;\n",
        encoding="utf-8",
    )

    consensus_tree, report = summarize_beast_posterior_trees(
        tree_path,
        burnin_fraction=0.25,
    )

    assert report.total_tree_count == 4
    assert report.burnin_tree_count == 1
    assert report.kept_tree_count == 3
    assert report.rooted_topology_count == 2
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.annotated_node_count == 2
    assert report.minimum_posterior_probability == pytest.approx(2 / 3)
    assert report.maximum_posterior_probability == pytest.approx(2 / 3)
    assert report.clade_frequency_count == 4
    assert (
        report.consensus_newick
        == "((A:0.1,B:0.1)0.666666666666667:0.35,(C:0.1,D:0.1)0.666666666666667:0.35);"
    )
    assert report.retained_tree_set_path.read_text(encoding="utf-8").count("\n") == 3
    assert consensus_tree.tip_names == ["A", "B", "C", "D"]


def test_summarize_beast_posterior_trees_reads_real_beast_fixture(
    tmp_path: Path,
) -> None:
    copied_path = tmp_path / "beast2_strict_yule_posterior.trees"
    copied_path.write_text(
        fixture("beast2_strict_yule_posterior.trees").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    _tree, report = summarize_beast_posterior_trees(
        copied_path,
        burnin_fraction=0.1,
    )

    assert report.kept_tree_count == 91
    assert report.rooted_topology_count >= 1
    assert report.annotated_node_count >= 1
    assert report.minimum_posterior_probability is not None
    assert report.maximum_posterior_probability is not None
    assert 0.0 < report.minimum_posterior_probability <= 1.0
    assert 0.0 < report.maximum_posterior_probability <= 1.0
    assert "A" in report.consensus_newick


def test_real_beast_reference_bundle_matches_expected_outputs(tmp_path: Path) -> None:
    reference = SHARED_BEAST_POSTERIOR.load_reference()
    for burnin_fraction, counts in sorted(reference.burnin_reference.items()):
        log_summary = summarize_beast_log(
            SHARED_BEAST_POSTERIOR.posterior_log_path,
            burnin_fraction=burnin_fraction,
        )
        tree_report = parse_beast_posterior_tree_samples(
            SHARED_BEAST_POSTERIOR.posterior_trees_path,
            burnin_fraction=burnin_fraction,
        )
        assert log_summary.burnin_row_count == counts.burnin_row_count
        assert log_summary.kept_row_count == counts.kept_row_count
        assert tree_report.burnin_tree_count == counts.burnin_tree_count
        assert tree_report.kept_tree_count == counts.kept_tree_count

    parameter_reference = reference.parameter_reference
    parameter_summary = summarize_beast_log(
        SHARED_BEAST_POSTERIOR.posterior_log_path,
        burnin_fraction=SHARED_BEAST_POSTERIOR.recommended_burnin_fraction,
    )
    parameters = {row.parameter: row for row in parameter_summary.parameter_summaries}
    for parameter, expected in parameter_reference.items():
        observed = parameters[parameter]
        assert observed.effective_sample_size == pytest.approx(
            expected.effective_sample_size
        )
        assert observed.mean == pytest.approx(expected.mean)
        assert observed.median == pytest.approx(expected.median)
        assert observed.hpd_95_lower == pytest.approx(expected.hpd_95_lower)
        assert observed.hpd_95_upper == pytest.approx(expected.hpd_95_upper)

    copied_path = tmp_path / "beast2_strict_yule_posterior.trees"
    copied_path.write_text(
        SHARED_BEAST_POSTERIOR.posterior_trees_path.read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    _consensus_tree, consensus = summarize_beast_posterior_trees(
        copied_path,
        burnin_fraction=SHARED_BEAST_POSTERIOR.recommended_burnin_fraction,
    )
    _mcc_tree, mcc = summarize_maximum_clade_credibility_tree(
        copied_path,
        burnin_fraction=SHARED_BEAST_POSTERIOR.recommended_burnin_fraction,
    )
    assert (
        consensus.consensus_newick
        == SHARED_BEAST_POSTERIOR.consensus_tree_path.read_text(
            encoding="utf-8"
        ).strip()
    )
    assert consensus.consensus_newick == reference.consensus_reference.newick
    assert (
        consensus.annotated_node_count
        == reference.consensus_reference.annotated_node_count
    )
    assert consensus.minimum_posterior_probability == pytest.approx(
        reference.consensus_reference.minimum_posterior_probability
    )
    assert consensus.maximum_posterior_probability == pytest.approx(
        reference.consensus_reference.maximum_posterior_probability
    )
    assert (
        mcc.mcc_newick
        == SHARED_BEAST_POSTERIOR.mcc_tree_path.read_text(encoding="utf-8").strip()
    )
    assert mcc.mcc_newick == reference.mcc_reference.newick
    assert mcc.selected_tree_index == reference.mcc_reference.selected_tree_index
    assert mcc.clade_credibility_score == pytest.approx(
        reference.mcc_reference.clade_credibility_score
    )


def test_summarize_beast_posterior_topology_diversity_reports_topology_metrics(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "posterior.trees"
    tree_path.write_text(
        "#NEXUS\n"
        "Begin trees;\n"
        "tree STATE_0 = ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2):0.0;\n"
        "tree STATE_10 = ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2):0.0;\n"
        "tree STATE_20 = ((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.4):0.0;\n"
        "tree STATE_30 = ((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3):0.0;\n"
        "End;\n",
        encoding="utf-8",
    )

    report = summarize_beast_posterior_topology_diversity(
        tree_path,
        burnin_fraction=0.25,
    )

    assert report.total_tree_count == 4
    assert report.burnin_tree_count == 1
    assert report.kept_tree_count == 3
    assert report.rooted_topology_count == 2
    assert report.dominant_topology_frequency == pytest.approx(2 / 3)
    assert report.pair_count == 3
    assert report.mean_normalized_robinson_foulds_distance > 0.0
    assert report.maximum_robinson_foulds_distance >= 0
    assert report.unstable_clade_count >= 1
    assert report.retained_tree_set_path.read_text(encoding="utf-8").count("\n") == 3


def test_summarize_beast_posterior_topology_diversity_reads_real_beast_fixture(
    tmp_path: Path,
) -> None:
    copied_path = tmp_path / "beast2_strict_yule_posterior.trees"
    copied_path.write_text(
        fixture("beast2_strict_yule_posterior.trees").read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    report = summarize_beast_posterior_topology_diversity(
        copied_path,
        burnin_fraction=0.1,
    )

    assert report.kept_tree_count == 91
    assert report.rooted_topology_count >= 1
    assert report.dominant_topology_frequency > 0.0
    assert report.effective_topology_count >= 1.0
    assert report.pair_count > 0
    assert report.unstable_clade_count >= 0


def test_validate_beast_posterior_log_reports_missing_columns_and_nonmonotonic_states(
    tmp_path: Path,
) -> None:
    broken_log = tmp_path / "broken-beast.log"
    broken_log.write_text(
        "# BEAST fixture log\n"
        "state\tposterior\tclockRate\n"
        "0\t-510.0\t0.0010\n"
        "0\t-509.0\tbad\n",
        encoding="utf-8",
    )

    report = validate_beast_posterior_log(broken_log)

    assert report.valid is False
    assert report.missing_columns == ["likelihood"]
    assert {issue.code for issue in report.issues} >= {
        "missing-required-column",
        "nonmonotonic-state",
        "invalid-parameter-value",
    }


def test_assess_beast_burnin_sensitivity_reports_tree_and_log_shifts() -> None:
    report = assess_beast_burnin_sensitivity(
        fixture("example_tree_set_left.nwk"),
        log_path=fixture("example_beast.log"),
        burnin_fractions=(0.0, 0.25, 0.5),
    )

    assert [row.burnin_fraction for row in report.slices] == [0.0, 0.25, 0.5]
    assert report.slices[0].posterior_mean == -503.0
    assert report.slices[1].posterior_mean == -500.666667
    assert report.slices[2].tree_height_mean == 13.1


def test_assess_beast_burnin_sensitivity_reports_parameter_and_clade_instability(
    tmp_path: Path,
) -> None:
    posterior_path = tmp_path / "burnin-sensitive.trees"
    log_path = tmp_path / "burnin-sensitive.log"
    posterior_path.write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n"
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n"
        "((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);\n"
        "((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);\n",
        encoding="utf-8",
    )
    lines = [
        "# burnin sensitivity fixture",
        "state\tposterior\tlikelihood\ttreeHeight",
    ]
    for index in range(19):
        lines.append(f"{index * 10}\t0.0\t{-100.0 + index}\t10.0")
    lines.append("190\t100.0\t-81.0\t10.0")
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    report = assess_beast_burnin_sensitivity(
        posterior_path,
        log_path=log_path,
        burnin_fractions=(0.0, 0.95),
    )

    ab_shift = next(shift for shift in report.clade_shifts if shift.clade == "A|B")
    posterior_shift = next(
        shift for shift in report.parameter_shifts if shift.parameter == "posterior"
    )

    assert report.changed_consensus_count == 1
    assert report.unstable_parameter_count >= 1
    assert report.unstable_clade_count >= 1
    assert report.slices[0].clade_frequency_count == 4
    assert report.slices[1].first_kept_state == 190
    assert report.slices[1].posterior_mean == 100.0
    assert ab_shift.crosses_majority_threshold is True
    assert ab_shift.maximum_posterior_probability == pytest.approx(0.5)
    assert posterior_shift.unstable is True
    assert posterior_shift.common_hpd_95_lower is None


def test_assess_beast_chain_mixing_flags_stuck_and_inconsistent_chains(
    tmp_path: Path,
) -> None:
    stable = tmp_path / "stable.log"
    stable.write_text(
        "# BEAST fixture log\n"
        "state\tposterior\tlikelihood\tclockRate\ttreeHeight\n"
        "0\t-510.0\t-490.0\t0.0010\t12.0\n"
        "1000\t-509.0\t-489.0\t0.0011\t12.2\n"
        "2000\t-508.0\t-488.0\t0.0012\t12.3\n"
        "3000\t-507.0\t-487.0\t0.0013\t12.4\n",
        encoding="utf-8",
    )
    stuck = tmp_path / "stuck.log"
    stuck.write_text(
        "# BEAST fixture log\n"
        "state\tposterior\tlikelihood\tclockRate\ttreeHeight\n"
        "0\t-500.0\t-480.0\t0.0010\t15.0\n"
        "1000\t-500.0\t-480.0\t0.0010\t15.0\n"
        "2000\t-500.0\t-480.0\t0.0010\t15.0\n"
        "3000\t-500.0\t-480.0\t0.0010\t15.0\n",
        encoding="utf-8",
    )

    report = assess_beast_chain_mixing(
        [stable, stuck],
        ess_threshold=2.0,
        mean_shift_threshold=0.5,
        cross_chain_mean_shift_threshold=1.0,
    )

    assert report.converged is False
    assert {issue.code for issue in report.issues} >= {
        "stuck-parameter",
        "inconsistent-chains",
    }


def test_render_calibration_audit_report_includes_calibration_and_tip_date_sections(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "calibration-audit.html"

    report = render_calibration_audit_report(
        tree_path=fixture("example_tree_named_clades.nwk"),
        calibration_path=fixture("example_calibrations.tsv"),
        tip_dates_path=fixture("example_tip_dates.tsv"),
        alignment_path=fixture("example_alignment.fasta"),
        out_path=output_path,
    )

    html = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert report.invalid_calibration_count == 0
    assert report.method_tier.tier == "parser-only"
    assert "method-tier" in html
    assert "parser-only" in html
    assert "fossil-calibrations" in html
    assert "impossible-constraints" in html
    assert "tip-dates" in html
    assert "limitations" in html
    assert "limitations" in report.machine_manifest["sections"]


def test_build_bayesian_evidence_package_bundles_inputs_outputs_and_reports(
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "analysis.xml"
    report_path = tmp_path / "calibration-audit.html"
    diagnostic_path = tmp_path / "diagnostics.json"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        config_path,
        tree_path=fixture("example_tree_named_clades.nwk"),
        calibration_path=fixture("example_calibrations.tsv"),
        tip_dates_path=fixture("example_tip_dates.tsv"),
    )
    render_calibration_audit_report(
        tree_path=fixture("example_tree_named_clades.nwk"),
        calibration_path=fixture("example_calibrations.tsv"),
        tip_dates_path=fixture("example_tip_dates.tsv"),
        alignment_path=fixture("example_alignment.fasta"),
        out_path=report_path,
    )
    diagnostic_path.write_text(
        json.dumps({"warning_count": 0}, indent=2) + "\n", encoding="utf-8"
    )

    bundle = build_bayesian_evidence_package(
        bundle_root=tmp_path / "bayesian-bundle",
        input_paths=[
            fixture("example_alignment.fasta"),
            fixture("example_calibrations.tsv"),
            fixture("example_tip_dates.tsv"),
        ],
        config_paths=[config_path],
        tree_paths=[fixture("example_tree_named_clades.nwk")],
        log_paths=[fixture("example_beast.log")],
        diagnostic_paths=[diagnostic_path],
        report_paths=[report_path],
    )

    assert bundle.valid is True
    assert bundle.file_count == 8
    assert bundle.config_count == 1
    assert bundle.tree_count == 1
    assert bundle.log_count == 1
    assert bundle.report_count == 1


def test_render_bayesian_diagnostics_report_includes_log_burnin_mixing_and_calibration_sections(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "bayesian-diagnostics.html"
    second_chain = tmp_path / "chain-2.log"
    second_chain.write_text(
        "# BEAST fixture log\n"
        "state\tposterior\tlikelihood\tclockRate\ttreeHeight\n"
        "0\t-501.0\t-481.0\t0.0010\t13.0\n"
        "1000\t-500.8\t-480.8\t0.0011\t13.1\n"
        "2000\t-500.6\t-480.6\t0.0012\t13.1\n"
        "3000\t-500.5\t-480.5\t0.0011\t13.2\n",
        encoding="utf-8",
    )

    report = render_bayesian_diagnostics_report(
        posterior_tree_path=fixture("example_tree_set_left.nwk"),
        primary_log_path=fixture("example_beast.log"),
        additional_log_paths=[second_chain],
        analysis_xml_path=fixture("beast2_strict_yule_posterior.xml"),
        tree_path=fixture("example_tree_named_clades.nwk"),
        calibration_path=fixture("example_calibrations.tsv"),
        tip_dates_path=fixture("example_tip_dates.tsv"),
        alignment_path=fixture("example_alignment.fasta"),
        out_path=output_path,
        burnin_fractions=(0.0, 0.25, 0.5),
        ess_threshold=2.0,
        mean_shift_threshold=1.0,
        cross_chain_mean_shift_threshold=5.0,
    )

    html = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert report.chain_count == 2
    assert report.method_tier.tier == "parser-only"
    assert "method-tier" in html
    assert "analysis-assumptions" in html
    assert "posterior-log-validation" in html
    assert "burnin-sensitivity" in html
    assert "chain-mixing" in html
    assert "supplementary-diagnostics-table" in html
    assert "methods-summary-text" in html
    assert "fossil-calibrations" in html
    assert "tip-dates" in html


def test_render_bayesian_diagnostics_report_marks_recorded_beast_runs_honestly(
    tmp_path: Path,
) -> None:
    executable = _fake_beast(tmp_path / "beast-fixture")
    analysis_path = tmp_path / "analysis.xml"
    output_path = tmp_path / "bayesian-diagnostics.html"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        analysis_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )
    workflow = run_beast_posterior_inference(
        analysis_path,
        executable=executable,
        seed=9,
    )

    report = render_bayesian_diagnostics_report(
        posterior_tree_path=workflow.output_paths["posterior_trees"],
        primary_log_path=workflow.output_paths["posterior_log"],
        analysis_xml_path=analysis_path,
        out_path=output_path,
        burnin_fractions=(0.0, 0.5),
        ess_threshold=2.0,
        mean_shift_threshold=1.0,
        cross_chain_mean_shift_threshold=5.0,
    )

    html = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert report.method_tier.tier == "parser-only"
    assert "method-tier" in html
    assert "workflow-evidence" in html
    assert "recorded-prior-beast-run" in html
    assert "did not execute BEAST itself" in html
    assert "inferred posterior log" in html


def test_render_bayesian_diagnostics_report_rejects_stale_beast_manifest_claims(
    tmp_path: Path,
) -> None:
    executable = _fake_beast(tmp_path / "beast-fixture")
    analysis_path = tmp_path / "analysis.xml"
    output_path = tmp_path / "bayesian-diagnostics.html"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        analysis_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )
    workflow = run_beast_posterior_inference(
        analysis_path,
        executable=executable,
        seed=9,
    )
    workflow.output_paths["posterior_log"].write_text(
        workflow.output_paths["posterior_log"].read_text(encoding="utf-8")
        + "40\t-117.0\t-78.0\t-39.0\t1.0\t0.012\t0.22\n",
        encoding="utf-8",
    )

    report = render_bayesian_diagnostics_report(
        posterior_tree_path=workflow.output_paths["posterior_trees"],
        primary_log_path=workflow.output_paths["posterior_log"],
        analysis_xml_path=analysis_path,
        out_path=output_path,
        burnin_fractions=(0.0, 0.5),
        ess_threshold=2.0,
        mean_shift_threshold=1.0,
        cross_chain_mean_shift_threshold=5.0,
    )

    html = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert "manifest-mismatch" in html
    assert "treated as parsed existing outputs" in html
    assert "recorded-prior-beast-run" not in html
