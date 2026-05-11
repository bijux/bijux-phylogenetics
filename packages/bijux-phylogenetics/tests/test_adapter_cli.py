from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.cli import main

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_mafft(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

if "--version" in sys.argv:
    print("mafft v7.999", file=sys.stderr)
    raise SystemExit(0)

input_path = Path(sys.argv[-1])
records = []
identifier = None
sequence = []
for raw_line in input_path.read_text(encoding="utf-8").splitlines():
    line = raw_line.strip()
    if not line:
        continue
    if line.startswith(">"):
        if identifier is not None:
            records.append((identifier, "".join(sequence)))
        identifier = line[1:]
        sequence = []
    else:
        sequence.append(line)
if identifier is not None:
    records.append((identifier, "".join(sequence)))
width = max(len(sequence) for _identifier, sequence in records)
for identifier, sequence in records:
    print(f">{identifier}")
    print(sequence.ljust(width, "-"))
print("WARNING: mafft fixture inserted alignment padding", file=sys.stderr)
""",
    )


def _fake_trimal(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

if "--version" in sys.argv:
    print("trimAl v2.0")
    raise SystemExit(0)

args = sys.argv[1:]
input_path = Path(args[args.index("-in") + 1])
output_path = Path(args[args.index("-out") + 1])
if "-strictplus" in args:
    trim_count = 3
    warning = "warning: trimal fixture strictplus trimmed three trailing sites"
elif "-strict" in args:
    trim_count = 2
    warning = "warning: trimal fixture strict trimmed two trailing sites"
elif "-automated1" in args:
    trim_count = 2
    warning = "warning: trimal fixture automated1 trimmed two trailing sites"
elif "-gappyout" in args:
    trim_count = 1
    warning = "warning: trimal fixture gappyout trimmed one trailing site"
else:
    trim_count = 1
    warning = "warning: trimal fixture gap-threshold trimmed one trailing site"
records = []
identifier = None
sequence = []
for raw_line in input_path.read_text(encoding="utf-8").splitlines():
    line = raw_line.strip()
    if not line:
        continue
    if line.startswith(">"):
        if identifier is not None:
            records.append((identifier, "".join(sequence)))
        identifier = line[1:]
        sequence = []
    else:
        sequence.append(line)
if identifier is not None:
    records.append((identifier, "".join(sequence)))
output_path.parent.mkdir(parents=True, exist_ok=True)
with output_path.open("w", encoding="utf-8") as handle:
    for identifier, sequence in records:
        handle.write(f">{identifier}\\n{sequence[:-trim_count]}\\n")
print(warning, file=sys.stderr)
""",
    )


def _fake_iqtree(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("IQ-TREE multicore version 2.9.9")
    raise SystemExit(0)

prefix = Path(args[args.index("-pre") + 1]) if "-pre" in args else Path("iqtree")
prefix.parent.mkdir(parents=True, exist_ok=True)
if "-m" in args and args[args.index("-m") + 1] == "MF":
    prefix.with_suffix(".iqtree").write_text(
        "Best-fit model according to BIC: GTR+G\\nWARNING: model search used a fixture backend\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".model").write_text("Best-fit model: GTR+G\\n", encoding="utf-8")
    print("warning: iqtree fixture model selection", file=sys.stderr)
    raise SystemExit(0)

if "-bb" in args:
    prefix.with_suffix(".treefile").write_text("((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n", encoding="utf-8")
    prefix.with_suffix(".ufboot").write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".iqtree").write_text("Bootstrap analysis completed\\n", encoding="utf-8")
    print("warning: iqtree fixture bootstrap", file=sys.stderr)
    raise SystemExit(0)

prefix.with_suffix(".treefile").write_text("((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n", encoding="utf-8")
prefix.with_suffix(".iqtree").write_text("Tree inference completed\\n", encoding="utf-8")
print("warning: iqtree fixture tree inference", file=sys.stderr)
""",
    )


def _fake_fasttree(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys

args = sys.argv[1:]
if not args or "-help" in args:
    print("FastTree Version 2.2 fixture")
    raise SystemExit(0)

print("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
print("warning: fasttree fixture approximate support only", file=sys.stderr)
""",
    )


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


def test_adapter_inspect_cli_reports_engine_version(tmp_path: Path, capsys) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")

    exit_code = main(
        ["adapter", "inspect", "iqtree", "--executable", str(executable), "--json"]
    )

    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["text"] == "IQ-TREE multicore version 2.9.9"


def test_adapter_align_cli_writes_aligned_fasta_and_manifest(
    tmp_path: Path, capsys
) -> None:
    executable = _fake_mafft(tmp_path / "mafft-fixture")
    input_path = tmp_path / "unaligned.fasta"
    input_path.write_text(">A\nACTG\n>B\nACTGA\n>C\nACT\n", encoding="utf-8")
    output_path = tmp_path / "aligned.fasta"

    exit_code = main(
        [
            "adapter",
            "align",
            str(input_path),
            "--out",
            str(output_path),
            "--executable",
            str(executable),
            "--json",
        ]
    )

    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["warnings"] == ["WARNING: mafft fixture inserted alignment padding"]
    assert output_path.exists()
    assert Path(payload["data"]["manifest_path"]).exists()


def test_adapter_align_cli_passes_named_mafft_mode_to_workflow(
    tmp_path: Path, capsys
) -> None:
    executable = _fake_mafft(tmp_path / "mafft-fixture")
    input_path = tmp_path / "unaligned-protein.fasta"
    input_path.write_text(">A\nMKTW\n>B\nMKTWA\n>C\nMKT\n", encoding="utf-8")
    output_path = tmp_path / "aligned-protein.fasta"

    exit_code = main(
        [
            "adapter",
            "align",
            str(input_path),
            "--out",
            str(output_path),
            "--executable",
            str(executable),
            "--mode",
            "linsi",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["mode"] == "linsi"
    assert payload["data"]["run"]["command"][1:-1] == [
        "--localpair",
        "--maxiterate",
        "1000",
    ]
    assert payload["data"]["notes"][0] == "mafft alignment mode: linsi"
    assert output_path.exists()


def test_adapter_trim_cli_reports_retained_and_removed_sites(
    tmp_path: Path, capsys
) -> None:
    executable = _fake_trimal(tmp_path / "trimal-fixture")
    input_path = fixture("alignments/example_alignment_trim.fasta")
    output_path = tmp_path / "trimmed.fasta"

    exit_code = main(
        [
            "adapter",
            "trim",
            str(input_path),
            "--out",
            str(output_path),
            "--executable",
            str(executable),
            "--gap-threshold",
            "0.2",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["mode"] == "gap-threshold"
    assert payload["metrics"]["retained_site_count"] == 5
    assert payload["metrics"]["removed_site_count"] == 1
    assert payload["data"]["trimming_summary"]["gap_threshold"] == 0.2
    assert payload["data"]["notes"][1] == "retained sites: 5 of 6"
    assert output_path.exists()


def test_adapter_trim_cli_passes_named_trimal_mode_to_workflow(
    tmp_path: Path, capsys
) -> None:
    executable = _fake_trimal(tmp_path / "trimal-fixture")
    input_path = fixture("alignments/example_alignment_trim.fasta")
    output_path = tmp_path / "trimmed-strictplus.fasta"

    exit_code = main(
        [
            "adapter",
            "trim",
            str(input_path),
            "--out",
            str(output_path),
            "--executable",
            str(executable),
            "--mode",
            "strictplus",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["mode"] == "strictplus"
    assert payload["data"]["run"]["command"][5:] == ["-strictplus"]
    assert payload["data"]["trimming_summary"]["removed_site_count"] == 3
    assert payload["warnings"] == [
        "warning: trimal fixture strictplus trimmed three trailing sites"
    ]
    assert output_path.exists()


def test_adapter_model_select_and_compare_cli_produce_outputs(
    tmp_path: Path, capsys
) -> None:
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    fasttree = _fake_fasttree(tmp_path / "FastTree-fixture")
    input_path = fixture("alignments/example_alignment.fasta")
    ml_dir = tmp_path / "ml"
    comparison_path = tmp_path / "comparison.html"
    fast_tree_path = tmp_path / "fasttree.nwk"

    model_exit = main(
        [
            "adapter",
            "model-select",
            str(input_path),
            "--out-dir",
            str(tmp_path / "model"),
            "--prefix",
            "example",
            "--executable",
            str(iqtree),
            "--json",
        ]
    )
    model_payload = json.loads(capsys.readouterr().out)
    assert model_exit == 0
    assert model_payload["metrics"]["selected_model"] == "GTR+G"

    ml_exit = main(
        [
            "adapter",
            "infer-ml",
            str(input_path),
            "--out-dir",
            str(ml_dir),
            "--model",
            "GTR+G",
            "--prefix",
            "example",
            "--executable",
            str(iqtree),
            "--json",
        ]
    )
    ml_payload = json.loads(capsys.readouterr().out)
    assert ml_exit == 0
    ml_tree_path = Path(ml_payload["data"]["output_paths"]["tree"])

    fast_exit = main(
        [
            "adapter",
            "infer-fast",
            str(input_path),
            "--out",
            str(fast_tree_path),
            "--executable",
            str(fasttree),
            "--json",
        ]
    )
    fast_payload = json.loads(capsys.readouterr().out)
    assert fast_exit == 0
    assert fast_payload["warnings"] == [
        "warning: fasttree fixture approximate support only"
    ]

    compare_exit = main(
        [
            "adapter",
            "compare",
            "--fast-tree",
            str(fast_tree_path),
            "--ml-tree",
            str(ml_tree_path),
            "--out",
            str(comparison_path),
            "--json",
        ]
    )
    compare_payload = json.loads(capsys.readouterr().out)
    assert compare_exit == 0
    assert compare_payload["metrics"]["shared_taxa"] == 4
    assert comparison_path.exists()


def test_adapter_fasta_to_tree_cli_materializes_pipeline_outputs(
    tmp_path: Path, capsys
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_sequences_raw.fasta")
    out_dir = tmp_path / "fasta-to-tree"

    exit_code = main(
        [
            "adapter",
            "fasta-to-tree",
            str(input_path),
            "--out-dir",
            str(out_dir),
            "--prefix",
            "example",
            "--mafft-executable",
            str(mafft),
            "--trimal-executable",
            str(trimal),
            "--iqtree-executable",
            str(iqtree),
            "--bootstrap-replicates",
            "200",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["selected_model"] == "GTR+G"
    assert payload["metrics"]["sequence_type"] == "dna"
    assert (
        "warning: trimal fixture gap-threshold trimmed one trailing site"
        in payload["warnings"]
    )
    assert Path(payload["data"]["engine_artifact_dir"]).name == "example"
    assert (
        Path(payload["data"]["engine_artifact_dir"]).parent.name == "engine-artifacts"
    )
    outputs = {Path(path).name for path in payload["outputs"]}
    assert {
        "example",
        "example.aln",
        "example.trimmed.aln",
        "example.tree",
        "example.log",
        "example.model.tsv",
        "example.support.tsv",
        "example.manifest.json",
    }.issubset(outputs)
    assert Path(payload["data"]["output_paths"]["tree"]).exists()
    assert Path(payload["data"]["manifest_path"]).exists()


def test_adapter_fasta_to_tree_cli_passes_named_mafft_mode_to_alignment_step(
    tmp_path: Path, capsys
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_sequences_raw.fasta")
    out_dir = tmp_path / "fasta-to-tree-linsi"

    exit_code = main(
        [
            "adapter",
            "fasta-to-tree",
            str(input_path),
            "--out-dir",
            str(out_dir),
            "--prefix",
            "example",
            "--mafft-executable",
            str(mafft),
            "--alignment-mode",
            "linsi",
            "--trimal-executable",
            str(trimal),
            "--iqtree-executable",
            str(iqtree),
            "--bootstrap-replicates",
            "200",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["alignment_mode"] == "linsi"
    assert payload["data"]["alignment_workflow"]["run"]["command"][1:-1] == [
        "--localpair",
        "--maxiterate",
        "1000",
    ]
    assert "mafft alignment mode: linsi" in payload["data"]["notes"]


def test_adapter_fasta_to_tree_cli_passes_named_trimal_mode_to_trimming_step(
    tmp_path: Path, capsys
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_sequences_raw.fasta")
    out_dir = tmp_path / "fasta-to-tree-strictplus"

    exit_code = main(
        [
            "adapter",
            "fasta-to-tree",
            str(input_path),
            "--out-dir",
            str(out_dir),
            "--prefix",
            "example",
            "--mafft-executable",
            str(mafft),
            "--trimal-executable",
            str(trimal),
            "--trimming-mode",
            "strictplus",
            "--iqtree-executable",
            str(iqtree),
            "--bootstrap-replicates",
            "200",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["trimming_mode"] == "strictplus"
    assert payload["metrics"]["removed_site_count"] == 3
    assert payload["data"]["trimming_workflow"]["run"]["command"][5:] == [
        "-strictplus"
    ]
    assert payload["data"]["trimming_workflow"]["trimming_summary"][
        "removed_site_count"
    ] == 3
    assert "trimal trimming mode: strictplus" in payload["data"]["notes"]


def test_adapter_fasta_to_tree_cli_repairs_invalid_input_when_requested(
    tmp_path: Path, capsys
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_sequences_invalid_input.fasta")
    out_dir = tmp_path / "fasta-to-tree-repaired"

    exit_code = main(
        [
            "adapter",
            "fasta-to-tree",
            str(input_path),
            "--out-dir",
            str(out_dir),
            "--prefix",
            "example",
            "--sequence-type",
            "dna",
            "--mafft-executable",
            str(mafft),
            "--trimal-executable",
            str(trimal),
            "--iqtree-executable",
            str(iqtree),
            "--normalize-identifiers",
            "--remove-invalid-records",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["normalized_identifier_count"] == 2
    assert payload["metrics"]["removed_record_count"] == 2
    prepared_input_path = Path(payload["data"]["prepared_input_path"])
    assert prepared_input_path.read_text(encoding="utf-8") == (
        ">Alpha_sample\nACTGACTG\n"
        ">rare_taxon\nACTGACTGACTGACTGACTGACTG\n"
    )


def test_adapter_fasta_to_tree_cli_rejects_mixed_input_without_explicit_type(
    tmp_path: Path, capsys
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    mixed_input = tmp_path / "mixed.fasta"
    mixed_input.write_text(
        ">dna_like\nACTGACTG\n>rna_like\nACUGACUG\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "adapter",
            "fasta-to-tree",
            str(mixed_input),
            "--out-dir",
            str(tmp_path / "mixed-output"),
            "--prefix",
            "mixed",
            "--mafft-executable",
            str(mafft),
            "--trimal-executable",
            str(trimal),
            "--iqtree-executable",
            str(iqtree),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert "thymine-bearing and uracil-bearing" in payload["errors"][0]["message"]


def test_adapter_fasta_to_tree_cli_can_force_declared_type_on_mixed_input(
    tmp_path: Path, capsys
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    mixed_input = tmp_path / "mixed.fasta"
    mixed_input.write_text(
        ">dna_a\nACTGACTG\n>dna_b\nACTGACTA\n>rna_like\nACUGACUG\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "adapter",
            "fasta-to-tree",
            str(mixed_input),
            "--out-dir",
            str(tmp_path / "forced-dna"),
            "--prefix",
            "mixed",
            "--sequence-type",
            "dna",
            "--mafft-executable",
            str(mafft),
            "--trimal-executable",
            str(trimal),
            "--iqtree-executable",
            str(iqtree),
            "--remove-invalid-records",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["sequence_type"] == "dna"
    assert payload["metrics"]["sequence_type_confidence"] == "medium"
    prepared_input_path = Path(payload["data"]["prepared_input_path"])
    assert prepared_input_path.read_text(encoding="utf-8") == (
        ">dna_a\nACTGACTG\n>dna_b\nACTGACTA\n"
    )


def test_adapter_mrbayes_cli_and_engine_report(tmp_path: Path, capsys) -> None:
    executable = _fake_mrbayes(tmp_path / "mb-fixture")
    alignment_path = fixture("alignments/example_alignment.fasta")
    nexus_path = tmp_path / "analysis.nex"

    prepare_exit = main(
        [
            "adapter",
            "mrbayes-prepare",
            str(alignment_path),
            "--out",
            str(nexus_path),
            "--ngen",
            "2000",
            "--samplefreq",
            "50",
            "--printfreq",
            "50",
            "--json",
        ]
    )
    prepare_payload = json.loads(capsys.readouterr().out)
    assert prepare_exit == 0
    assert prepare_payload["metrics"]["taxon_count"] == 4

    run_exit = main(
        [
            "adapter",
            "mrbayes-run",
            str(nexus_path),
            "--executable",
            str(executable),
            "--json",
        ]
    )
    run_payload = json.loads(capsys.readouterr().out)
    assert run_exit == 0
    manifest_path = Path(run_payload["data"]["manifest_path"])
    tree_path = Path(run_payload["data"]["output_paths"]["posterior_trees"])
    trace_path = Path(run_payload["data"]["output_paths"]["parameter_traces"])
    assert run_payload["warnings"] == ["warning: mrbayes fixture posterior run"]

    summarize_exit = main(
        [
            "adapter",
            "mrbayes-summarize",
            str(tree_path),
            "--burnin-fraction",
            "0.25",
            "--json",
        ]
    )
    summarize_payload = json.loads(capsys.readouterr().out)
    assert summarize_exit == 0
    assert summarize_payload["metrics"]["kept_tree_count"] == 3

    traces_exit = main(["adapter", "mrbayes-traces", str(trace_path), "--json"])
    traces_payload = json.loads(capsys.readouterr().out)
    assert traces_exit == 0
    assert traces_payload["metrics"]["row_count"] == 4

    ess_exit = main(["adapter", "mrbayes-ess", str(trace_path), "--json"])
    ess_payload = json.loads(capsys.readouterr().out)
    assert ess_exit == 0
    assert ess_payload["metrics"]["parameter_count"] == 3

    report_path = tmp_path / "inference-report.html"
    report_exit = main(
        [
            "adapter",
            "report",
            str(manifest_path),
            "--out",
            str(report_path),
            "--json",
        ]
    )
    report_payload = json.loads(capsys.readouterr().out)
    assert report_exit == 0
    assert report_payload["metrics"]["warning_count"] == 1
    assert report_path.exists()


def test_adapter_mrbayes_convergence_and_posterior_report_cli_emit_metrics(
    tmp_path: Path, capsys
) -> None:
    executable = _fake_mrbayes(tmp_path / "mb-fixture")
    alignment_path = fixture("alignments/example_alignment.fasta")
    nexus_path = tmp_path / "analysis.nex"
    posterior_report_path = tmp_path / "posterior-report.html"

    assert (
        main(
            [
                "adapter",
                "mrbayes-prepare",
                str(alignment_path),
                "--out",
                str(nexus_path),
                "--json",
            ]
        )
        == 0
    )
    capsys.readouterr()

    run_exit = main(
        [
            "adapter",
            "mrbayes-run",
            str(nexus_path),
            "--executable",
            str(executable),
            "--json",
        ]
    )
    run_payload = json.loads(capsys.readouterr().out)
    assert run_exit == 0
    trace_path = Path(run_payload["data"]["output_paths"]["parameter_traces"])
    tree_path = Path(run_payload["data"]["output_paths"]["posterior_trees"])

    convergence_exit = main(
        [
            "adapter",
            "mrbayes-convergence",
            str(trace_path),
            "--ess-threshold",
            "5",
            "--mean-shift-threshold",
            "0.1",
            "--json",
        ]
    )
    convergence_payload = json.loads(capsys.readouterr().out)
    assert convergence_exit == 0
    assert convergence_payload["metrics"]["warning_count"] >= 1

    report_exit = main(
        [
            "adapter",
            "mrbayes-report",
            str(tree_path),
            "--traces",
            str(trace_path),
            "--out",
            str(posterior_report_path),
            "--ess-threshold",
            "5",
            "--mean-shift-threshold",
            "0.1",
            "--json",
        ]
    )
    report_payload = json.loads(capsys.readouterr().out)
    assert report_exit == 0
    assert report_payload["metrics"]["kept_tree_count"] == 3
    assert posterior_report_path.exists()


def test_adapter_beast_surface_and_bayesian_evidence_cli_write_outputs(
    tmp_path: Path, capsys
) -> None:
    analysis_path = tmp_path / "analysis.xml"
    report_path = tmp_path / "calibration-audit.html"
    diagnostics_path = tmp_path / "diagnostics.json"
    bundle_root = tmp_path / "bayesian-bundle"

    prepare_exit = main(
        [
            "adapter",
            "beast-prepare",
            str(fixture("alignments/example_alignment.fasta")),
            "--out",
            str(analysis_path),
            "--tree",
            str(fixture("trees/example_tree_named_clades.nwk")),
            "--calibrations",
            str(fixture("metadata/example_calibrations.tsv")),
            "--tip-dates",
            str(fixture("metadata/example_tip_dates.tsv")),
            "--clock-model",
            "relaxed-lognormal",
            "--tree-prior",
            "birth-death",
            "--json",
        ]
    )
    prepare_payload = json.loads(capsys.readouterr().out)
    assert prepare_exit == 0
    assert prepare_payload["metrics"]["calibration_count"] == 2

    calibrations_exit = main(
        [
            "adapter",
            "beast-calibrations",
            str(fixture("trees/example_tree_named_clades.nwk")),
            str(fixture("metadata/example_calibrations.tsv")),
            "--json",
        ]
    )
    calibrations_payload = json.loads(capsys.readouterr().out)
    assert calibrations_exit == 0
    assert calibrations_payload["metrics"]["invalid_calibration_count"] == 0

    tip_dates_exit = main(
        [
            "adapter",
            "beast-tip-dates",
            str(fixture("trees/example_tree_named_clades.nwk")),
            str(fixture("metadata/example_tip_dates.tsv")),
            "--alignment",
            str(fixture("alignments/example_alignment.fasta")),
            "--json",
        ]
    )
    tip_dates_payload = json.loads(capsys.readouterr().out)
    assert tip_dates_exit == 0
    assert tip_dates_payload["metrics"]["valid_tip_count"] == 4

    log_exit = main(
        ["adapter", "beast-log", str(fixture("metadata/example_beast.log")), "--json"]
    )
    log_payload = json.loads(capsys.readouterr().out)
    assert log_exit == 0
    assert log_payload["metrics"]["row_count"] == 4

    convergence_exit = main(
        [
            "adapter",
            "beast-convergence",
            str(fixture("metadata/example_beast.log")),
            "--ess-threshold",
            "5",
            "--mean-shift-threshold",
            "0.1",
            "--json",
        ]
    )
    convergence_payload = json.loads(capsys.readouterr().out)
    assert convergence_exit == 0
    assert convergence_payload["metrics"]["warning_count"] >= 1

    report_exit = main(
        [
            "adapter",
            "beast-calibration-report",
            str(fixture("trees/example_tree_named_clades.nwk")),
            str(fixture("metadata/example_calibrations.tsv")),
            "--tip-dates",
            str(fixture("metadata/example_tip_dates.tsv")),
            "--alignment",
            str(fixture("alignments/example_alignment.fasta")),
            "--out",
            str(report_path),
            "--json",
        ]
    )
    report_payload = json.loads(capsys.readouterr().out)
    assert report_exit == 0
    assert report_payload["metrics"]["invalid_calibration_count"] == 0
    assert report_path.exists()

    diagnostics_path.write_text(
        json.dumps({"warning_count": 0}, indent=2) + "\n", encoding="utf-8"
    )
    evidence_exit = main(
        [
            "adapter",
            "bayesian-evidence",
            "--out-dir",
            str(bundle_root),
            "--inputs",
            str(fixture("alignments/example_alignment.fasta")),
            str(fixture("metadata/example_calibrations.tsv")),
            str(fixture("metadata/example_tip_dates.tsv")),
            "--configs",
            str(analysis_path),
            "--trees",
            str(fixture("trees/example_tree_named_clades.nwk")),
            "--logs",
            str(fixture("metadata/example_beast.log")),
            "--diagnostics",
            str(diagnostics_path),
            "--reports",
            str(report_path),
            "--json",
        ]
    )
    evidence_payload = json.loads(capsys.readouterr().out)
    assert evidence_exit == 0
    assert evidence_payload["metrics"]["valid"] is True
    assert bundle_root.exists()


def test_tree_set_uncertainty_cli_surfaces_modes_conflicts_and_package(
    tmp_path: Path, capsys
) -> None:
    package_dir = tmp_path / "posterior-uncertainty-package"

    diversity_exit = main(
        [
            "tree-set",
            "diversity-compare",
            str(fixture("trees/example_tree_set_left.nwk")),
            str(fixture("trees/example_tree_set_right.nwk")),
            "--json",
        ]
    )
    diversity_payload = json.loads(capsys.readouterr().out)
    assert diversity_exit == 0
    assert diversity_payload["metrics"]["left_rooted_topology_count"] == 2

    multimodality_exit = main(
        [
            "tree-set",
            "multimodality",
            str(fixture("trees/example_tree_set_left.nwk")),
            "--min-mode-frequency",
            "0.3",
            "--json",
        ]
    )
    multimodality_payload = json.loads(capsys.readouterr().out)
    assert multimodality_exit == 0
    assert multimodality_payload["metrics"]["multimodal"] is True

    conflicts_exit = main(
        [
            "tree-set",
            "clade-conflicts",
            str(fixture("trees/example_tree_set_left.nwk")),
            "--credibility-threshold",
            "0.3",
            "--json",
        ]
    )
    conflicts_payload = json.loads(capsys.readouterr().out)
    assert conflicts_exit == 0
    assert conflicts_payload["metrics"]["conflict_count"] == 2

    package_exit = main(
        [
            "tree-set",
            "package",
            str(fixture("trees/example_tree_set_left.nwk")),
            "--out-dir",
            str(package_dir),
            "--json",
        ]
    )
    package_payload = json.loads(capsys.readouterr().out)
    assert package_exit == 0
    assert package_payload["metrics"]["artifact_count"] == 7
    assert (package_dir / "consensus-tree.svg").exists()
    assert (package_dir / "clade-frequency-plot.svg").exists()


def test_adapter_bayesian_uncertainty_cli_writes_table_and_methods_text(
    tmp_path: Path, capsys
) -> None:
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
    table_path = tmp_path / "bayesian-diagnostics.tsv"
    methods_path = tmp_path / "bayesian-methods.md"

    table_exit = main(
        [
            "adapter",
            "bayesian-diagnostics-table",
            str(fixture("trees/example_tree_set_left.nwk")),
            "--log",
            str(fixture("metadata/example_beast.log")),
            "--additional-logs",
            str(second_chain),
            "--out",
            str(table_path),
            "--burnin-fractions",
            "0.0",
            "0.25",
            "--ess-threshold",
            "2",
            "--mean-shift-threshold",
            "1.0",
            "--cross-chain-mean-shift-threshold",
            "5.0",
            "--json",
        ]
    )
    table_payload = json.loads(capsys.readouterr().out)
    assert table_exit == 0
    assert table_payload["metrics"]["row_count"] >= 1
    assert table_path.exists()

    methods_exit = main(
        [
            "adapter",
            "bayesian-methods",
            str(fixture("trees/example_tree_set_left.nwk")),
            "--log",
            str(fixture("metadata/example_beast.log")),
            "--additional-logs",
            str(second_chain),
            "--out",
            str(methods_path),
            "--tree-prior",
            "birth-death",
            "--clock-model",
            "relaxed-lognormal",
            "--calibration-path",
            str(fixture("metadata/example_calibrations.tsv")),
            "--tip-dates-path",
            str(fixture("metadata/example_tip_dates.tsv")),
            "--burnin-fractions",
            "0.0",
            "0.25",
            "--ess-threshold",
            "2",
            "--mean-shift-threshold",
            "1.0",
            "--cross-chain-mean-shift-threshold",
            "5.0",
            "--json",
        ]
    )
    methods_payload = json.loads(capsys.readouterr().out)
    assert methods_exit == 0
    assert methods_payload["metrics"]["warning_count"] >= 0
    assert "relaxed-lognormal" in methods_path.read_text(encoding="utf-8")
