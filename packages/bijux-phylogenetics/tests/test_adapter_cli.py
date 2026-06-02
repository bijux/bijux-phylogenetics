from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main

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
    is_protein = "-st" in args and args[args.index("-st") + 1] == "AA"
    selected_model = "LG+G4" if is_protein else "GTR+G"
    criteria_lines = (
        " No. Model         -LnL         df  AIC          AICc         BIC\\n"
        + (
            "  1  LG+G4         120.100      9   258.200      310.200      260.200\\n"
            "  2  LG+I+G4       119.900      10  259.800      323.800      262.000\\n"
            "  3  WAG+G4        121.500      9   261.000      313.000      263.000\\n"
            "Akaike Information Criterion:           LG+G4\\n"
            "Corrected Akaike Information Criterion: LG+G4\\n"
            "Bayesian Information Criterion:         LG+G4\\n"
            "Best-fit model: LG+G4 chosen according to BIC\\n"
            if is_protein
            else
            "  1  GTR+G         123.456      12  270.912      330.912      272.912\\n"
            "  2  HKY+G         124.000      10  268.000      320.000      269.000\\n"
            "  3  JC            130.500      5   271.000      300.000      271.500\\n"
            "Akaike Information Criterion:           HKY+G\\n"
            "Corrected Akaike Information Criterion: JC\\n"
            "Bayesian Information Criterion:         GTR+G\\n"
            "Best-fit model according to BIC: GTR+G\\n"
        )
    )
    prefix.with_suffix(".iqtree").write_text(
        criteria_lines
        + f"Log-likelihood of the tree: -123.456\\nWARNING: model search used a fixture backend\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture model-selection log\\nBEST SCORE FOUND : -123.456\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".model").write_text(
        f"Best-fit model: {selected_model}\\n",
        encoding="utf-8",
    )
    print("warning: iqtree fixture model selection", file=sys.stderr)
    raise SystemExit(0)

if "-alrt" in args:
    prefix.with_suffix(".treefile").write_text("((A:0.1,B:0.1)82/97:0.2,(C:0.1,D:0.1)79/96:0.2);\\n", encoding="utf-8")
    prefix.with_suffix(".ufboot").write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".iqtree").write_text(
        "Best-fit model: GTR+G\\nLog-likelihood of the tree: -222.222\\nSH-aLRT and ultrafast bootstrap analysis completed\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture sh-alrt log\\nBEST SCORE FOUND : -222.222\\n",
        encoding="utf-8",
    )
    print("warning: iqtree fixture sh-alrt", file=sys.stderr)
    raise SystemExit(0)

if "-bb" in args:
    prefix.with_suffix(".contree").write_text("((A:0.1,B:0.1)95:0.2,(C:0.1,D:0.1)88:0.2);\\n", encoding="utf-8")
    prefix.with_suffix(".ufboot").write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".iqtree").write_text(
        "Best-fit model: GTR+G\\nLog-likelihood of the tree: -234.567\\nBootstrap analysis completed\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture bootstrap log\\nBEST SCORE FOUND : -234.567\\n",
        encoding="utf-8",
    )
    print("warning: iqtree fixture bootstrap", file=sys.stderr)
    raise SystemExit(0)

prefix.with_suffix(".treefile").write_text("((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n", encoding="utf-8")
prefix.with_suffix(".iqtree").write_text(
    "Best-fit model: GTR+G\\nLog-likelihood of the tree: -345.678\\nTree inference completed\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture inference log\\nBEST SCORE FOUND : -345.678\\n",
    encoding="utf-8",
)
print("warning: iqtree fixture tree inference", file=sys.stderr)
""",
    )


def _unsupported_iqtree(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys

if "--version" in sys.argv:
    print("IQ-TREE multicore version 1.6.0")
    raise SystemExit(0)

raise SystemExit(1)
""",
    )


def _fake_iqtree_without_support_labels(path: Path) -> Path:
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
prefix.with_suffix(".contree").write_text(
    "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
    encoding="utf-8",
)
prefix.with_suffix(".ufboot").write_text(
    "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n"
    "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
    encoding="utf-8",
)
prefix.with_suffix(".iqtree").write_text(
    "Best-fit model: GTR+G\\nLog-likelihood of the tree: -222.222\\nSupport analysis completed\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture support log\\nBEST SCORE FOUND : -222.222\\n",
    encoding="utf-8",
)
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

print("((A:0.1,B:0.1)0.98:0.3,(C:0.1,D:0.1)0.62:0.3);")
print("warning: fasttree fixture approximate support only", file=sys.stderr)
""",
    )


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


def test_adapter_align_cli_can_run_codon_aware_alignment(
    tmp_path: Path, capsys
) -> None:
    executable = _fake_mafft(tmp_path / "mafft-fixture")
    input_path = tmp_path / "coding-raw.fasta"
    input_path.write_text(
        ">short_good\nATGGAATGG\n"
        ">long_good\nATGGAATGGAAA\n"
        ">frameshift\nATGGAATG\n"
        ">internal_stop\nATGTAGTGG\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "codon-aligned.fasta"

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
            "--codon-aware",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["codon_aware"] is True
    assert payload["metrics"]["sequence_type"] == "dna"
    assert payload["metrics"]["genetic_code_id"] == 1
    assert payload["metrics"]["accepted_sequence_count"] == 2
    assert payload["metrics"]["excluded_sequence_count"] == 2
    assert payload["metrics"]["invalid_codon_sequence_count"] == 0
    assert output_path.read_text(encoding="utf-8") == (
        ">long_good\nATGGAATGGAAA\n>short_good\nATGGAATGG---\n"
    )
    assert Path(payload["data"]["output_paths"]["guide_input"]).exists()
    assert Path(payload["data"]["output_paths"]["excluded_sequences"]).exists()
    assert Path(payload["data"]["output_paths"]["coding_summary"]).exists()


def test_adapter_align_cli_reports_codon_aware_failures(tmp_path: Path, capsys) -> None:
    executable = _fake_mafft(tmp_path / "mafft-fixture")
    input_path = tmp_path / "coding-invalid.fasta"
    input_path.write_text(
        ">frameshift\nATGGAATG\n>internal_stop\nATGTAGTGG\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "codon-aligned.fasta"

    exit_code = main(
        [
            "adapter",
            "align",
            str(input_path),
            "--out",
            str(output_path),
            "--executable",
            str(executable),
            "--codon-aware",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert "excluded every sequence" in payload["errors"][0]["message"]


def test_adapter_align_cli_explains_invalid_fasta_records(
    tmp_path: Path, capsys
) -> None:
    executable = _fake_mafft(tmp_path / "mafft-fixture")
    output_path = tmp_path / "aligned.fasta"

    exit_code = main(
        [
            "adapter",
            "align",
            str(fixture("alignments/example_sequences_invalid_input.fasta")),
            "--out",
            str(output_path),
            "--executable",
            str(executable),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "invalid_alignment_error"
    assert payload["errors"][0]["details"]["failure_reason"] == "invalid_fasta_input"
    assert (
        payload["errors"][0]["details"]["evidence"]["duplicate_identifier_count"] == 1
    )
    assert payload["errors"][0]["details"]["evidence"]["illegal_character_count"] == 1
    assert payload["errors"][0]["details"]["evidence"]["empty_sequence_count"] == 1


def test_adapter_align_cli_honors_configurable_genetic_code(
    tmp_path: Path, capsys
) -> None:
    executable = _fake_mafft(tmp_path / "mafft-fixture")
    input_path = tmp_path / "coding-mito.fasta"
    input_path.write_text(
        ">shared_good\nATGGAATGG\n>mito_triplet\nATGTGAGGG\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "codon-aligned.fasta"

    exit_code = main(
        [
            "adapter",
            "align",
            str(input_path),
            "--out",
            str(output_path),
            "--executable",
            str(executable),
            "--codon-aware",
            "--genetic-code",
            "2",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["genetic_code_id"] == 2
    assert payload["metrics"]["accepted_sequence_count"] == 2
    assert payload["metrics"]["excluded_sequence_count"] == 0


def test_adapter_align_cli_supports_resume_timeout_and_incomplete_policy(
    tmp_path: Path, capsys
) -> None:
    executable = _fake_mafft(tmp_path / "mafft-fixture")
    input_path = tmp_path / "unaligned.fasta"
    input_path.write_text(">A\nACTG\n>B\nACTGA\n>C\nACT\n", encoding="utf-8")
    output_path = tmp_path / "aligned.fasta"

    first_exit = main(
        [
            "adapter",
            "align",
            str(input_path),
            "--out",
            str(output_path),
            "--executable",
            str(executable),
            "--timeout-seconds",
            "30",
            "--incomplete-run-policy",
            "clean",
            "--json",
        ]
    )
    first_payload = json.loads(capsys.readouterr().out)
    assert first_exit == 0
    assert first_payload["metrics"]["resumed"] is False
    assert first_payload["metrics"]["timeout_seconds"] == 30.0

    second_exit = main(
        [
            "adapter",
            "align",
            str(input_path),
            "--out",
            str(output_path),
            "--executable",
            str(executable),
            "--resume",
            "--timeout-seconds",
            "30",
            "--incomplete-run-policy",
            "clean",
            "--json",
        ]
    )
    second_payload = json.loads(capsys.readouterr().out)
    assert second_exit == 0
    assert second_payload["metrics"]["resumed"] is True
    assert second_payload["metrics"]["timeout_seconds"] == 30.0
    assert Path(second_payload["data"]["manifest_path"]).exists()


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
    assert model_payload["metrics"]["selected_criterion"] == "BIC"
    assert model_payload["metrics"]["candidate_model_count"] == 3
    assert model_payload["metrics"]["best_model_aic"] == "HKY+G"
    assert model_payload["metrics"]["best_model_aicc"] == "JC"
    assert model_payload["metrics"]["best_model_bic"] == "GTR+G"
    assert model_payload["metrics"]["log_likelihood"] == -123.456
    assert Path(model_payload["data"]["output_paths"]["iqtree_log"]).exists()
    assert Path(model_payload["data"]["output_paths"]["model_candidates"]).exists()

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
    assert ml_payload["metrics"]["selected_model"] == "GTR+G"
    assert ml_payload["metrics"]["log_likelihood"] == -345.678
    assert ml_payload["metrics"]["support_value_count"] == 0
    assert Path(ml_payload["data"]["output_paths"]["iqtree_log"]).exists()
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
    assert fast_payload["metrics"]["approximate_method"] is True
    assert fast_payload["metrics"]["support_label_kind"] == "sh-like-local-support"
    assert fast_payload["metrics"]["support_scale"] == "proportion-0-to-1"
    assert fast_payload["metrics"]["annotated_node_count"] == 2
    assert fast_payload["metrics"]["minimum_local_support"] == 0.62
    assert fast_payload["metrics"]["maximum_local_support"] == 0.98
    assert fast_payload["metrics"]["weakly_supported_clade_count"] == 1
    assert Path(fast_payload["data"]["output_paths"]["support_table"]).exists()
    assert Path(fast_payload["data"]["output_paths"]["low_support_branches"]).exists()
    assert Path(fast_payload["data"]["output_paths"]["support_histogram"]).exists()

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


def test_adapter_model_select_cli_blocks_unsupported_iqtree_before_output_setup(
    tmp_path: Path, capsys
) -> None:
    iqtree = _unsupported_iqtree(tmp_path / "iqtree-unsupported")
    input_path = fixture("alignments/example_alignment.fasta")
    out_dir = tmp_path / "blocked-model-selection"

    exit_code = main(
        [
            "adapter",
            "model-select",
            str(input_path),
            "--out-dir",
            str(out_dir),
            "--executable",
            str(iqtree),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "engine_preflight_workflow_blocked"
    assert payload["errors"][0]["details"]["workflow_id"] == "model-selection"
    assert payload["errors"][0]["details"]["blocking_engines"] == ["IQ-TREE"]
    assert out_dir.exists() is False


def test_adapter_compare_engines_cli_reports_conflicts_and_outputs(
    tmp_path: Path, capsys
) -> None:
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    fasttree = _fake_fasttree(tmp_path / "FastTree-fixture")
    input_path = fixture("alignments/example_alignment.fasta")
    out_dir = tmp_path / "engine-comparison"

    exit_code = main(
        [
            "adapter",
            "compare-engines",
            str(input_path),
            "--out-dir",
            str(out_dir),
            "--prefix",
            "example",
            "--sequence-type",
            "dna",
            "--iqtree-executable",
            str(iqtree),
            "--fasttree-executable",
            str(fasttree),
            "--bootstrap-replicates",
            "1000",
            "--timeout-seconds",
            "30",
            "--incomplete-run-policy",
            "clean",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["selected_model"] == "GTR+G"
    assert payload["metrics"]["shared_clade_count"] == 2
    assert payload["metrics"]["conflicting_clade_count"] == 1
    assert payload["metrics"]["stable_clade_count"] == 1
    assert payload["metrics"]["unstable_clade_count"] == 1
    assert payload["metrics"]["engine_specific_clade_count"] == 0
    assert payload["metrics"]["support_disagreement_count"] == 1
    assert payload["metrics"]["high_support_conflict_count"] == 0
    assert payload["metrics"]["low_support_disagreement_count"] == 0
    assert payload["metrics"]["serious_conflict_count"] == 0
    assert payload["metrics"]["resumed"] is False
    assert payload["metrics"]["timeout_seconds"] == 30.0
    assert Path(payload["data"]["output_paths"]["comparison_report"]).exists()
    assert Path(payload["data"]["output_paths"]["stability_summary"]).exists()
    assert Path(payload["data"]["output_paths"]["conclusion_table"]).exists()
    assert Path(payload["data"]["output_paths"]["support_weighted_conflicts"]).exists()
    assert Path(payload["data"]["output_paths"]["taxon_influence"]).exists()
    assert Path(payload["data"]["output_paths"]["comparison_table"]).exists()
    assert Path(payload["data"]["output_paths"]["shared_clades"]).exists()
    assert Path(payload["data"]["output_paths"]["conflicting_clades"]).exists()
    assert any(
        "normalized to fractions only for side-by-side review" in note
        for note in payload["data"]["notes"]
    )

    resumed_exit = main(
        [
            "adapter",
            "compare-engines",
            str(input_path),
            "--out-dir",
            str(out_dir),
            "--prefix",
            "example",
            "--sequence-type",
            "dna",
            "--iqtree-executable",
            str(iqtree),
            "--fasttree-executable",
            str(fasttree),
            "--bootstrap-replicates",
            "1000",
            "--resume",
            "--timeout-seconds",
            "30",
            "--incomplete-run-policy",
            "clean",
            "--json",
        ]
    )
    resumed_payload = json.loads(capsys.readouterr().out)
    assert resumed_exit == 0
    assert resumed_payload["metrics"]["resumed"] is True
    assert resumed_payload["metrics"]["timeout_seconds"] == 30.0


def test_adapter_infer_large_cli_reports_streamed_size_and_resource_metrics(
    tmp_path: Path, capsys
) -> None:
    fasttree = _fake_fasttree(tmp_path / "FastTree-fixture")
    input_path = fixture("alignments/example_alignment.fasta")
    out_dir = tmp_path / "large-inference"

    exit_code = main(
        [
            "adapter",
            "infer-large",
            str(input_path),
            "--out-dir",
            str(out_dir),
            "--prefix",
            "example",
            "--sequence-type",
            "dna",
            "--executable",
            str(fasttree),
            "--timeout-seconds",
            "30",
            "--incomplete-run-policy",
            "clean",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["sequence_count"] == 4
    assert payload["metrics"]["alignment_length"] == 8
    assert payload["metrics"]["total_site_cells"] == 32
    assert payload["metrics"]["sequence_type"] == "dna"
    assert payload["metrics"]["resumed"] is False
    assert payload["metrics"]["timeout_seconds"] == 30.0
    assert payload["metrics"]["peak_memory_bytes"] is not None
    assert Path(payload["data"]["output_paths"]["tree"]).exists()
    assert Path(payload["data"]["output_paths"]["resource_table"]).exists()
    assert any(
        "scanned linearly before inference" in note for note in payload["data"]["notes"]
    )


def test_adapter_reproducibility_cli_reports_deterministic_outputs(
    tmp_path: Path, capsys
) -> None:
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_alignment.fasta")
    out_dir = tmp_path / "reproducibility"

    exit_code = main(
        [
            "adapter",
            "reproducibility",
            str(input_path),
            "--out-dir",
            str(out_dir),
            "--prefix",
            "example",
            "--sequence-type",
            "dna",
            "--iqtree-executable",
            str(iqtree),
            "--bootstrap-replicates",
            "1000",
            "--repeats",
            "3",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["selected_model"] == "GTR+G"
    assert payload["metrics"]["overall_status"] == "deterministic"
    assert payload["metrics"]["repeat_count"] == 3
    assert payload["metrics"]["unstable_comparison_count"] == 0
    assert payload["metrics"]["equivalent_comparison_count"] == 0
    assert Path(payload["data"]["output_paths"]["runs_table"]).exists()
    assert Path(payload["data"]["output_paths"]["comparison_table"]).exists()
    assert Path(payload["data"]["output_paths"]["support_delta_table"]).exists()
    assert any(
        "classify each rerun as deterministic, equivalent, or unstable" in note
        for note in payload["data"]["notes"]
    )


def test_alignment_partition_summary_cli_writes_summary_table(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "partition-summary.tsv"

    exit_code = main(
        [
            "alignment",
            "partition-summary",
            str(fixture("alignments/example_multilocus_alignment.fasta")),
            str(fixture("alignments/example_multilocus_partitions.txt")),
            "--out",
            str(output_path),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["partition_count"] == 3
    assert payload["metrics"]["unassigned_site_count"] == 0
    assert output_path.exists()
    assert "gene_gamma" in output_path.read_text(encoding="utf-8")


def test_adapter_model_select_cli_supports_partitioned_alignment(
    tmp_path: Path, capsys
) -> None:
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    exit_code = main(
        [
            "adapter",
            "model-select",
            str(fixture("alignments/example_multilocus_alignment.fasta")),
            "--partitions",
            str(fixture("alignments/example_multilocus_partitions.txt")),
            "--out-dir",
            str(tmp_path / "model"),
            "--prefix",
            "partitioned",
            "--executable",
            str(iqtree),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["partitioned"] is True
    assert payload["metrics"]["selected_criterion"] == "BIC"
    assert payload["metrics"]["candidate_model_count"] == 3
    assert payload["metrics"]["best_model_aic"] == "HKY+G"
    assert payload["metrics"]["best_model_aicc"] == "JC"
    assert payload["metrics"]["best_model_bic"] == "GTR+G"
    assert payload["metrics"]["log_likelihood"] == -123.456
    assert Path(payload["data"]["output_paths"]["partition_summary"]).exists()
    assert Path(payload["data"]["output_paths"]["iqtree_log"]).exists()
    assert Path(payload["data"]["output_paths"]["model_candidates"]).exists()


def test_adapter_model_select_cli_supports_protein_alignment(
    tmp_path: Path, capsys
) -> None:
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    exit_code = main(
        [
            "adapter",
            "model-select",
            str(fixture("alignments/example_alignment_protein.fasta")),
            "--out-dir",
            str(tmp_path / "protein-model"),
            "--prefix",
            "protein",
            "--sequence-type",
            "protein",
            "--executable",
            str(iqtree),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["selected_model"] == "LG+G4"
    assert payload["metrics"]["selected_criterion"] == "BIC"
    assert payload["metrics"]["candidate_model_count"] == 3
    assert payload["metrics"]["best_model_bic"] == "LG+G4"
    assert Path(payload["data"]["output_paths"]["model_candidates"]).exists()


def test_adapter_infer_ml_cli_supports_mixed_partition_inputs(
    tmp_path: Path, capsys
) -> None:
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = tmp_path / "mixed-alignment.fasta"
    input_path.write_text(
        ">A\nACGTACMKTW\n>B\nACGTACMKTA\n>C\nACGTACMKTF\n>D\nACGTACMKTY\n",
        encoding="utf-8",
    )
    partition_path = tmp_path / "mixed.partitions"
    partition_path.write_text(
        "DNA,gene_alpha = 1-6\nPROTEIN,gene_beta = 7-10\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "adapter",
            "infer-ml",
            str(input_path),
            "--partitions",
            str(partition_path),
            "--out-dir",
            str(tmp_path / "ml"),
            "--model",
            "MFP",
            "--prefix",
            "mixed",
            "--executable",
            str(iqtree),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["selected_model"] == "GTR+G"
    assert payload["metrics"]["log_likelihood"] == -345.678
    assert payload["metrics"]["support_value_count"] == 0
    assert payload["metrics"]["partitioned"] is True
    assert payload["data"]["run"]["command"].count("-s") == 0
    assert payload["data"]["output_paths"]["partition_scheme"].endswith(".nex")
    assert payload["data"]["output_paths"]["iqtree_log"].endswith(".log")


def test_adapter_bootstrap_cli_reports_support_metrics(tmp_path: Path, capsys) -> None:
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_alignment.fasta")

    exit_code = main(
        [
            "adapter",
            "bootstrap",
            str(input_path),
            "--out-dir",
            str(tmp_path / "bootstrap"),
            "--model",
            "GTR+G",
            "--prefix",
            "example",
            "--replicates",
            "1000",
            "--executable",
            str(iqtree),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["bootstrap_replicates"] == 1000
    assert payload["metrics"]["selected_model"] == "GTR+G"
    assert payload["metrics"]["log_likelihood"] == -234.567
    assert payload["metrics"]["support_value_count"] == 2
    assert payload["metrics"]["minimum_support"] == 88.0
    assert payload["metrics"]["maximum_support"] == 95.0
    assert payload["metrics"]["weakly_supported_clade_count"] == 0
    assert payload["metrics"]["weak_backbone_node_count"] == 0
    assert payload["metrics"]["support_histogram"] == {
        "lt50": 0,
        "50to69": 0,
        "70to89": 1,
        "ge90": 1,
    }
    assert Path(payload["data"]["output_paths"]["iqtree_log"]).exists()
    assert Path(payload["data"]["output_paths"]["bootstrap_trees"]).exists()
    assert Path(payload["data"]["output_paths"]["support_table"]).exists()
    assert Path(payload["data"]["output_paths"]["low_support_branches"]).exists()
    assert Path(payload["data"]["output_paths"]["support_histogram"]).exists()


def test_adapter_bootstrap_cli_surfaces_structured_support_validation_error(
    tmp_path: Path, capsys
) -> None:
    iqtree = _fake_iqtree_without_support_labels(tmp_path / "iqtree-no-support")
    input_path = fixture("alignments/example_alignment.fasta")

    exit_code = main(
        [
            "adapter",
            "bootstrap",
            str(input_path),
            "--out-dir",
            str(tmp_path / "bootstrap"),
            "--model",
            "GTR+G",
            "--prefix",
            "example",
            "--replicates",
            "1000",
            "--executable",
            str(iqtree),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 2
    assert payload["errors"][0]["code"] == "engine_support_values_missing"
    assert payload["errors"][0]["details"]["workflow"] == "bootstrap-support"
    assert payload["errors"][0]["details"]["support_kind"] == "bootstrap support"


def test_adapter_sh_alrt_cli_reports_combined_support_metrics(
    tmp_path: Path, capsys
) -> None:
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_alignment.fasta")

    exit_code = main(
        [
            "adapter",
            "sh-alrt",
            str(input_path),
            "--out-dir",
            str(tmp_path / "sh-alrt"),
            "--model",
            "GTR+G",
            "--prefix",
            "example",
            "--alrt-replicates",
            "1000",
            "--bootstrap-replicates",
            "1000",
            "--executable",
            str(iqtree),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["sh_alrt_replicates"] == 1000
    assert payload["metrics"]["bootstrap_replicates"] == 1000
    assert payload["metrics"]["selected_model"] == "GTR+G"
    assert payload["metrics"]["log_likelihood"] == -222.222
    assert payload["metrics"]["support_value_count"] == 2
    assert payload["metrics"]["sh_alrt_supported_node_count"] == 2
    assert payload["metrics"]["conflicting_support_signal_count"] == 1
    assert payload["metrics"]["minimum_sh_alrt_support"] == 79.0
    assert payload["metrics"]["maximum_sh_alrt_support"] == 82.0
    assert payload["metrics"]["minimum_ufboot_support"] == 96.0
    assert payload["metrics"]["maximum_ufboot_support"] == 97.0
    assert Path(payload["data"]["output_paths"]["support_table"]).exists()
    assert Path(
        payload["data"]["output_paths"]["conflicting_support_branches"]
    ).exists()


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
            "1000",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["selected_model"] == "GTR+G"
    assert payload["metrics"]["sequence_type"] == "dna"
    assert payload["metrics"]["method_tier"] == "supported"
    assert payload["metrics"]["method_inference_mode"] == "inference"
    assert payload["metrics"]["iqtree_seed"] == 1
    assert payload["metrics"]["iqtree_threads"] == 1
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


def test_adapter_fasta_to_tree_cli_blocks_missing_iqtree_before_alignment_starts(
    tmp_path: Path, capsys
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    missing_iqtree = tmp_path / "missing-iqtree"
    input_path = fixture("alignments/example_sequences_raw.fasta")
    out_dir = tmp_path / "blocked-fasta-to-tree"

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
            str(missing_iqtree),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "engine_preflight_workflow_blocked"
    assert payload["errors"][0]["details"]["workflow_id"] == "fasta-to-tree"
    assert payload["errors"][0]["details"]["blocking_engines"] == ["IQ-TREE"]
    assert out_dir.exists() is False


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
            "1000",
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
            "1000",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["trimming_mode"] == "strictplus"
    assert payload["metrics"]["removed_site_count"] == 3
    assert payload["data"]["trimming_workflow"]["run"]["command"][5:] == ["-strictplus"]
    assert (
        payload["data"]["trimming_workflow"]["trimming_summary"]["removed_site_count"]
        == 3
    )
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
        ">Alpha_sample\nACTGACTG\n>rare_taxon\nACTGACTGACTGACTGACTGACTG\n"
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


def test_adapter_fasta_to_tree_cli_rejects_ufboot_replicates_below_iqtree_minimum(
    tmp_path: Path, capsys
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_sequences_raw.fasta")

    exit_code = main(
        [
            "adapter",
            "fasta-to-tree",
            str(input_path),
            "--out-dir",
            str(tmp_path / "too-few-bootstraps"),
            "--prefix",
            "example",
            "--mafft-executable",
            str(mafft),
            "--trimal-executable",
            str(trimal),
            "--iqtree-executable",
            str(iqtree),
            "--bootstrap-replicates",
            "999",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert (
        "ultrafast bootstrap requires at least 1000 replicates"
        in payload["errors"][0]["message"]
    )


def test_adapter_fasta_to_tree_cli_passes_deterministic_iqtree_controls(
    tmp_path: Path, capsys
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_sequences_raw.fasta")

    exit_code = main(
        [
            "adapter",
            "fasta-to-tree",
            str(input_path),
            "--out-dir",
            str(tmp_path / "deterministic-controls"),
            "--prefix",
            "example",
            "--mafft-executable",
            str(mafft),
            "--trimal-executable",
            str(trimal),
            "--iqtree-executable",
            str(iqtree),
            "--iqtree-seed",
            "7",
            "--iqtree-threads",
            "3",
            "--bootstrap-replicates",
            "1000",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["iqtree_seed"] == 7
    assert payload["metrics"]["iqtree_threads"] == 3
    for key in (
        "model_selection_workflow",
        "maximum_likelihood_workflow",
        "bootstrap_workflow",
    ):
        command = payload["data"][key]["run"]["command"]
        seed_index = command.index("-seed")
        thread_index = command.index("-nt")
        assert command[seed_index : seed_index + 2] == ["-seed", "7"]
        assert command[thread_index : thread_index + 2] == ["-nt", "3"]


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
            "--timeout-seconds",
            "30",
            "--incomplete-run-policy",
            "clean",
            "--json",
        ]
    )
    run_payload = json.loads(capsys.readouterr().out)
    assert run_exit == 0
    manifest_path = Path(run_payload["data"]["manifest_path"])
    tree_path = Path(run_payload["data"]["output_paths"]["posterior_trees"])
    trace_path = Path(run_payload["data"]["output_paths"]["parameter_traces"])
    mcmc_path = Path(run_payload["data"]["output_paths"]["mcmc_diagnostics"])
    consensus_path = Path(run_payload["data"]["output_paths"]["consensus_tree"])
    assert run_payload["warnings"] == ["warning: mrbayes fixture posterior run"]
    assert run_payload["metrics"]["resumed"] is False
    assert run_payload["metrics"]["timeout_seconds"] == 30.0

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

    trees_exit = main(["adapter", "mrbayes-trees", str(tree_path), "--json"])
    trees_payload = json.loads(capsys.readouterr().out)
    assert trees_exit == 0
    assert trees_payload["metrics"]["tree_count"] == 4

    mcmc_exit = main(["adapter", "mrbayes-mcmc", str(mcmc_path), "--json"])
    mcmc_payload = json.loads(capsys.readouterr().out)
    assert mcmc_exit == 0
    assert mcmc_payload["metrics"]["comment_count"] == 2

    consensus_exit = main(
        ["adapter", "mrbayes-consensus", str(consensus_path), "--json"]
    )
    consensus_payload = json.loads(capsys.readouterr().out)
    assert consensus_exit == 0
    assert consensus_payload["metrics"]["annotated_node_count"] == 6

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


def test_adapter_mrbayes_prepare_cli_supports_partitioned_alignment(
    tmp_path: Path, capsys
) -> None:
    alignment_path = fixture("alignments/example_multilocus_alignment.fasta")
    partition_path = fixture("alignments/example_multilocus_partitions.txt")
    nexus_path = tmp_path / "partitioned-analysis.nex"

    exit_code = main(
        [
            "adapter",
            "mrbayes-prepare",
            str(alignment_path),
            "--out",
            str(nexus_path),
            "--partitions",
            str(partition_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["partitioned"] is True
    assert payload["metrics"]["partition_count"] == 3
    assert payload["metrics"]["partition_warning_count"] == 0
    assert payload["data"]["partition_names"] == [
        "gene_alpha",
        "gene_beta",
        "gene_gamma",
    ]
    assert payload["data"]["partition_data_types"] == ["DNA"]
    text = nexus_path.read_text(encoding="utf-8")
    assert "partition loci = 3: gene_alpha, gene_beta, gene_gamma;" in text
    assert "set partition=loci;" in text


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
    assert report_payload["metrics"]["warning_count"] >= 1
    assert report_payload["metrics"]["method_tier"] == "parser-only"
    assert report_payload["metrics"]["method_inference_mode"] == "parser-only"
    assert report_payload["warnings"] == [
        "parser-only method tier: this report parses or audits external Bayesian outputs and does not itself infer posterior trees."
    ]
    assert posterior_report_path.exists()


def test_adapter_mrbayes_run_cli_rejects_or_cleans_incomplete_outputs(
    tmp_path: Path, capsys
) -> None:
    malformed = _fake_mrbayes_malformed_outputs(tmp_path / "mb-malformed")
    valid = _fake_mrbayes(tmp_path / "mb-valid")
    alignment_path = fixture("alignments/example_alignment.fasta")
    nexus_path = tmp_path / "analysis.nex"

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

    malformed_exit = main(
        [
            "adapter",
            "mrbayes-run",
            str(nexus_path),
            "--executable",
            str(malformed),
            "--json",
        ]
    )
    malformed_payload = json.loads(capsys.readouterr().out)
    manifest_path = nexus_path.with_suffix("").with_suffix(".manifest.json")
    marker_path = manifest_path.with_suffix(".incomplete.json")

    assert malformed_exit == 2
    assert malformed_payload["status"] == "error"
    assert (
        malformed_payload["errors"][0]["code"]
        == "mrbayes_trace_invalid_parameter_value"
    )
    assert malformed_payload["errors"][0]["details"]["artifact_kind"] == "mrbayes-trace"
    assert (
        malformed_payload["errors"][0]["details"]["failure_reason"]
        == "mrbayes_artifact_malformed"
    )
    assert marker_path.exists()

    rejected_exit = main(
        [
            "adapter",
            "mrbayes-run",
            str(nexus_path),
            "--executable",
            str(valid),
            "--resume",
            "--incomplete-run-policy",
            "reject",
            "--json",
        ]
    )
    rejected_payload = json.loads(capsys.readouterr().out)

    assert rejected_exit == 2
    assert rejected_payload["status"] == "error"
    assert rejected_payload["errors"][0]["code"] == "engine_incomplete_outputs_present"
    assert "incomplete outputs" in rejected_payload["errors"][0]["message"]
    assert (
        rejected_payload["errors"][0]["details"]["failure_reason"]
        == "mrbayes_trace_invalid_parameter_value"
    )
    assert rejected_payload["errors"][0]["details"]["available_actions"] == [
        "resume",
        "clean",
    ]
    assert rejected_payload["errors"][0]["details"]["missing_output_names"] == []
    assert marker_path.exists()

    cleaned_exit = main(
        [
            "adapter",
            "mrbayes-run",
            str(nexus_path),
            "--executable",
            str(valid),
            "--resume",
            "--incomplete-run-policy",
            "clean",
            "--json",
        ]
    )
    cleaned_payload = json.loads(capsys.readouterr().out)

    assert cleaned_exit == 0
    assert cleaned_payload["metrics"]["resumed"] is False
    assert Path(cleaned_payload["data"]["output_paths"]["parameter_traces"]).exists()
    assert marker_path.exists() is False


def test_adapter_mrbayes_run_cli_surfaces_structured_missing_output_error(
    tmp_path: Path, capsys
) -> None:
    executable = _fake_mrbayes_missing_consensus_output(
        tmp_path / "mb-missing-consensus"
    )
    alignment_path = fixture("alignments/example_alignment.fasta")
    nexus_path = tmp_path / "analysis.nex"

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

    exit_code = main(
        [
            "adapter",
            "mrbayes-run",
            str(nexus_path),
            "--executable",
            str(executable),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["errors"][0]["code"] == "engine_required_output_missing"
    assert payload["errors"][0]["details"]["engine_name"] == "MrBayes"
    assert (
        payload["errors"][0]["details"]["failure_reason"]
        == "posterior_artifact_missing"
    )
    assert payload["errors"][0]["details"]["missing_outputs"] == [
        {
            "output_name": "consensus_tree",
            "path": str(Path(f"{nexus_path}.con.tre")),
        }
    ]


def test_adapter_mrbayes_parameters_cli_writes_burnin_aware_summary_table(
    tmp_path: Path, capsys
) -> None:
    trace_path = tmp_path / "diagnostics.run1.p"
    summary_path = tmp_path / "diagnostics-summary.tsv"
    rows = ["Gen\tLnL\tTL"]
    posterior_values = [*range(24), 100.0]
    for index, value in enumerate(posterior_values):
        rows.append(f"{index * 10}\t{value}\t0.5")
    trace_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    exit_code = main(
        [
            "adapter",
            "mrbayes-parameters",
            str(trace_path),
            "--burnin-fraction",
            "0.2",
            "--summary-out",
            str(summary_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    text = summary_path.read_text(encoding="utf-8")

    assert exit_code == 0
    assert payload["metrics"]["kept_row_count"] == 20
    assert payload["metrics"]["parameter_count"] == 2
    assert payload["data"]["first_kept_generation"] == 50
    assert payload["data"]["parameter_summaries"][0]["median"] >= 0
    assert "median\tstandard_deviation" in text
    assert "hpd_95_lower\thpd_95_upper" in text


def test_adapter_mrbayes_burnin_sensitivity_cli_writes_shift_tables(
    tmp_path: Path, capsys
) -> None:
    posterior_path = tmp_path / "burnin-sensitive.run1.t"
    trace_path = tmp_path / "burnin-sensitive.run1.p"
    slice_path = tmp_path / "mrbayes-burnin-slices.tsv"
    parameter_path = tmp_path / "mrbayes-burnin-parameters.tsv"
    clade_path = tmp_path / "mrbayes-burnin-clades.tsv"
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

    exit_code = main(
        [
            "adapter",
            "mrbayes-burnin-sensitivity",
            str(posterior_path),
            "--traces",
            str(trace_path),
            "--slice-out",
            str(slice_path),
            "--parameter-out",
            str(parameter_path),
            "--clade-out",
            str(clade_path),
            "--burnin-fractions",
            "0.0",
            "0.95",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["slice_count"] == 2
    assert payload["metrics"]["unstable_parameter_count"] >= 1
    assert payload["metrics"]["unstable_clade_count"] >= 1
    assert "burnin_fraction\tburnin_tree_count" in slice_path.read_text(
        encoding="utf-8"
    )
    assert "parameter\tminimum_mean" in parameter_path.read_text(encoding="utf-8")
    assert "crosses_majority_threshold" in clade_path.read_text(encoding="utf-8")


def test_adapter_mrbayes_subsample_cli_writes_retained_tree_set_and_metadata(
    tmp_path: Path, capsys
) -> None:
    posterior_path = tmp_path / "posterior.run1.t"
    retained_path = tmp_path / "posterior-retained.nwk"
    table_path = tmp_path / "posterior-retained.tsv"
    posterior_path.write_text(
        "#NEXUS\n"
        "begin trees;\n"
        "tree gen1 = [&R] ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n"
        "tree gen2 = [&R] ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);\n"
        "tree gen3 = [&R] ((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);\n"
        "tree gen4 = [&R] ((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.4);\n"
        "end;\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "adapter",
            "mrbayes-subsample",
            str(posterior_path),
            "--method",
            "random",
            "--burnin-fraction",
            "0.25",
            "--sample-count",
            "2",
            "--seed",
            "11",
            "--tree-set-out",
            str(retained_path),
            "--sample-table-out",
            str(table_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["total_tree_count"] == 4
    assert payload["metrics"]["burnin_tree_count"] == 1
    assert payload["metrics"]["pre_subsampling_tree_count"] == 3
    assert payload["metrics"]["retained_tree_count"] == 2
    assert payload["metrics"]["selection_method"] == "random"
    assert payload["data"]["retained_source_indices"] == [3, 4]
    assert retained_path.read_text(encoding="utf-8").count("\n") == 2
    assert "gen3\t\t3\ttrue" in table_path.read_text(encoding="utf-8")


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
    assert prepare_payload["metrics"]["tip_date_count"] == 4
    assert prepare_payload["metrics"]["warning_count"] == 3
    assert prepare_payload["metrics"]["starting_tree_source"] == "provided-tree"
    assert prepare_payload["metrics"]["beast_data_type"] == "nucleotide"
    assert analysis_path.exists()

    xml_exit = main(
        [
            "adapter",
            "beast-xml",
            str(analysis_path),
            "--json",
        ]
    )
    xml_payload = json.loads(capsys.readouterr().out)
    assert xml_exit == 0
    assert xml_payload["metrics"]["valid"] is True
    assert xml_payload["metrics"]["chain_length"] == 1000000
    assert xml_payload["metrics"]["calibration_count"] == 2
    assert xml_payload["metrics"]["tip_date_count"] == 4

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
        [
            "adapter",
            "beast-log",
            str(fixture("metadata/example_beast.log")),
            "--burnin-fraction",
            "0.25",
            "--json",
        ]
    )
    log_payload = json.loads(capsys.readouterr().out)
    assert log_exit == 0
    assert log_payload["metrics"]["row_count"] == 4
    assert log_payload["metrics"]["kept_row_count"] == 3

    convergence_exit = main(
        [
            "adapter",
            "beast-convergence",
            str(fixture("metadata/example_beast.log")),
            "--burnin-fraction",
            "0.25",
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
    assert convergence_payload["metrics"]["sample_count"] == 3

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
    assert report_payload["metrics"]["warning_count"] >= 1
    assert report_payload["metrics"]["method_tier"] == "parser-only"
    assert report_payload["warnings"] == [
        "parser-only method tier: this report parses or audits external Bayesian outputs and does not itself infer posterior trees."
    ]
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


def test_adapter_beast_run_cli_supports_resume_timeout_and_overwrite_controls(
    tmp_path: Path, capsys
) -> None:
    executable = _fake_beast(tmp_path / "beast-fixture")
    analysis_path = tmp_path / "analysis.xml"
    analysis_path.write_text(
        '<beast version="2.7" namespace="beast.base.evolution.alignment:beast.base.evolution.tree"/>\n',
        encoding="utf-8",
    )

    first_exit = main(
        [
            "adapter",
            "beast-run",
            str(analysis_path),
            "--executable",
            str(executable),
            "--threads",
            "2",
            "--seed",
            "7",
            "--timeout-seconds",
            "30",
            "--incomplete-run-policy",
            "clean",
            "--json",
        ]
    )
    first_payload = json.loads(capsys.readouterr().out)
    assert first_exit == 0
    assert first_payload["metrics"]["threads"] == 2
    assert first_payload["metrics"]["seed"] == 7
    assert first_payload["metrics"]["overwrite"] is True
    assert first_payload["metrics"]["resumed"] is False
    assert first_payload["metrics"]["timeout_seconds"] == 30.0
    assert Path(first_payload["data"]["output_paths"]["posterior_log"]).exists()
    assert Path(first_payload["data"]["output_paths"]["posterior_trees"]).exists()

    second_exit = main(
        [
            "adapter",
            "beast-run",
            str(analysis_path),
            "--executable",
            str(executable),
            "--threads",
            "2",
            "--seed",
            "7",
            "--resume",
            "--no-overwrite",
            "--timeout-seconds",
            "30",
            "--incomplete-run-policy",
            "clean",
            "--json",
        ]
    )
    second_payload = json.loads(capsys.readouterr().out)
    assert second_exit == 0
    assert second_payload["metrics"]["overwrite"] is False
    assert second_payload["metrics"]["resumed"] is False

    resumed_exit = main(
        [
            "adapter",
            "beast-run",
            str(analysis_path),
            "--executable",
            str(executable),
            "--threads",
            "2",
            "--seed",
            "7",
            "--resume",
            "--no-overwrite",
            "--timeout-seconds",
            "30",
            "--incomplete-run-policy",
            "clean",
            "--json",
        ]
    )
    resumed_payload = json.loads(capsys.readouterr().out)
    assert resumed_exit == 0
    assert resumed_payload["metrics"]["overwrite"] is False
    assert resumed_payload["metrics"]["resumed"] is True
    assert resumed_payload["metrics"]["timeout_seconds"] == 30.0


def test_adapter_beast_run_cli_rejects_or_cleans_incomplete_outputs(
    tmp_path: Path, capsys
) -> None:
    malformed = _fake_beast_malformed_outputs(tmp_path / "beast-malformed")
    valid = _fake_beast(tmp_path / "beast-valid")
    analysis_path = tmp_path / "analysis.xml"
    analysis_path.write_text(
        '<beast version="2.7" namespace="beast.base.evolution.alignment:beast.base.evolution.tree"/>\n',
        encoding="utf-8",
    )

    malformed_exit = main(
        [
            "adapter",
            "beast-run",
            str(analysis_path),
            "--executable",
            str(malformed),
            "--seed",
            "7",
            "--json",
        ]
    )
    malformed_payload = json.loads(capsys.readouterr().out)
    manifest_path = analysis_path.with_suffix(".manifest.json")
    marker_path = manifest_path.with_suffix(".incomplete.json")

    assert malformed_exit == 2
    assert malformed_payload["status"] == "error"
    assert malformed_payload["errors"][0]["code"] == "beast_log_invalid_parameter_value"
    assert malformed_payload["errors"][0]["details"]["artifact_kind"] == "beast-log"
    assert (
        malformed_payload["errors"][0]["details"]["failure_reason"]
        == "beast_artifact_malformed"
    )
    assert marker_path.exists()

    rejected_exit = main(
        [
            "adapter",
            "beast-run",
            str(analysis_path),
            "--executable",
            str(valid),
            "--seed",
            "7",
            "--resume",
            "--incomplete-run-policy",
            "reject",
            "--json",
        ]
    )
    rejected_payload = json.loads(capsys.readouterr().out)

    assert rejected_exit == 2
    assert rejected_payload["status"] == "error"
    assert rejected_payload["errors"][0]["code"] == "engine_incomplete_outputs_present"
    assert "incomplete outputs" in rejected_payload["errors"][0]["message"]
    assert (
        rejected_payload["errors"][0]["details"]["failure_reason"]
        == "beast_log_invalid_parameter_value"
    )
    assert rejected_payload["errors"][0]["details"]["available_actions"] == [
        "resume",
        "clean",
    ]
    assert rejected_payload["errors"][0]["details"]["missing_output_names"] == []
    assert marker_path.exists()

    cleaned_exit = main(
        [
            "adapter",
            "beast-run",
            str(analysis_path),
            "--executable",
            str(valid),
            "--seed",
            "7",
            "--resume",
            "--incomplete-run-policy",
            "clean",
            "--json",
        ]
    )
    cleaned_payload = json.loads(capsys.readouterr().out)

    assert cleaned_exit == 0
    assert cleaned_payload["metrics"]["resumed"] is False
    assert Path(cleaned_payload["data"]["output_paths"]["posterior_log"]).exists()
    assert marker_path.exists() is False


def test_adapter_beast_run_cli_surfaces_structured_missing_output_error(
    tmp_path: Path, capsys
) -> None:
    executable = _fake_beast_missing_tree_output(tmp_path / "beast-missing-tree")
    analysis_path = tmp_path / "analysis.xml"
    analysis_path.write_text(
        '<beast version="2.7" namespace="beast.base.evolution.alignment:beast.base.evolution.tree"/>\n',
        encoding="utf-8",
    )

    exit_code = main(
        [
            "adapter",
            "beast-run",
            str(analysis_path),
            "--executable",
            str(executable),
            "--seed",
            "7",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["errors"][0]["code"] == "engine_required_output_missing"
    assert payload["errors"][0]["details"]["engine_name"] == "BEAST"
    assert (
        payload["errors"][0]["details"]["failure_reason"]
        == "posterior_artifact_missing"
    )
    assert payload["errors"][0]["details"]["missing_outputs"] == [
        {
            "output_name": "posterior_trees",
            "path": str(analysis_path.with_name("analysis.7.trees")),
        }
    ]


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
    assert conflicts_payload["metrics"]["conflict_count"] == 4

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
    assert package_payload["metrics"]["artifact_count"] == 15
    assert package_payload["metrics"]["publication_ready"] is True
    assert package_payload["metrics"]["support_labels_validated"] is True
    assert (package_dir / "consensus-tree.svg").exists()
    assert (package_dir / "clade-support-plot.svg").exists()
    assert (package_dir / "unstable-taxa-plot.svg").exists()
    assert (package_dir / "topology-clusters-plot.svg").exists()
    assert (package_dir / "figure-legend.tsv").exists()
    assert (package_dir / "figure-caption.md").exists()
    assert (package_dir / "uncertainty-review.html").exists()


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
            "--analysis-xml",
            str(fixture("metadata/beast2_strict_yule_posterior.xml")),
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
    methods_text = methods_path.read_text(encoding="utf-8")
    assert "strict" in methods_text
    assert "beast2_strict_yule_posterior.xml" in methods_text


def test_adapter_beast_log_cli_writes_summary_table_with_parameter_categories(
    tmp_path: Path, capsys
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

    exit_code = main(
        [
            "adapter",
            "beast-log",
            str(log_path),
            "--burnin-fraction",
            "0.25",
            "--summary-out",
            str(summary_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    text = summary_path.read_text(encoding="utf-8")

    assert exit_code == 0
    assert payload["metrics"]["prior_parameter_count"] == 1
    assert payload["metrics"]["clock_parameter_count"] == 1
    assert payload["metrics"]["tree_parameter_count"] == 2
    assert payload["data"]["summary"]["kept_row_count"] == 3
    assert "median\tstandard_deviation" in text
    assert "prior\tprior\t3" in text
    assert "tree\tbirthRate\t3" in text


def test_adapter_beast_log_cli_returns_structured_error_payload_for_malformed_log(
    tmp_path: Path, capsys
) -> None:
    log_path = tmp_path / "malformed-beast.log"
    log_path.write_text(
        "# malformed BEAST log\n"
        "state\tposterior\tlikelihood\n"
        "0\t-100.0\t-90.0\n"
        "10\tbad\t-89.0\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "adapter",
            "beast-log",
            str(log_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "beast_log_invalid_parameter_value"
    assert payload["errors"][0]["details"]["artifact_kind"] == "beast-log"
    assert payload["errors"][0]["details"]["column"] == "posterior"


def test_adapter_beast_parameters_cli_reports_posterior_diagnostics(
    tmp_path: Path, capsys
) -> None:
    log_path = tmp_path / "posterior-diagnostics.log"
    summary_path = tmp_path / "posterior-diagnostics.tsv"
    lines = [
        "# posterior diagnostics fixture",
        "state\tposterior\tclockRate",
    ]
    for state, posterior in enumerate([*range(19), 100.0]):
        lines.append(f"{state}\t{posterior}\t0.001")
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    exit_code = main(
        [
            "adapter",
            "beast-parameters",
            str(log_path),
            "--summary-out",
            str(summary_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    text = summary_path.read_text(encoding="utf-8")

    assert exit_code == 0
    assert payload["metrics"]["parameter_count"] == 2
    assert payload["metrics"]["posterior_parameter_count"] == 1
    assert payload["data"]["parameter_summaries"][0]["hpd_95_upper"] >= 0
    assert "median\tstandard_deviation" in text
    assert "hpd_95_lower\thpd_95_upper" in text


def test_adapter_beast_burnin_sensitivity_cli_writes_shift_tables(
    tmp_path: Path, capsys
) -> None:
    posterior_path = tmp_path / "burnin-sensitive.trees"
    log_path = tmp_path / "burnin-sensitive.log"
    slice_path = tmp_path / "beast-burnin-slices.tsv"
    parameter_path = tmp_path / "beast-burnin-parameters.tsv"
    clade_path = tmp_path / "beast-burnin-clades.tsv"
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

    exit_code = main(
        [
            "adapter",
            "beast-burnin-sensitivity",
            str(posterior_path),
            "--log",
            str(log_path),
            "--slice-out",
            str(slice_path),
            "--parameter-out",
            str(parameter_path),
            "--clade-out",
            str(clade_path),
            "--burnin-fractions",
            "0.0",
            "0.95",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["slice_count"] == 2
    assert payload["metrics"]["unstable_parameter_count"] >= 1
    assert payload["metrics"]["unstable_clade_count"] >= 1
    assert "consensus_newick" in slice_path.read_text(encoding="utf-8")
    assert "common_hpd_95_lower" in parameter_path.read_text(encoding="utf-8")
    assert "crosses_majority_threshold" in clade_path.read_text(encoding="utf-8")


def test_adapter_beast_trees_cli_writes_normalized_tree_set(
    tmp_path: Path, capsys
) -> None:
    posterior_path = tmp_path / "posterior.trees"
    normalized_path = tmp_path / "posterior.nwk"
    posterior_path.write_text(
        "#NEXUS\n"
        "Begin trees;\n"
        "  Translate\n"
        "    1 A,\n"
        "    2 B,\n"
        "    3 C,\n"
        "    4 D\n"
        "  ;\n"
        "tree STATE_0 = ((1:0.1,2:0.1):0.2,(3:0.1,4:0.1):0.2):0.0;\n"
        "tree STATE_10 = ((1:0.1,3:0.1):0.2,(2:0.1,4:0.1):0.2):0.0;\n"
        "tree STATE_20 = ((1:0.1,2:0.1):0.2,(3:0.1,4:0.1):0.2):0.0;\n"
        "tree STATE_30 = ((1:0.1,2:0.1):0.2,(3:0.1,4:0.1):0.2):0.0;\n"
        "End;\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "adapter",
            "beast-trees",
            str(posterior_path),
            "--burnin-fraction",
            "0.25",
            "--tree-set-out",
            str(normalized_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    text = normalized_path.read_text(encoding="utf-8")

    assert exit_code == 0
    assert payload["metrics"]["total_tree_count"] == 4
    assert payload["metrics"]["kept_tree_count"] == 3
    assert payload["metrics"]["clade_count"] >= 1
    assert payload["data"]["sampled_states"] == [10, 20, 30]
    assert "((A:0.1,B:0.1):0.2,(C:0.1,D:0.2)" not in text
    assert "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);" in text


def test_adapter_beast_subsample_cli_writes_retained_tree_set_and_metadata(
    tmp_path: Path, capsys
) -> None:
    posterior_path = tmp_path / "posterior.trees"
    retained_path = tmp_path / "posterior-retained.nwk"
    table_path = tmp_path / "posterior-retained.tsv"
    posterior_path.write_text(
        "#NEXUS\n"
        "Begin trees;\n"
        "tree STATE_0 = ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2):0.0;\n"
        "tree STATE_10 = ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2):0.0;\n"
        "tree STATE_20 = ((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3):0.0;\n"
        "tree STATE_30 = ((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.4):0.0;\n"
        "End;\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "adapter",
            "beast-subsample",
            str(posterior_path),
            "--method",
            "evenly-spaced",
            "--burnin-fraction",
            "0.25",
            "--thinning-interval",
            "2",
            "--tree-set-out",
            str(retained_path),
            "--sample-table-out",
            str(table_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["total_tree_count"] == 4
    assert payload["metrics"]["burnin_tree_count"] == 1
    assert payload["metrics"]["pre_subsampling_tree_count"] == 3
    assert payload["metrics"]["retained_tree_count"] == 2
    assert payload["metrics"]["selection_method"] == "evenly-spaced"
    assert payload["data"]["retained_source_indices"] == [2, 4]
    assert retained_path.read_text(encoding="utf-8").count("\n") == 2
    assert "STATE_10\t10\t\ttrue" in table_path.read_text(encoding="utf-8")


def test_adapter_beast_consensus_cli_writes_consensus_and_clade_ledger(
    tmp_path: Path, capsys
) -> None:
    posterior_path = tmp_path / "posterior.trees"
    consensus_path = tmp_path / "posterior-consensus.nwk"
    retained_path = tmp_path / "posterior-retained.nwk"
    clade_table_path = tmp_path / "posterior-clades.tsv"
    posterior_path.write_text(
        "#NEXUS\n"
        "Begin trees;\n"
        "tree STATE_0 = ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2):0.0;\n"
        "tree STATE_10 = ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2):0.0;\n"
        "tree STATE_20 = ((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.4):0.0;\n"
        "tree STATE_30 = ((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3):0.0;\n"
        "End;\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "adapter",
            "beast-consensus",
            str(posterior_path),
            "--burnin-fraction",
            "0.25",
            "--out",
            str(consensus_path),
            "--tree-set-out",
            str(retained_path),
            "--clade-table-out",
            str(clade_table_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["kept_tree_count"] == 3
    assert payload["metrics"]["annotated_node_count"] == 2
    assert payload["metrics"]["clade_frequency_count"] == 4
    assert payload["data"]["maximum_posterior_probability"] == pytest.approx(2 / 3)
    assert consensus_path.read_text(encoding="utf-8").strip() == (
        "((A:0.1,B:0.1)0.666666666666667:0.35,(C:0.1,D:0.1)0.666666666666667:0.35);"
    )
    assert retained_path.read_text(encoding="utf-8").count("\n") == 3
    assert "A|B\t2\t0.666666666666667" in clade_table_path.read_text(encoding="utf-8")


def test_adapter_beast_diversity_cli_writes_distance_and_instability_ledgers(
    tmp_path: Path, capsys
) -> None:
    posterior_path = tmp_path / "posterior.trees"
    retained_path = tmp_path / "posterior-retained.nwk"
    distance_path = tmp_path / "posterior-distances.tsv"
    topology_path = tmp_path / "posterior-topologies.tsv"
    unstable_path = tmp_path / "posterior-unstable-clades.tsv"
    posterior_path.write_text(
        "#NEXUS\n"
        "Begin trees;\n"
        "tree STATE_0 = ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2):0.0;\n"
        "tree STATE_10 = ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2):0.0;\n"
        "tree STATE_20 = ((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.4):0.0;\n"
        "tree STATE_30 = ((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3):0.0;\n"
        "End;\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "adapter",
            "beast-diversity",
            str(posterior_path),
            "--burnin-fraction",
            "0.25",
            "--tree-set-out",
            str(retained_path),
            "--distance-out",
            str(distance_path),
            "--topology-out",
            str(topology_path),
            "--unstable-clade-out",
            str(unstable_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["kept_tree_count"] == 3
    assert payload["metrics"]["rooted_topology_count"] == 2
    assert payload["metrics"]["pair_count"] == 3
    assert payload["metrics"]["unstable_clade_count"] >= 1
    assert payload["data"]["dominant_topology_frequency"] == pytest.approx(2 / 3)
    assert retained_path.read_text(encoding="utf-8").count("\n") == 3
    assert "normalized_robinson_foulds" in distance_path.read_text(encoding="utf-8")
    assert "representative_newick" in topology_path.read_text(encoding="utf-8")
    assert "support_classification" in unstable_path.read_text(encoding="utf-8")
