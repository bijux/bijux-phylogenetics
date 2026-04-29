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


def test_adapter_inspect_cli_reports_engine_version(tmp_path: Path, capsys) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")

    exit_code = main(["adapter", "inspect", "iqtree", "--executable", str(executable), "--json"])

    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["text"] == "IQ-TREE multicore version 2.9.9"


def test_adapter_align_cli_writes_aligned_fasta_and_manifest(tmp_path: Path, capsys) -> None:
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


def test_adapter_model_select_and_compare_cli_produce_outputs(tmp_path: Path, capsys) -> None:
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
    assert fast_payload["warnings"] == ["warning: fasttree fixture approximate support only"]

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
