from __future__ import annotations
from pathlib import Path

from bijux_phylogenetics.engines import (
    compare_fast_and_ml_trees,
    run_alignment_trimming,
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_fast_tree_inference,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
    run_multiple_sequence_alignment,
)
from bijux_phylogenetics.io.fasta import load_fasta_alignment


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
        handle.write(f">{identifier}\\n{sequence[:-1]}\\n")
print("warning: trimal fixture trimmed one trailing site", file=sys.stderr)
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

if "-con" in args:
    prefix.with_suffix(".contree").write_text("((A:0.1,B:0.1)90:0.2,(C:0.1,D:0.1)85:0.2);\\n", encoding="utf-8")
    print("warning: iqtree fixture consensus", file=sys.stderr)
    raise SystemExit(0)

if "-bb" in args:
    prefix.with_suffix(".treefile").write_text("((A:0.1,B:0.1)95:0.2,(C:0.1,D:0.1)88:0.2);\\n", encoding="utf-8")
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


def test_run_multiple_sequence_alignment_captures_logs_version_and_manifest(tmp_path: Path) -> None:
    executable = _fake_mafft(tmp_path / "mafft-fixture")
    input_path = tmp_path / "unaligned.fasta"
    input_path.write_text(">A\nACTG\n>B\nACTGA\n>C\nACT\n", encoding="utf-8")
    output_path = tmp_path / "aligned.fasta"

    report = run_multiple_sequence_alignment(input_path, output_path, executable=executable)

    records = load_fasta_alignment(output_path)
    assert [len(record.sequence) for record in records] == [5, 5, 5]
    assert "mafft v7.999" in report.run.version.text
    assert report.run.command[0] == str(executable)
    assert report.run.warning_lines == ["WARNING: mafft fixture inserted alignment padding"]
    assert report.manifest_path.exists()


def test_run_alignment_trimming_writes_trimmed_alignment_and_warning_manifest(tmp_path: Path) -> None:
    executable = _fake_trimal(tmp_path / "trimal-fixture")
    input_path = fixture("alignments/example_alignment_trim.fasta")
    output_path = tmp_path / "trimmed.fasta"

    report = run_alignment_trimming(input_path, output_path, executable=executable, gap_threshold=0.2)

    records = load_fasta_alignment(output_path)
    assert len(records[0].sequence) == len(load_fasta_alignment(input_path)[0].sequence) - 1
    assert report.run.warning_lines == ["warning: trimal fixture trimmed one trailing site"]
    assert report.manifest_path.exists()


def test_run_model_selection_parses_best_fit_model_and_writes_manifest(tmp_path: Path) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_alignment.fasta")

    report = run_model_selection(input_path, out_dir=tmp_path / "model", executable=executable, prefix="example")

    assert report.selected_model == "GTR+G"
    selected_model_path = report.output_paths["selected_model"]
    assert selected_model_path.read_text(encoding="utf-8").strip() == "GTR+G"
    assert report.run.warning_lines == ["warning: iqtree fixture model selection"]
    assert report.manifest_path.exists()


def test_run_ml_bootstrap_consensus_and_fast_tree_workflows(tmp_path: Path) -> None:
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    fasttree = _fake_fasttree(tmp_path / "FastTree-fixture")
    input_path = fixture("alignments/example_alignment.fasta")

    ml_report = run_maximum_likelihood_tree_inference(
        input_path,
        out_dir=tmp_path / "ml",
        model="GTR+G",
        executable=iqtree,
        prefix="example",
    )
    bootstrap_report = run_bootstrap_support_estimation(
        input_path,
        out_dir=tmp_path / "bootstrap",
        model="GTR+G",
        executable=iqtree,
        prefix="example",
        replicates=1000,
    )
    consensus_report = run_bootstrap_consensus_tree(
        bootstrap_report.output_paths["bootstrap_trees"],
        out_dir=tmp_path / "consensus",
        executable=iqtree,
        prefix="example",
    )
    fast_report = run_fast_tree_inference(input_path, tmp_path / "fasttree.nwk", executable=fasttree)

    assert ml_report.output_paths["tree"].exists()
    assert bootstrap_report.output_paths["bootstrap_trees"].exists()
    assert consensus_report.output_paths["consensus_tree"].exists()
    assert fast_report.output_paths["tree"].exists()
    assert ml_report.run.warning_lines == ["warning: iqtree fixture tree inference"]
    assert bootstrap_report.run.warning_lines == ["warning: iqtree fixture bootstrap"]
    assert consensus_report.run.warning_lines == ["warning: iqtree fixture consensus"]
    assert fast_report.run.warning_lines == ["warning: fasttree fixture approximate support only"]


def test_compare_fast_and_ml_trees_builds_html_report(tmp_path: Path) -> None:
    fast_tree_path = tmp_path / "fast.nwk"
    ml_tree_path = tmp_path / "ml.nwk"
    fast_tree_path.write_text("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);\n", encoding="utf-8")
    ml_tree_path.write_text("((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n", encoding="utf-8")

    comparison = compare_fast_and_ml_trees(fast_tree_path, ml_tree_path, out_path=tmp_path / "comparison.html")

    assert comparison.comparison_report.output_path.exists()
    assert comparison.comparison_report.topology.shared_taxa == ["A", "B", "C", "D"]
