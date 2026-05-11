from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.engines.fasta_to_tree import (
    FastaToTreeModelRow,
    FastaToTreeSupportRow,
    infer_unaligned_sequence_type,
    run_fasta_to_tree_workflow,
    write_fasta_to_tree_model_table,
    write_fasta_to_tree_support_table,
)

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
is_protein = "-st" in args and args[args.index("-st") + 1] == "AA"
selected_model = "LG+G" if is_protein else "GTR+G"
if "-m" in args and args[args.index("-m") + 1] == "MF":
    prefix.with_suffix(".iqtree").write_text(
        f"Best-fit model according to BIC: {selected_model}\\nWARNING: model search used a fixture backend\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".model").write_text(
        f"Best-fit model: {selected_model}\\n",
        encoding="utf-8",
    )
    print("warning: iqtree fixture model selection", file=sys.stderr)
    raise SystemExit(0)

if "-bb" in args:
    prefix.with_suffix(".treefile").write_text(
        "((A:0.1,B:0.1)95:0.2,(C:0.1,D:0.1)88:0.2);\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".ufboot").write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".iqtree").write_text(
        "Bootstrap analysis completed\\n",
        encoding="utf-8",
    )
    print("warning: iqtree fixture bootstrap", file=sys.stderr)
    raise SystemExit(0)

prefix.with_suffix(".treefile").write_text(
    "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
    encoding="utf-8",
)
prefix.with_suffix(".iqtree").write_text(
    "Tree inference completed\\n",
    encoding="utf-8",
)
print("warning: iqtree fixture tree inference", file=sys.stderr)
""",
    )


def test_infer_unaligned_sequence_type_detects_dna_and_protein_inputs() -> None:
    assert (
        infer_unaligned_sequence_type(
            [("A", "ACTGACTG"), ("B", "ACTGACTA"), ("C", "ACTGACGG")]
        )
        == "dna"
    )
    assert (
        infer_unaligned_sequence_type(
            [("P1", "MKTWFLIM"), ("P2", "MKTWYLIM"), ("P3", "MRTWFLVM")]
        )
        == "protein"
    )


def test_write_fasta_to_tree_tables_emits_expected_tsv(tmp_path: Path) -> None:
    model_path = tmp_path / "example.model.tsv"
    support_path = tmp_path / "example.support.tsv"

    write_fasta_to_tree_model_table(
        model_path,
        [
            FastaToTreeModelRow(
                workflow="model-selection",
                engine_name="iqtree",
                sequence_type="dna",
                selected_model="GTR+G",
                report_selected_model="GTR+G",
                artifact_selected_model="GTR+G",
                model_consistent=True,
                alignment_path=Path("example.aln"),
                trimmed_alignment_path=Path("example.trimmed.aln"),
                manifest_path=Path("engine-artifacts/example/model.manifest.json"),
            )
        ],
    )
    write_fasta_to_tree_support_table(
        support_path,
        [
            FastaToTreeSupportRow(
                node="A|B",
                descendant_taxa=("A", "B"),
                support=95.0,
                support_fraction=0.95,
                is_backbone=True,
            )
        ],
    )

    assert model_path.read_text(encoding="utf-8") == (
        "workflow\tengine_name\tsequence_type\tselected_model\t"
        "report_selected_model\tartifact_selected_model\tmodel_consistent\t"
        "alignment_path\ttrimmed_alignment_path\tmanifest_path\n"
        "model-selection\tiqtree\tdna\tGTR+G\tGTR+G\tGTR+G\ttrue\t"
        "example.aln\texample.trimmed.aln\t"
        "engine-artifacts/example/model.manifest.json\n"
    )
    assert support_path.read_text(encoding="utf-8") == (
        "node\tdescendant_taxa\tsupport\tsupport_fraction\tis_backbone\n"
        "A|B\tA,B\t95\t0.95\ttrue\n"
    )


def test_run_fasta_to_tree_workflow_materializes_expected_outputs_for_three_datasets(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    cases = [
        (
            fixture("alignments/example_sequences_raw.fasta"),
            "raw-dna",
            "dna",
            "GTR+G",
        ),
        (
            fixture("alignments/example_alignment.fasta"),
            "aligned-dna",
            "dna",
            "GTR+G",
        ),
        (
            fixture("alignments/example_alignment_protein.fasta"),
            "protein",
            "protein",
            "LG+G",
        ),
    ]

    for input_path, prefix, expected_sequence_type, expected_model in cases:
        report = run_fasta_to_tree_workflow(
            input_path,
            out_dir=tmp_path / prefix,
            prefix=prefix,
            mafft_executable=mafft,
            trimal_executable=trimal,
            iqtree_executable=iqtree,
            bootstrap_replicates=200,
        )

        assert report.sequence_type == expected_sequence_type
        assert report.selected_model == expected_model
        assert report.output_paths["alignment"].suffix == ".aln"
        assert report.output_paths["trimmed_alignment"].name.endswith(".trimmed.aln")
        assert report.output_paths["tree"].suffix == ".tree"
        assert report.output_paths["log"].suffix == ".log"
        assert report.output_paths["model_table"].name.endswith(".model.tsv")
        assert report.output_paths["support_table"].name.endswith(".support.tsv")
        for path in report.output_paths.values():
            assert path.exists()

        model_rows = list(
            csv.DictReader(
                report.output_paths["model_table"].open(encoding="utf-8"),
                delimiter="\t",
            )
        )
        support_rows = list(
            csv.DictReader(
                report.output_paths["support_table"].open(encoding="utf-8"),
                delimiter="\t",
            )
        )
        assert model_rows[0]["selected_model"] == expected_model
        assert model_rows[0]["model_consistent"] == "true"
        assert support_rows[0]["support"] == "95"
        log_text = report.output_paths["log"].read_text(encoding="utf-8")
        assert "selected_model:" in log_text
        assert "warning: iqtree fixture bootstrap" in log_text
