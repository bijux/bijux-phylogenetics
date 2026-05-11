from __future__ import annotations

import csv
from pathlib import Path

import pytest

from bijux_phylogenetics.engines.fasta_to_tree import (
    FastaToTreeModelRow,
    FastaToTreeSupportRow,
    infer_unaligned_sequence_type,
    run_fasta_to_tree_workflow,
    write_fasta_to_tree_model_table,
    write_fasta_to_tree_support_table,
)
from bijux_phylogenetics.errors import InvalidAlignmentError

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


def test_infer_unaligned_sequence_type_returns_unknown_for_mixed_inputs() -> None:
    assert (
        infer_unaligned_sequence_type(
            [("dna_like", "ACTGACTG"), ("rna_like", "ACUGACUG")]
        )
        == "unknown"
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


def test_run_fasta_to_tree_workflow_rejects_invalid_raw_input_without_repair(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    with pytest.raises(InvalidAlignmentError, match="duplicate identifiers"):
        run_fasta_to_tree_workflow(
            fixture("alignments/example_sequences_invalid_input.fasta"),
            out_dir=tmp_path / "strict",
            prefix="strict",
            sequence_type="dna",
            mafft_executable=mafft,
            trimal_executable=trimal,
            iqtree_executable=iqtree,
        )


def test_run_fasta_to_tree_workflow_rejects_mixed_raw_input_without_explicit_type(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    mixed_input = tmp_path / "mixed.fasta"
    mixed_input.write_text(
        ">dna_like\nACTGACTG\n>rna_like\nACUGACUG\n",
        encoding="utf-8",
    )

    with pytest.raises(
        InvalidAlignmentError, match="thymine-bearing and uracil-bearing"
    ):
        run_fasta_to_tree_workflow(
            mixed_input,
            out_dir=tmp_path / "mixed",
            prefix="mixed",
            mafft_executable=mafft,
            trimal_executable=trimal,
            iqtree_executable=iqtree,
        )


def test_run_fasta_to_tree_workflow_repairs_invalid_raw_input_when_requested(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    report = run_fasta_to_tree_workflow(
        fixture("alignments/example_sequences_invalid_input.fasta"),
        out_dir=tmp_path / "repaired",
        prefix="repaired",
        sequence_type="dna",
        mafft_executable=mafft,
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        normalize_identifiers=True,
        remove_invalid_records=True,
    )

    assert report.prepared_input_path.name == "input-curation.fasta"
    assert report.input_repair is not None
    assert report.repaired_input_validation is not None
    assert report.input_repair.output_path == report.prepared_input_path
    assert report.prepared_input_path.read_text(encoding="utf-8") == (
        ">Alpha_sample\nACTGACTG\n"
        ">rare_taxon\nACTGACTGACTGACTGACTGACTG\n"
    )
    assert any("duplicate sequence identifiers" in warning for warning in report.warnings)
    assert "prepared_input_path:" in report.output_paths["log"].read_text(
        encoding="utf-8"
    )


def test_run_fasta_to_tree_workflow_can_force_declared_type_on_mixed_input(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    mixed_input = tmp_path / "mixed.fasta"
    mixed_input.write_text(
        ">dna_a\nACTGACTG\n>dna_b\nACTGACTA\n>rna_like\nACUGACUG\n",
        encoding="utf-8",
    )

    report = run_fasta_to_tree_workflow(
        mixed_input,
        out_dir=tmp_path / "forced-dna",
        prefix="forced-dna",
        sequence_type="dna",
        mafft_executable=mafft,
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        normalize_identifiers=False,
        remove_invalid_records=True,
    )

    assert report.sequence_type == "dna"
    assert report.input_repair is not None
    assert [row.identifier for row in report.input_repair.removed_records] == [
        "rna_like"
    ]
    assert report.prepared_input_path.read_text(encoding="utf-8") == (
        ">dna_a\nACTGACTG\n>dna_b\nACTGACTA\n"
    )
    assert "raw sequence type detection: dna (medium)" in report.notes


def test_run_fasta_to_tree_workflow_passes_named_mafft_mode_to_alignment_step(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    report = run_fasta_to_tree_workflow(
        fixture("alignments/example_sequences_raw.fasta"),
        out_dir=tmp_path / "workflow-linsi",
        prefix="workflow-linsi",
        mafft_executable=mafft,
        alignment_mode="linsi",
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        bootstrap_replicates=200,
    )

    assert report.alignment_workflow.run.command[1:-1] == [
        "--localpair",
        "--maxiterate",
        "1000",
    ]
    assert report.alignment_workflow.notes[0] == "mafft alignment mode: linsi"
    assert "mafft alignment mode: linsi" in report.notes
    assert "command:" in report.output_paths["log"].read_text(encoding="utf-8")


def test_run_fasta_to_tree_workflow_passes_named_trimal_mode_to_trimming_step(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    report = run_fasta_to_tree_workflow(
        fixture("alignments/example_sequences_raw.fasta"),
        out_dir=tmp_path / "workflow-strictplus",
        prefix="workflow-strictplus",
        mafft_executable=mafft,
        trimal_executable=trimal,
        trimming_mode="strictplus",
        iqtree_executable=iqtree,
        bootstrap_replicates=200,
    )

    assert report.trimming_workflow.run.command[5:] == ["-strictplus"]
    assert report.trimming_workflow.trimming_summary is not None
    assert report.trimming_workflow.trimming_summary.mode == "strictplus"
    assert report.trimming_workflow.trimming_summary.removed_site_count == 3
    assert "trimal trimming mode: strictplus" in report.notes
