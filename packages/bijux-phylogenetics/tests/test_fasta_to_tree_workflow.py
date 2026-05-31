from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from bijux_phylogenetics.engines.inference import (
    FastaToTreeModelRow,
    FastaToTreeSupportRow,
    infer_unaligned_sequence_type,
    run_fasta_to_tree_workflow,
    write_fasta_to_tree_log,
    write_fasta_to_tree_model_table,
    write_fasta_to_tree_support_table,
)
from bijux_phylogenetics.runtime.errors import (
    EngineWorkflowError,
    InvalidAlignmentError,
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


def _fake_trimal_empty_output(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

if "--version" in sys.argv:
    print("trimAl v2.0")
    raise SystemExit(0)

args = sys.argv[1:]
output_path = Path(args[args.index("-out") + 1])
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text("", encoding="utf-8")
print("warning: trimal fixture wrote an empty alignment", file=sys.stderr)
""",
    )


def _fake_trimal_failure(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys

if "--version" in sys.argv:
    print("trimAl v2.0")
    raise SystemExit(0)

print("trimAl fixture failed deliberately", file=sys.stderr)
raise SystemExit(9)
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
    criteria_lines = (
        " No. Model         -LnL         df  AIC          AICc         BIC\\n"
        + (
            "  1  LG+G          123.456      9   264.912      324.912      266.912\\n"
            "  2  WAG+G         124.000      9   266.000      326.000      268.000\\n"
            "  3  JTT+G         125.500      9   269.000      329.000      271.000\\n"
            "Akaike Information Criterion:           LG+G\\n"
            "Corrected Akaike Information Criterion: LG+G\\n"
            "Bayesian Information Criterion:         LG+G\\n"
            "Best-fit model: LG+G chosen according to BIC\\n"
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
        + "Log-likelihood of the tree: -123.456\\nWARNING: model search used a fixture backend\\n",
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
        f"Best-fit model: {selected_model}\\nLog-likelihood of the tree: -234.567\\nBootstrap analysis completed\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture bootstrap log\\nBEST SCORE FOUND : -234.567\\n",
        encoding="utf-8",
    )
    print("warning: iqtree fixture bootstrap", file=sys.stderr)
    raise SystemExit(0)

prefix.with_suffix(".treefile").write_text(
    "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
    encoding="utf-8",
)
prefix.with_suffix(".iqtree").write_text(
    f"Best-fit model: {selected_model}\\nLog-likelihood of the tree: -345.678\\nTree inference completed\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture inference log\\nBEST SCORE FOUND : -345.678\\n",
    encoding="utf-8",
)
print("warning: iqtree fixture tree inference", file=sys.stderr)
""",
    )


def _fake_iqtree_failure(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys

if "--version" in sys.argv:
    print("IQ-TREE multicore version 2.9.9")
    raise SystemExit(0)

print("iqtree fixture failed deliberately", file=sys.stderr)
raise SystemExit(11)
""",
    )


def _fake_iqtree_missing_support_tree(path: Path) -> Path:
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
        "  1  GTR+G         123.456      12  270.912      330.912      272.912\\n"
        "Akaike Information Criterion:           GTR+G\\n"
        "Corrected Akaike Information Criterion: GTR+G\\n"
        "Bayesian Information Criterion:         GTR+G\\n"
        "Best-fit model according to BIC: GTR+G\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture model-selection log\\nBEST SCORE FOUND : -123.456\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".model").write_text(
        "Best-fit model: GTR+G\\n",
        encoding="utf-8",
    )
    raise SystemExit(0)

if "-bb" in args:
    prefix.with_suffix(".ufboot").write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".iqtree").write_text(
        "Best-fit model: GTR+G\\nLog-likelihood of the tree: -234.567\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture bootstrap log\\nBEST SCORE FOUND : -234.567\\n",
        encoding="utf-8",
    )
    print("warning: iqtree fixture omitted the support tree", file=sys.stderr)
    raise SystemExit(0)

prefix.with_suffix(".treefile").write_text(
    "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
    encoding="utf-8",
)
prefix.with_suffix(".iqtree").write_text(
    "Best-fit model: GTR+G\\nLog-likelihood of the tree: -345.678\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture inference log\\nBEST SCORE FOUND : -345.678\\n",
    encoding="utf-8",
)
""",
    )


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


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


def test_write_fasta_to_tree_model_table_can_render_paths_relative_to_root(
    tmp_path: Path,
) -> None:
    workflow_root = tmp_path / "workflow"
    model_path = workflow_root / "example.model.tsv"

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
                alignment_path=workflow_root / "example.aln",
                trimmed_alignment_path=workflow_root / "example.trimmed.aln",
                manifest_path=workflow_root
                / "engine-artifacts"
                / "example"
                / "model-selection"
                / "model-selection.manifest.json",
            )
        ],
        root_dir=workflow_root,
    )

    assert model_path.read_text(encoding="utf-8") == (
        "workflow\tengine_name\tsequence_type\tselected_model\t"
        "report_selected_model\tartifact_selected_model\tmodel_consistent\t"
        "alignment_path\ttrimmed_alignment_path\tmanifest_path\n"
        "model-selection\tiqtree\tdna\tGTR+G\tGTR+G\tGTR+G\ttrue\t"
        "example.aln\texample.trimmed.aln\t"
        "engine-artifacts/example/model-selection/model-selection.manifest.json\n"
    )


@pytest.mark.slow
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
            bootstrap_replicates=1000,
        )

        assert report.sequence_type == expected_sequence_type
        assert report.selected_model == expected_model
        assert report.method_tier.tier == "supported"
        assert any(
            basis.startswith("real-engine-validation:")
            for basis in report.method_tier.validation_basis
        )
        assert report.output_paths["alignment"].suffix == ".aln"
        assert report.output_paths["trimmed_alignment"].name.endswith(".trimmed.aln")
        assert report.output_paths["tree"].suffix == ".tree"
        assert report.output_paths["log"].suffix == ".log"
        assert report.output_paths["methods_summary"].name.endswith(
            ".methods-summary.md"
        )
        assert report.output_paths["model_table"].name.endswith(".model.tsv")
        assert report.output_paths["support_table"].name.endswith(".support.tsv")
        assert report.output_paths["manifest"].name.endswith(".manifest.json")
        assert report.output_paths["run_manifest"].name.endswith(".run.json")
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
        methods_text = report.output_paths["methods_summary"].read_text(
            encoding="utf-8"
        )
        assert "selected_model:" in log_text
        assert "iqtree random seed: 1" in log_text
        assert "iqtree threads: 1" in log_text
        assert "warning: iqtree fixture bootstrap" in log_text
        assert "Tree Inference Methods Summary" in methods_text
        assert "- selected substitution model:" in methods_text
        model_text = report.output_paths["model_table"].read_text(encoding="utf-8")
        assert "engine-artifacts/" in model_text
        assert str(report.out_dir) not in model_text
        workflow_manifest = _load_json(report.manifest_path)
        run_manifest = _load_json(report.run_manifest_path)
        assert report.workflow == "fasta-to-tree"
        assert workflow_manifest["selected_model"] == expected_model
        assert workflow_manifest["workflow"] == "fasta-to-tree"
        assert workflow_manifest["iqtree_seed"] == 1
        assert workflow_manifest["bootstrap_replicates"] == 1000
        assert workflow_manifest["config"]["iqtree_threads"] == 1
        assert workflow_manifest["output_paths"]["methods_summary"].endswith(
            ".methods-summary.md"
        )
        assert "alignment" in workflow_manifest["commands"]
        assert "iqtree_model_selection" in workflow_manifest["engine_versions"]
        assert report.run_manifest_path == report.output_paths["run_manifest"]
        assert run_manifest["command"] == "run_fasta_to_tree_workflow"
        assert str(report.manifest_path) in run_manifest["output_paths"]
        assert str(input_path) in run_manifest["input_paths"]
        assert str(report.output_paths["tree"]) in run_manifest["output_checksums"]


def test_run_fasta_to_tree_workflow_records_stage_fingerprints(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    report = run_fasta_to_tree_workflow(
        fixture("alignments/example_sequences_raw.fasta"),
        out_dir=tmp_path / "stage-fingerprints",
        prefix="stage-fingerprints",
        mafft_executable=mafft,
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        bootstrap_replicates=1000,
    )

    assert list(report.stage_fingerprints) == [
        "fasta_validation",
        "alignment",
        "trimming",
        "model_selection",
        "inference",
        "support",
        "report",
    ]
    assert report.stage_fingerprints["alignment"].engine_versions == {
        "mafft": "mafft v7.999"
    }
    assert report.stage_fingerprints["support"].engine_versions == {
        "iqtree_bootstrap_support": "IQ-TREE multicore version 2.9.9"
    }
    assert (
        report.stage_fingerprints["report"].upstream_fingerprints["support"]
        == report.stage_fingerprints["support"].fingerprint
    )
    assert all(
        len(stage.fingerprint) == 64 for stage in report.stage_fingerprints.values()
    )
    manifest_payload = _load_json(report.manifest_path)
    assert manifest_payload["stage_fingerprints"]["fasta_validation"]["stage"] == (
        "fasta_validation"
    )
    assert (
        manifest_payload["stage_fingerprints"]["report"]["upstream_fingerprints"][
            "inference"
        ]
        == report.stage_fingerprints["inference"].fingerprint
    )


def test_run_fasta_to_tree_workflow_reruns_only_support_stage_when_replicates_change(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_sequences_raw.fasta")

    first = run_fasta_to_tree_workflow(
        input_path,
        out_dir=tmp_path / "support-rerun",
        prefix="support-rerun",
        mafft_executable=mafft,
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        bootstrap_replicates=1000,
    )
    second = run_fasta_to_tree_workflow(
        input_path,
        out_dir=tmp_path / "support-rerun",
        prefix="support-rerun",
        mafft_executable=mafft,
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        bootstrap_replicates=2000,
        resume=True,
    )

    assert first.alignment_workflow.resumed is False
    assert second.alignment_workflow.resumed is True
    assert second.trimming_workflow.resumed is True
    assert second.model_selection_workflow.resumed is True
    assert second.maximum_likelihood_workflow.resumed is True
    assert second.bootstrap_workflow.resumed is False
    assert (
        first.stage_fingerprints["support"].fingerprint
        != second.stage_fingerprints["support"].fingerprint
    )
    assert (
        first.stage_fingerprints["inference"].fingerprint
        == second.stage_fingerprints["inference"].fingerprint
    )


@pytest.mark.slow
def test_run_fasta_to_tree_workflow_input_change_invalidates_downstream_stages(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = tmp_path / "changed-input.fasta"
    input_path.write_text(">A\nACTG\n>B\nACTGA\n>C\nACT\n", encoding="utf-8")

    first = run_fasta_to_tree_workflow(
        input_path,
        out_dir=tmp_path / "changed-input",
        prefix="changed-input",
        mafft_executable=mafft,
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        bootstrap_replicates=1000,
    )
    input_path.write_text(">A\nTTTTT\n>B\nACTGA\n>C\nGGG\n", encoding="utf-8")
    second = run_fasta_to_tree_workflow(
        input_path,
        out_dir=tmp_path / "changed-input",
        prefix="changed-input",
        mafft_executable=mafft,
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        bootstrap_replicates=1000,
        resume=True,
    )

    assert second.alignment_workflow.resumed is False
    assert second.trimming_workflow.resumed is False
    assert second.model_selection_workflow.resumed is False
    assert second.maximum_likelihood_workflow.resumed is False
    assert second.bootstrap_workflow.resumed is False
    assert (
        first.stage_fingerprints["fasta_validation"].fingerprint
        != second.stage_fingerprints["fasta_validation"].fingerprint
    )
    assert (
        first.stage_fingerprints["alignment"].fingerprint
        != second.stage_fingerprints["alignment"].fingerprint
    )


def test_write_fasta_to_tree_log_renders_workflow_outputs_relative_to_root(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    report = run_fasta_to_tree_workflow(
        fixture("alignments/example_sequences_raw.fasta"),
        out_dir=tmp_path / "portable-log",
        prefix="portable-log",
        mafft_executable=mafft,
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        bootstrap_replicates=1000,
    )

    log_path = tmp_path / "portable-log.rendered.log"
    write_fasta_to_tree_log(log_path, report, root_dir=report.out_dir)
    log_text = log_path.read_text(encoding="utf-8")

    assert (
        "manifest: engine-artifacts/portable-log/model-selection/model-selection.manifest.json"
        in log_text
    )
    assert (
        "output.alignment: engine-artifacts/portable-log/alignment/alignment.aln"
        in log_text
    )
    assert "run_manifest_path: portable-log.run.json" in log_text
    assert str(report.out_dir) not in log_text


def test_run_fasta_to_tree_workflow_log_preserves_relative_input_spelling(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo_root = Path(__file__).resolve().parents[3]
    relative_input = Path(
        "packages/bijux-phylogenetics/tests/fixtures/alignments/example_sequences_raw.fasta"
    )
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    monkeypatch.chdir(repo_root)

    report = run_fasta_to_tree_workflow(
        relative_input,
        out_dir=tmp_path / "relative-input-log",
        prefix="relative-input-log",
        mafft_executable=mafft,
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        bootstrap_replicates=1000,
    )

    log_text = report.output_paths["log"].read_text(encoding="utf-8")
    assert f"input_path: {relative_input}" in log_text
    assert f"command: {mafft} --auto {relative_input}" in log_text


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
        ">Alpha_sample\nACTGACTG\n>rare_taxon\nACTGACTGACTGACTGACTGACTG\n"
    )
    assert any(
        "duplicate sequence identifiers" in warning for warning in report.warnings
    )
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


def test_run_fasta_to_tree_workflow_rejects_ufboot_replicates_below_iqtree_minimum(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    with pytest.raises(
        EngineWorkflowError,
        match="ultrafast bootstrap requires at least 1000 replicates",
    ):
        run_fasta_to_tree_workflow(
            fixture("alignments/example_sequences_raw.fasta"),
            out_dir=tmp_path / "workflow-too-few-bootstraps",
            prefix="workflow-too-few-bootstraps",
            mafft_executable=mafft,
            trimal_executable=trimal,
            iqtree_executable=iqtree,
            bootstrap_replicates=999,
        )


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
        bootstrap_replicates=1000,
    )

    assert report.alignment_workflow.run.command[1:-1] == [
        "--localpair",
        "--maxiterate",
        "1000",
    ]
    assert report.alignment_workflow.notes[0] == "mafft alignment mode: linsi"
    assert "mafft alignment mode: linsi" in report.notes
    assert "command:" in report.output_paths["log"].read_text(encoding="utf-8")


def test_run_fasta_to_tree_workflow_rejects_empty_trimmed_alignment(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal_empty_output(tmp_path / "trimal-empty-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    with pytest.raises(EngineWorkflowError) as error:
        run_fasta_to_tree_workflow(
            fixture("alignments/example_sequences_raw.fasta"),
            out_dir=tmp_path / "empty-trimmed-alignment",
            prefix="empty-trimmed-alignment",
            mafft_executable=mafft,
            trimal_executable=trimal,
            iqtree_executable=iqtree,
            bootstrap_replicates=1000,
        )
    assert error.value.code == "engine_output_empty"
    assert error.value.details["output_name"] == "trimmed_alignment"


def test_run_fasta_to_tree_workflow_surfaces_trimming_failure(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal_failure(tmp_path / "trimal-failure-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    with pytest.raises(
        EngineWorkflowError, match="trimal alignment-trimming failed with exit code 9"
    ):
        run_fasta_to_tree_workflow(
            fixture("alignments/example_sequences_raw.fasta"),
            out_dir=tmp_path / "trimming-failure",
            prefix="trimming-failure",
            mafft_executable=mafft,
            trimal_executable=trimal,
            iqtree_executable=iqtree,
            bootstrap_replicates=1000,
        )


def test_run_fasta_to_tree_workflow_surfaces_tree_inference_failure(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree_failure(tmp_path / "iqtree-failure-fixture")

    with pytest.raises(
        EngineWorkflowError, match="iqtree model-selection failed with exit code 11"
    ):
        run_fasta_to_tree_workflow(
            fixture("alignments/example_sequences_raw.fasta"),
            out_dir=tmp_path / "inference-failure",
            prefix="inference-failure",
            mafft_executable=mafft,
            trimal_executable=trimal,
            iqtree_executable=iqtree,
            bootstrap_replicates=1000,
        )


def test_run_fasta_to_tree_workflow_rejects_missing_support_tree_output(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree_missing_support_tree(
        tmp_path / "iqtree-missing-support-tree-fixture"
    )

    with pytest.raises(
        EngineWorkflowError,
        match="iqtree bootstrap-support did not produce expected outputs",
    ):
        run_fasta_to_tree_workflow(
            fixture("alignments/example_sequences_raw.fasta"),
            out_dir=tmp_path / "missing-support-tree",
            prefix="missing-support-tree",
            mafft_executable=mafft,
            trimal_executable=trimal,
            iqtree_executable=iqtree,
            bootstrap_replicates=1000,
        )


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
        bootstrap_replicates=1000,
    )

    assert report.trimming_workflow.run.command[5:] == ["-strictplus"]
    assert report.trimming_workflow.trimming_summary is not None
    assert report.trimming_workflow.trimming_summary.mode == "strictplus"
    assert report.trimming_workflow.trimming_summary.removed_site_count == 3
    assert "trimal trimming mode: strictplus" in report.notes


def test_run_fasta_to_tree_workflow_passes_deterministic_iqtree_controls(
    tmp_path: Path,
) -> None:
    mafft = _fake_mafft(tmp_path / "mafft-fixture")
    trimal = _fake_trimal(tmp_path / "trimal-fixture")
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")

    report = run_fasta_to_tree_workflow(
        fixture("alignments/example_sequences_raw.fasta"),
        out_dir=tmp_path / "workflow-deterministic-controls",
        prefix="workflow-deterministic-controls",
        mafft_executable=mafft,
        trimal_executable=trimal,
        iqtree_executable=iqtree,
        iqtree_seed=7,
        iqtree_threads=3,
        bootstrap_replicates=1000,
    )

    for workflow in (
        report.model_selection_workflow,
        report.maximum_likelihood_workflow,
        report.bootstrap_workflow,
    ):
        seed_index = workflow.run.command.index("-seed")
        thread_index = workflow.run.command.index("-nt")
        assert workflow.run.command[seed_index : seed_index + 2] == ["-seed", "7"]
        assert workflow.run.command[thread_index : thread_index + 2] == ["-nt", "3"]
    assert "iqtree random seed: 7" in report.notes
    assert "iqtree threads: 3" in report.notes
