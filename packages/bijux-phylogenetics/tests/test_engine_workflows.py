from __future__ import annotations

import os
from pathlib import Path
import shutil

import pytest

from bijux_phylogenetics.engines import (
    build_inference_sensitivity_report,
    build_model_selection_limitations_report,
    bundle_inference_workflow_evidence,
    compare_fast_and_ml_trees,
    list_mafft_alignment_modes,
    list_trimal_trimming_modes,
    render_inference_sensitivity_report,
    render_inference_workflow_report,
    render_model_selection_limitations_report,
    resolve_mafft_alignment_mode,
    resolve_trimal_trimming_mode,
    run_alignment_trimming,
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_codon_aware_multiple_sequence_alignment,
    run_fast_tree_inference,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
    run_multiple_sequence_alignment,
    run_sh_alrt_support_estimation,
)
from bijux_phylogenetics.errors import EngineWorkflowError, InvalidAlignmentError
from bijux_phylogenetics.io.fasta import load_fasta_alignment

FIXTURES = Path(__file__).parent / "fixtures"
REPOSITORY_ROOT = Path(__file__).resolve().parents[3]


def fixture(path: str) -> Path:
    return FIXTURES / path


def _real_mafft_executable() -> Path | None:
    configured = os.environ.get("BIJUX_PHYLOGENETICS_MAFFT_EXECUTABLE")
    if configured:
        candidate = Path(configured)
        if candidate.exists():
            return candidate
    resolved = shutil.which("mafft")
    if resolved is not None:
        return Path(resolved)
    artifact_candidate = (
        REPOSITORY_ROOT / "artifacts" / "mafft" / "mafft-mac" / "mafft.bat"
    )
    if artifact_candidate.exists():
        return artifact_candidate
    return None


def _real_trimal_executable() -> Path | None:
    configured = os.environ.get("BIJUX_PHYLOGENETICS_TRIMAL_EXECUTABLE")
    if configured:
        candidate = Path(configured)
        if candidate.exists():
            return candidate
    resolved = shutil.which("trimal")
    if resolved is not None:
        return Path(resolved)
    artifact_candidate = (
        REPOSITORY_ROOT / "artifacts" / "trimal" / "trimal" / "source" / "trimal"
    )
    if artifact_candidate.exists():
        return artifact_candidate
    return None


def _real_iqtree_executable() -> Path | None:
    configured = os.environ.get("BIJUX_PHYLOGENETICS_IQTREE_EXECUTABLE")
    if configured:
        candidate = Path(configured)
        if candidate.exists():
            return candidate
    resolved = shutil.which("iqtree2")
    if resolved is not None:
        return Path(resolved)
    return None


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


def _fake_trimal_empty(path: Path) -> Path:
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
input_path = Path(args[args.index("-in") + 1])
records = []
identifier = None
for raw_line in input_path.read_text(encoding="utf-8").splitlines():
    line = raw_line.strip()
    if line.startswith(">"):
        identifier = line[1:]
        records.append(identifier)
output_path.parent.mkdir(parents=True, exist_ok=True)
with output_path.open("w", encoding="utf-8") as handle:
    for identifier in records:
        handle.write(f">{identifier}\\n\\n")
print("warning: trimal fixture removed every site", file=sys.stderr)
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

if "-con" in args:
    prefix.with_suffix(".contree").write_text("((A:0.1,B:0.1)90:0.2,(C:0.1,D:0.1)85:0.2);\\n", encoding="utf-8")
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture consensus log\\n",
        encoding="utf-8",
    )
    print("warning: iqtree fixture consensus", file=sys.stderr)
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
    prefix.with_suffix(".treefile").write_text("((A:0.1,B:0.1)95:0.2,(C:0.1,D:0.1)88:0.2);\\n", encoding="utf-8")
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


def _fake_iqtree_fail(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys

if "--version" in sys.argv[1:]:
    print("IQ-TREE multicore version 2.9.9")
    raise SystemExit(0)

print("warning: fixture is about to fail", file=sys.stderr)
print("fatal: inference failed", file=sys.stderr)
raise SystemExit(3)
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


def _fake_iqtree_tree_variant(path: Path, *, tree_newick: str) -> Path:
    return _write_executable(
        path,
        f"""#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("IQ-TREE multicore version 2.9.9")
    raise SystemExit(0)

prefix = Path(args[args.index("-pre") + 1])
prefix.parent.mkdir(parents=True, exist_ok=True)
prefix.with_suffix(".treefile").write_text({tree_newick!r} + "\\n", encoding="utf-8")
prefix.with_suffix(".iqtree").write_text(
    "Best-fit model: GTR+G\\nLog-likelihood of the tree: -456.789\\nTree inference completed\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture inference log\\nBEST SCORE FOUND : -456.789\\n",
    encoding="utf-8",
)
print("warning: iqtree fixture tree inference", file=sys.stderr)
""",
    )


def test_run_multiple_sequence_alignment_captures_logs_version_and_manifest(
    tmp_path: Path,
) -> None:
    executable = _fake_mafft(tmp_path / "mafft-fixture")
    input_path = tmp_path / "unaligned.fasta"
    input_path.write_text(">A\nACTG\n>B\nACTGA\n>C\nACT\n", encoding="utf-8")
    output_path = tmp_path / "aligned.fasta"

    report = run_multiple_sequence_alignment(
        input_path, output_path, executable=executable
    )

    records = load_fasta_alignment(output_path)
    assert [len(record.sequence) for record in records] == [5, 5, 5]
    assert "mafft v7.999" in report.run.version.text
    assert report.run.command[0] == str(executable)
    assert report.run.warning_lines == [
        "WARNING: mafft fixture inserted alignment padding"
    ]
    assert report.manifest_path.exists()


def test_run_codon_aware_multiple_sequence_alignment_preserves_triplet_gaps(
    tmp_path: Path,
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

    report = run_codon_aware_multiple_sequence_alignment(
        input_path,
        output_path,
        executable=executable,
        mode="linsi",
    )

    records = load_fasta_alignment(output_path)
    assert [(record.identifier, record.sequence) for record in records] == [
        ("long_good", "ATGGAATGGAAA"),
        ("short_good", "ATGGAATGG---"),
    ]
    assert report.sequence_type == "dna"
    assert report.accepted_sequence_count == 2
    assert [row.identifier for row in report.excluded_sequences] == [
        "frameshift",
        "internal_stop",
    ]
    assert report.output_paths["guide_input"].read_text(encoding="utf-8") == (
        ">long_good\nMEWK\n>short_good\nMEW\n"
    )
    assert report.output_paths["excluded_sequences"].read_text(encoding="utf-8") == (
        "identifier\tcomparable_length\treason\tpremature_stop_count\t"
        "terminal_stop_count\ttrailing_bases\tnote\n"
        "frameshift\t8\tframe-error\t0\t0\t2\t"
        "sequence is not frame-consistent after gaps and missing data are removed\n"
        "internal_stop\t9\tinternal-stop-codon\t1\t0\t0\t"
        "sequence contains one or more premature stop codons\n"
    )
    assert report.run.command[1:-1] == ["--localpair", "--maxiterate", "1000"]
    assert report.notes[0].startswith("codon-aware alignment preserved")
    assert report.manifest_path.exists()


def test_run_codon_aware_multiple_sequence_alignment_fails_when_every_sequence_is_excluded(
    tmp_path: Path,
) -> None:
    executable = _fake_mafft(tmp_path / "mafft-fixture")
    input_path = tmp_path / "coding-invalid.fasta"
    input_path.write_text(
        ">frameshift\nATGGAATG\n>internal_stop\nATGTAGTGG\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "codon-aligned.fasta"

    with pytest.raises(
        InvalidAlignmentError, match="excluded every sequence"
    ) as error_info:
        run_codon_aware_multiple_sequence_alignment(
            input_path,
            output_path,
            executable=executable,
        )

    assert "excluded every sequence" in str(error_info.value)


def test_mafft_alignment_modes_resolve_to_explicit_documented_arguments() -> None:
    assert list_mafft_alignment_modes() == (
        "auto",
        "linsi",
        "ginsi",
        "einsi",
        "fast",
    )
    assert resolve_mafft_alignment_mode("auto") == ("--auto",)
    assert resolve_mafft_alignment_mode("linsi") == (
        "--localpair",
        "--maxiterate",
        "1000",
    )
    assert resolve_mafft_alignment_mode("ginsi") == (
        "--globalpair",
        "--maxiterate",
        "1000",
    )
    assert resolve_mafft_alignment_mode("einsi") == (
        "--ep",
        "0",
        "--genafpair",
        "--maxiterate",
        "1000",
    )
    assert resolve_mafft_alignment_mode("fast") == (
        "--retree",
        "2",
        "--maxiterate",
        "0",
    )


def test_run_multiple_sequence_alignment_supports_all_named_mafft_modes(
    tmp_path: Path,
) -> None:
    executable = _fake_mafft(tmp_path / "mafft-fixture")
    input_path = tmp_path / "unaligned.fasta"
    input_path.write_text(">A\nACTG\n>B\nACTGA\n>C\nACT\n", encoding="utf-8")

    expected_prefixes = {
        "auto": ["--auto"],
        "linsi": ["--localpair", "--maxiterate", "1000"],
        "ginsi": ["--globalpair", "--maxiterate", "1000"],
        "einsi": ["--ep", "0", "--genafpair", "--maxiterate", "1000"],
        "fast": ["--retree", "2", "--maxiterate", "0"],
    }

    for mode, expected_args in expected_prefixes.items():
        output_path = tmp_path / f"{mode}.fasta"
        report = run_multiple_sequence_alignment(
            input_path,
            output_path,
            executable=executable,
            mode=mode,
        )

        assert report.run.command[1:-1] == expected_args
        assert report.notes[0] == f"mafft alignment mode: {mode}"
        assert load_fasta_alignment(output_path)


def test_run_multiple_sequence_alignment_with_real_mafft_on_small_dataset(
    tmp_path: Path,
) -> None:
    executable = _real_mafft_executable()
    if executable is None:
        pytest.skip("real MAFFT executable is not available for integration coverage")

    output_path = tmp_path / "real-mafft-alignment.fasta"
    report = run_multiple_sequence_alignment(
        fixture("alignments/example_sequences_raw.fasta"),
        output_path,
        executable=executable,
        mode="linsi",
    )

    alignment = load_fasta_alignment(output_path)
    widths = {len(record.sequence) for record in alignment}

    assert report.run.exit_code == 0
    assert report.run.command[1:-1] == ["--localpair", "--maxiterate", "1000"]
    assert "v7." in report.run.version.text
    assert report.manifest_path.exists()
    assert len(alignment) == 4
    assert len(widths) == 1


def test_run_alignment_trimming_with_real_trimal_on_small_dataset(
    tmp_path: Path,
) -> None:
    executable = _real_trimal_executable()
    if executable is None:
        pytest.skip("real trimAl executable is not available for integration coverage")

    input_path = fixture("alignments/example_alignment_trim.fasta")
    output_path = tmp_path / "real-trimal-trimmed.fasta"
    report = run_alignment_trimming(
        input_path,
        output_path,
        executable=executable,
        mode="gap-threshold",
        gap_threshold=0.8,
    )

    input_records = load_fasta_alignment(input_path)
    trimmed_records = load_fasta_alignment(output_path)

    assert report.run.exit_code == 0
    assert report.run.command[5:] == ["-gt", "0.800000"]
    assert "trimal v" in report.run.version.text.lower()
    assert report.manifest_path.exists()
    assert report.trimming_summary is not None
    assert report.trimming_summary.removed_site_count >= 1
    assert report.trimming_summary.trimmed_alignment_length == len(
        trimmed_records[0].sequence
    )
    assert len(trimmed_records[0].sequence) < len(input_records[0].sequence)


def test_run_iqtree_backend_with_real_executable_on_small_dataset(
    tmp_path: Path,
) -> None:
    executable = _real_iqtree_executable()
    if executable is None:
        pytest.skip("real IQ-TREE executable is not available for integration coverage")

    input_path = fixture("alignments/example_alignment.fasta")
    model_report = run_model_selection(
        input_path,
        out_dir=tmp_path / "model",
        executable=executable,
        prefix="real",
        sequence_type="dna",
    )
    assert model_report.selected_model is not None
    ml_report = run_maximum_likelihood_tree_inference(
        input_path,
        out_dir=tmp_path / "ml",
        model=model_report.selected_model,
        executable=executable,
        prefix="real",
        sequence_type="dna",
    )
    bootstrap_report = run_bootstrap_support_estimation(
        input_path,
        out_dir=tmp_path / "bootstrap",
        model=model_report.selected_model,
        executable=executable,
        prefix="real",
        sequence_type="dna",
        replicates=1000,
    )
    sh_alrt_report = run_sh_alrt_support_estimation(
        input_path,
        out_dir=tmp_path / "sh-alrt",
        model=model_report.selected_model,
        executable=executable,
        prefix="real",
        sequence_type="dna",
        sh_alrt_replicates=1000,
        bootstrap_replicates=1000,
    )
    consensus_report = run_bootstrap_consensus_tree(
        bootstrap_report.output_paths["bootstrap_trees"],
        out_dir=tmp_path / "consensus",
        executable=executable,
        prefix="real",
    )

    assert model_report.output_paths["iqtree_report"].exists()
    assert model_report.output_paths["iqtree_log"].exists()
    assert model_report.output_paths["model_candidates"].exists()
    assert model_report.log_likelihood is not None
    assert model_report.model_selection_summary is not None
    assert model_report.model_selection_summary.candidate_count >= 1
    assert model_report.model_selection_summary.best_model_bic is not None
    assert ml_report.output_paths["tree"].exists()
    assert ml_report.output_paths["iqtree_log"].exists()
    assert ml_report.log_likelihood is not None
    assert bootstrap_report.output_paths["bootstrap_trees"].exists()
    assert bootstrap_report.output_paths["support_tree"].exists()
    assert bootstrap_report.output_paths["support_table"].exists()
    assert bootstrap_report.output_paths["low_support_branches"].exists()
    assert bootstrap_report.output_paths["support_histogram"].exists()
    assert bootstrap_report.output_paths["iqtree_log"].exists()
    assert bootstrap_report.log_likelihood is not None
    assert bootstrap_report.iqtree_summary is not None
    assert bootstrap_report.iqtree_summary.support_value_count >= 1
    assert bootstrap_report.bootstrap_support_summary is not None
    assert bootstrap_report.weak_backbone_report is not None
    assert sh_alrt_report.output_paths["support_tree"].exists()
    assert sh_alrt_report.output_paths["bootstrap_trees"].exists()
    assert sh_alrt_report.output_paths["support_table"].exists()
    assert sh_alrt_report.output_paths["conflicting_support_branches"].exists()
    assert sh_alrt_report.iqtree_summary is not None
    assert sh_alrt_report.iqtree_summary.support_value_count >= 1
    assert sh_alrt_report.sh_alrt_support_summary is not None
    assert sh_alrt_report.sh_alrt_support_summary.annotated_node_count >= 1
    assert consensus_report.output_paths["consensus_tree"].exists()
    assert consensus_report.output_paths["iqtree_log"].exists()
    assert consensus_report.iqtree_summary is not None
    assert consensus_report.iqtree_summary.support_value_count >= 1


def test_run_bootstrap_support_estimation_exports_branch_ledgers_and_histogram(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")

    report = run_bootstrap_support_estimation(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "bootstrap",
        model="GTR+G",
        executable=executable,
        prefix="example",
        replicates=1000,
    )

    assert report.output_paths["support_table"].read_text(encoding="utf-8") == (
        "node\tdescendant_taxa\tsupport\tsupport_fraction\tis_backbone\t"
        "support_bucket\tlow_support\n"
        "A|B\tA,B\t95\t0.95\ttrue\tge90\tfalse\n"
        "C|D\tC,D\t88\t0.88\ttrue\t70to89\tfalse\n"
    )
    assert (
        report.output_paths["low_support_branches"].read_text(encoding="utf-8")
        == "node\tdescendant_taxa\tsupport\tsupport_fraction\tis_backbone\t"
        "support_bucket\tlow_support\n"
    )
    assert report.output_paths["support_histogram"].read_text(encoding="utf-8") == (
        "support_bucket\tminimum_support\tmaximum_support\tnode_count\n"
        "lt50\t\t50\t0\n"
        "50to69\t50\t70\t0\n"
        "70to89\t70\t90\t1\n"
        "ge90\t90\t\t1\n"
    )
    assert report.bootstrap_support_summary is not None
    assert report.bootstrap_support_summary.weakly_supported_clade_count == 0
    assert report.weak_backbone_report is not None
    assert report.weak_backbone_report.weak_backbone_node_count == 0


def test_run_sh_alrt_support_estimation_exports_combined_support_and_conflicts(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")

    report = run_sh_alrt_support_estimation(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "sh-alrt",
        model="GTR+G",
        executable=executable,
        prefix="example",
        sh_alrt_replicates=1000,
        bootstrap_replicates=1000,
    )

    assert report.output_paths["support_tree"].exists()
    assert report.output_paths["bootstrap_trees"].exists()
    assert report.output_paths["support_table"].read_text(encoding="utf-8") == (
        "node\tdescendant_taxa\tsh_alrt_support\tsh_alrt_support_fraction\t"
        "ufboot_support\tufboot_support_fraction\tis_backbone\tsh_alrt_strong\t"
        "ufboot_strong\tconflicting_support_signal\tsupport_agreement\n"
        "A|B\tA,B\t82\t0.82\t97\t0.97\ttrue\ttrue\ttrue\tfalse\tboth_strong\n"
        "C|D\tC,D\t79\t0.79\t96\t0.96\ttrue\tfalse\ttrue\ttrue\tufboot_only\n"
    )
    assert (
        report.output_paths["conflicting_support_branches"].read_text(
            encoding="utf-8"
        )
        == "node\tdescendant_taxa\tsh_alrt_support\tsh_alrt_support_fraction\t"
        "ufboot_support\tufboot_support_fraction\tis_backbone\tsh_alrt_strong\t"
        "ufboot_strong\tconflicting_support_signal\tsupport_agreement\n"
        "C|D\tC,D\t79\t0.79\t96\t0.96\ttrue\tfalse\ttrue\ttrue\tufboot_only\n"
    )
    assert report.iqtree_summary is not None
    assert report.iqtree_summary.support_value_count == 2
    assert report.sh_alrt_support_summary is not None
    assert report.sh_alrt_support_summary.conflicting_support_signal_count == 1


def test_trimal_trimming_modes_resolve_to_explicit_documented_arguments() -> None:
    assert list_trimal_trimming_modes() == (
        "gap-threshold",
        "gappyout",
        "strict",
        "strictplus",
        "automated1",
    )
    assert resolve_trimal_trimming_mode("gap-threshold", gap_threshold=0.2) == (
        "-gt",
        "0.200000",
    )
    assert resolve_trimal_trimming_mode("gappyout", gap_threshold=0.2) == (
        "-gappyout",
    )
    assert resolve_trimal_trimming_mode("strict", gap_threshold=0.2) == ("-strict",)
    assert resolve_trimal_trimming_mode("strictplus", gap_threshold=0.2) == (
        "-strictplus",
    )
    assert resolve_trimal_trimming_mode("automated1", gap_threshold=0.2) == (
        "-automated1",
    )


def test_run_alignment_trimming_writes_trimmed_alignment_and_warning_manifest(
    tmp_path: Path,
) -> None:
    executable = _fake_trimal(tmp_path / "trimal-fixture")
    input_path = fixture("alignments/example_alignment_trim.fasta")
    output_path = tmp_path / "trimmed.fasta"

    report = run_alignment_trimming(
        input_path, output_path, executable=executable, gap_threshold=0.2
    )

    records = load_fasta_alignment(output_path)
    input_alignment_length = len(load_fasta_alignment(input_path)[0].sequence)
    assert (
        len(records[0].sequence)
        == input_alignment_length - 1
    )
    assert report.run.warning_lines == [
        "warning: trimal fixture gap-threshold trimmed one trailing site"
    ]
    assert report.trimming_summary is not None
    assert report.trimming_summary.mode == "gap-threshold"
    assert report.trimming_summary.input_alignment_length == input_alignment_length
    assert report.trimming_summary.trimmed_alignment_length == input_alignment_length - 1
    assert report.trimming_summary.retained_site_count == input_alignment_length - 1
    assert report.trimming_summary.removed_site_count == 1
    assert report.notes[1] == (
        f"retained sites: {input_alignment_length - 1} of {input_alignment_length}"
    )
    assert report.manifest_path.exists()


def test_run_alignment_trimming_supports_all_named_trimal_modes(tmp_path: Path) -> None:
    executable = _fake_trimal(tmp_path / "trimal-fixture")
    input_path = fixture("alignments/example_alignment_trim.fasta")
    expected = {
        "gap-threshold": (["-gt", "0.200000"], 1),
        "gappyout": (["-gappyout"], 1),
        "strict": (["-strict"], 2),
        "strictplus": (["-strictplus"], 3),
        "automated1": (["-automated1"], 2),
    }

    for mode, (mode_args, removed_sites) in expected.items():
        output_path = tmp_path / f"{mode}.fasta"
        report = run_alignment_trimming(
            input_path,
            output_path,
            executable=executable,
            mode=mode,
            gap_threshold=0.2,
        )

        assert report.run.command[5:] == mode_args
        assert report.notes[0] == f"trimal trimming mode: {mode}"
        assert report.trimming_summary is not None
        assert report.trimming_summary.mode == mode
        trimmed = load_fasta_alignment(output_path)
        assert len(trimmed[0].sequence) == (
            len(load_fasta_alignment(input_path)[0].sequence) - removed_sites
        )


def test_run_model_selection_parses_best_fit_model_and_writes_manifest(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_alignment.fasta")

    report = run_model_selection(
        input_path, out_dir=tmp_path / "model", executable=executable, prefix="example"
    )

    assert report.selected_model == "GTR+G"
    seed_index = report.run.command.index("-seed")
    thread_index = report.run.command.index("-nt")
    assert report.run.command[seed_index : seed_index + 2] == ["-seed", "1"]
    assert report.run.command[thread_index : thread_index + 2] == ["-nt", "1"]
    selected_model_path = report.output_paths["selected_model"]
    assert selected_model_path.read_text(encoding="utf-8").strip() == "GTR+G"
    assert report.output_paths["model_candidates"].exists()
    assert "HKY+G" in report.output_paths["model_candidates"].read_text(
        encoding="utf-8"
    )
    assert report.run.warning_lines == ["warning: iqtree fixture model selection"]
    assert report.manifest_path.exists()
    assert report.model_selection_summary is not None
    assert report.model_selection_summary.selected_criterion == "BIC"
    assert report.model_selection_summary.best_model_aic == "HKY+G"
    assert report.model_selection_summary.best_model_aicc == "JC"
    assert report.model_selection_summary.best_model_bic == "GTR+G"
    assert report.model_selection_summary.candidate_count == 3


def test_run_model_selection_supports_partitioned_alignment_and_writes_summary(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")

    report = run_model_selection(
        fixture("alignments/example_multilocus_alignment.fasta"),
        out_dir=tmp_path / "model",
        executable=executable,
        prefix="partitioned",
        partition_path=fixture("alignments/example_multilocus_partitions.txt"),
    )

    assert "-p" in report.run.command
    assert str(
        fixture("alignments/example_multilocus_alignment.fasta").resolve()
    ) in report.run.command
    assert report.output_paths["iqtree_log"].exists()
    assert report.selected_model == "GTR+G"
    assert report.log_likelihood == pytest.approx(-123.456)
    assert report.iqtree_summary is not None
    assert report.iqtree_summary.support_value_count == 0
    assert report.model_selection_summary is not None
    assert report.model_selection_summary.candidate_count == 3
    summary_path = report.output_paths["partition_summary"]
    assert summary_path.exists()
    assert "gene_beta" in summary_path.read_text(encoding="utf-8")


def test_run_model_selection_supports_protein_alignment_candidates(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")

    report = run_model_selection(
        fixture("alignments/example_alignment_protein.fasta"),
        out_dir=tmp_path / "protein-model",
        executable=executable,
        prefix="protein",
        sequence_type="protein",
    )

    assert report.selected_model == "LG+G4"
    assert report.model_selection_summary is not None
    assert report.model_selection_summary.selected_criterion == "BIC"
    assert report.model_selection_summary.best_model_bic == "LG+G4"
    assert report.model_selection_summary.candidate_count == 3
    assert "LG+G4" in report.output_paths["model_candidates"].read_text(
        encoding="utf-8"
    )


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
    fast_report = run_fast_tree_inference(
        input_path, tmp_path / "fasttree.nwk", executable=fasttree
    )

    assert ml_report.output_paths["tree"].exists()
    assert bootstrap_report.output_paths["bootstrap_trees"].exists()
    assert consensus_report.output_paths["consensus_tree"].exists()
    assert fast_report.output_paths["tree"].exists()
    ml_seed_index = ml_report.run.command.index("-seed")
    ml_thread_index = ml_report.run.command.index("-nt")
    bootstrap_seed_index = bootstrap_report.run.command.index("-seed")
    bootstrap_thread_index = bootstrap_report.run.command.index("-nt")
    assert ml_report.run.command[ml_seed_index : ml_seed_index + 2] == ["-seed", "1"]
    assert ml_report.run.command[ml_thread_index : ml_thread_index + 2] == ["-nt", "1"]
    assert (
        bootstrap_report.run.command[bootstrap_seed_index : bootstrap_seed_index + 2]
        == ["-seed", "1"]
    )
    assert (
        bootstrap_report.run.command[
            bootstrap_thread_index : bootstrap_thread_index + 2
        ]
        == ["-nt", "1"]
    )
    assert ml_report.run.warning_lines == ["warning: iqtree fixture tree inference"]
    assert bootstrap_report.run.warning_lines == ["warning: iqtree fixture bootstrap"]
    assert consensus_report.run.warning_lines == ["warning: iqtree fixture consensus"]
    assert fast_report.run.warning_lines == [
        "warning: fasttree fixture approximate support only"
    ]
    assert fast_report.fasttree_support_summary is not None
    assert fast_report.fasttree_support_summary.annotated_node_count == 2
    assert fast_report.fasttree_support_summary.minimum_local_support == pytest.approx(
        0.62
    )
    assert fast_report.fasttree_support_summary.maximum_local_support == pytest.approx(
        0.98
    )
    assert fast_report.fasttree_support_summary.support_label_kind == (
        "sh-like-local-support"
    )
    assert fast_report.fasttree_support_summary.support_scale == "proportion-0-to-1"
    assert fast_report.output_paths["support_table"].exists()
    assert fast_report.output_paths["low_support_branches"].exists()
    assert fast_report.output_paths["support_histogram"].exists()
    assert "0.98" in fast_report.output_paths["support_table"].read_text(
        encoding="utf-8"
    )
    assert "0.62" in fast_report.output_paths["low_support_branches"].read_text(
        encoding="utf-8"
    )
    assert any(
        "approximately maximum-likelihood" in note for note in fast_report.notes
    )
    assert ml_report.output_paths["iqtree_log"].exists()
    assert bootstrap_report.output_paths["iqtree_log"].exists()
    assert consensus_report.output_paths["iqtree_log"].exists()
    assert ml_report.log_likelihood == pytest.approx(-345.678)
    assert bootstrap_report.log_likelihood == pytest.approx(-234.567)
    assert consensus_report.log_likelihood is None
    assert ml_report.iqtree_summary is not None
    assert ml_report.iqtree_summary.support_value_count == 0
    assert bootstrap_report.iqtree_summary is not None
    assert bootstrap_report.iqtree_summary.support_value_count == 2
    assert consensus_report.iqtree_summary is not None
    assert consensus_report.iqtree_summary.support_value_count == 2


def test_run_fast_tree_inference_supports_nucleotide_and_protein_modes(
    tmp_path: Path,
) -> None:
    executable = _fake_fasttree(tmp_path / "fasttree-fixture")

    dna_report = run_fast_tree_inference(
        fixture("alignments/example_alignment.fasta"),
        tmp_path / "dna-fasttree.nwk",
        executable=executable,
        sequence_type="dna",
    )
    protein_report = run_fast_tree_inference(
        fixture("alignments/example_alignment_protein.fasta"),
        tmp_path / "protein-fasttree.nwk",
        executable=executable,
        sequence_type="protein",
    )

    assert dna_report.run.command[1:3] == ["-gtr", "-nt"]
    assert dna_report.fasttree_support_summary is not None
    assert dna_report.fasttree_support_summary.annotated_node_count == 2
    assert protein_report.run.command[1] == "-lg"
    assert protein_report.fasttree_support_summary is not None
    assert protein_report.fasttree_support_summary.annotated_node_count == 2


def test_run_ml_and_bootstrap_support_mixed_partition_datatypes(tmp_path: Path) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")
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

    ml_report = run_maximum_likelihood_tree_inference(
        input_path,
        out_dir=tmp_path / "ml",
        model="MFP",
        executable=executable,
        prefix="mixed",
        partition_path=partition_path,
    )
    bootstrap_report = run_bootstrap_support_estimation(
        input_path,
        out_dir=tmp_path / "bootstrap",
        model="MFP",
        executable=executable,
        prefix="mixed",
        partition_path=partition_path,
    )

    assert "-p" in ml_report.run.command
    assert "-s" not in ml_report.run.command
    assert ml_report.output_paths["partition_scheme"].suffix == ".nex"
    assert ml_report.output_paths["partition_summary"].exists()
    assert any(
        key.startswith("partition_alignment_") for key in ml_report.output_paths
    )
    assert "-p" in bootstrap_report.run.command
    assert "-s" not in bootstrap_report.run.command
    assert ml_report.selected_model == "GTR+G"
    assert bootstrap_report.selected_model == "GTR+G"
    assert ml_report.log_likelihood == pytest.approx(-345.678)
    assert bootstrap_report.log_likelihood == pytest.approx(-234.567)


def test_run_ml_rejects_fixed_model_for_mixed_partition_datatypes(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")
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

    with pytest.raises(
        EngineWorkflowError,
        match="mixed DNA/protein partition analyses require a model-selection keyword",
    ):
        run_maximum_likelihood_tree_inference(
            input_path,
            out_dir=tmp_path / "ml",
            model="GTR+G",
            executable=executable,
            prefix="mixed",
            partition_path=partition_path,
        )


def test_inference_workflows_detect_failed_runs_and_empty_filtered_alignments(
    tmp_path: Path,
) -> None:
    failing_iqtree = _fake_iqtree_fail(tmp_path / "iqtree-fail")
    empty_trimal = _fake_trimal_empty(tmp_path / "trimal-empty")
    input_path = fixture("alignments/example_alignment.fasta")

    try:
        run_maximum_likelihood_tree_inference(
            input_path,
            out_dir=tmp_path / "failed-ml",
            model="GTR+G",
            executable=failing_iqtree,
            prefix="example",
        )
    except EngineWorkflowError as error:
        assert "failed with exit code 3" in error.message
    else:  # pragma: no cover - defensive assertion
        raise AssertionError(
            "expected failing IQ-TREE fixture to raise EngineWorkflowError"
        )

    try:
        run_alignment_trimming(
            input_path, tmp_path / "empty.fasta", executable=empty_trimal
        )
    except EngineWorkflowError as error:
        assert "inference alignment is empty after filtering" in error.message
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected empty trim output to raise EngineWorkflowError")


def test_inference_workflow_resume_and_html_report(tmp_path: Path) -> None:
    iqtree = _fake_iqtree(tmp_path / "iqtree-fixture")
    input_path = fixture("alignments/example_alignment.fasta")

    first = run_maximum_likelihood_tree_inference(
        input_path,
        out_dir=tmp_path / "ml",
        model="GTR+G",
        executable=iqtree,
        prefix="example",
        resume=False,
    )
    second = run_maximum_likelihood_tree_inference(
        input_path,
        out_dir=tmp_path / "ml",
        model="GTR+G",
        executable=iqtree,
        prefix="example",
        resume=True,
    )

    assert first.resumed is False
    assert second.resumed is True
    assert second.input_checksums == first.input_checksums
    assert second.output_checksums == first.output_checksums

    html_report = render_inference_workflow_report(
        manifest_path=first.manifest_path,
        out_path=tmp_path / "workflow-report.html",
    )
    assert html_report.output_path.exists()
    assert html_report.report_kind == "inference-workflow"
    assert "workflow-consistency" in html_report.supplement_sections


def test_compare_fast_and_ml_trees_builds_html_report(tmp_path: Path) -> None:
    fast_tree_path = tmp_path / "fast.nwk"
    ml_tree_path = tmp_path / "ml.nwk"
    fast_tree_path.write_text(
        "((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);\n", encoding="utf-8"
    )
    ml_tree_path.write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n", encoding="utf-8"
    )

    comparison = compare_fast_and_ml_trees(
        fast_tree_path, ml_tree_path, out_path=tmp_path / "comparison.html"
    )

    assert comparison.comparison_report.output_path.exists()
    assert comparison.comparison_report.topology.shared_taxa == ["A", "B", "C", "D"]


def test_model_selection_limitations_report_records_interpretation_boundaries(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")
    workflow = run_model_selection(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "model",
        executable=executable,
        prefix="example",
    )

    report = build_model_selection_limitations_report(workflow.manifest_path)
    rendered = render_model_selection_limitations_report(
        manifest_path=workflow.manifest_path,
        out_path=tmp_path / "model-limitations.html",
    )

    assert report.selected_model == "GTR+G"
    assert report.selected_criterion == "BIC"
    assert report.candidate_model_count == 3
    assert report.best_model_aic == "HKY+G"
    assert report.best_model_aicc == "JC"
    assert report.best_model_bic == "GTR+G"
    assert report.limitations
    assert report.interpretation_limits
    assert rendered.output_path.exists()
    assert rendered.report_kind == "model-selection-limitations"


def test_inference_sensitivity_report_summarizes_filter_model_engine_and_bootstrap_changes(
    tmp_path: Path,
) -> None:
    baseline_exec = _fake_iqtree_tree_variant(
        tmp_path / "baseline-iqtree",
        tree_newick="((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);",
    )
    filtered_exec = _fake_iqtree_tree_variant(
        tmp_path / "filtered-iqtree",
        tree_newick="((A:0.1,C:0.1):0.25,(B:0.1,D:0.1):0.25);",
    )
    model_exec = _fake_iqtree_tree_variant(
        tmp_path / "model-iqtree",
        tree_newick="((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);",
    )
    fast_exec = _fake_fasttree(tmp_path / "fasttree-fixture")
    bootstrap_exec = _fake_iqtree(tmp_path / "bootstrap-iqtree")

    baseline = run_maximum_likelihood_tree_inference(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "baseline",
        model="GTR+G",
        executable=baseline_exec,
        prefix="baseline",
    )
    filtered = run_maximum_likelihood_tree_inference(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "filtered",
        model="GTR+G",
        executable=filtered_exec,
        prefix="filtered",
    )
    compared_model = run_maximum_likelihood_tree_inference(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "model",
        model="HKY+G",
        executable=model_exec,
        prefix="model",
    )
    fast = run_fast_tree_inference(
        fixture("alignments/example_alignment.fasta"),
        tmp_path / "fasttree.nwk",
        executable=fast_exec,
    )
    bootstrap = run_bootstrap_support_estimation(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "bootstrap",
        model="GTR+G",
        executable=bootstrap_exec,
        prefix="bootstrap",
    )

    report = build_inference_sensitivity_report(
        baseline.manifest_path,
        filtered_manifest_path=filtered.manifest_path,
        compare_model_manifest_path=compared_model.manifest_path,
        compare_engine_manifest_path=fast.manifest_path,
        bootstrap_manifest_path=bootstrap.manifest_path,
    )
    rendered = render_inference_sensitivity_report(
        baseline_manifest_path=baseline.manifest_path,
        filtered_manifest_path=filtered.manifest_path,
        compare_model_manifest_path=compared_model.manifest_path,
        compare_engine_manifest_path=fast.manifest_path,
        bootstrap_manifest_path=bootstrap.manifest_path,
        out_path=tmp_path / "inference-sensitivity.html",
    )

    assert report.alignment_filtering_sensitivity is not None
    assert report.model_sensitivity is not None
    assert report.engine_sensitivity is not None
    assert report.bootstrap_support is not None
    assert report.weak_backbone is not None
    assert rendered.output_path.exists()
    assert rendered.machine_manifest["has_bootstrap_support"] is True


def test_bootstrap_workflow_report_includes_support_and_backbone_sections(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "bootstrap-iqtree")
    workflow = run_bootstrap_support_estimation(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "bootstrap",
        model="GTR+G",
        executable=executable,
        prefix="bootstrap",
    )

    rendered = render_inference_workflow_report(
        manifest_path=workflow.manifest_path,
        out_path=tmp_path / "bootstrap-report.html",
    )

    assert "bootstrap-support-summary" in rendered.supplement_sections
    assert "bootstrap-support-histogram" in rendered.supplement_sections
    assert "low-support-branches" in rendered.supplement_sections
    assert "weak-backbone" in rendered.supplement_sections


def test_fasttree_workflow_report_includes_support_and_approximation_sections(
    tmp_path: Path,
) -> None:
    executable = _fake_fasttree(tmp_path / "fasttree-fixture")
    workflow = run_fast_tree_inference(
        fixture("alignments/example_alignment.fasta"),
        tmp_path / "fasttree.nwk",
        executable=executable,
    )

    rendered = render_inference_workflow_report(
        manifest_path=workflow.manifest_path,
        out_path=tmp_path / "fasttree-report.html",
    )

    assert "fasttree-approximation-limits" in rendered.supplement_sections
    assert "fasttree-support-summary" in rendered.supplement_sections
    assert "fasttree-support-histogram" in rendered.supplement_sections
    assert "fasttree-low-support-branches" in rendered.supplement_sections


def test_sh_alrt_workflow_report_includes_combined_support_sections(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "sh-alrt-iqtree")
    workflow = run_sh_alrt_support_estimation(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "sh-alrt",
        model="GTR+G",
        executable=executable,
        prefix="sh-alrt",
        sh_alrt_replicates=1000,
        bootstrap_replicates=1000,
    )

    rendered = render_inference_workflow_report(
        manifest_path=workflow.manifest_path,
        out_path=tmp_path / "sh-alrt-report.html",
    )

    assert "bootstrap-tree-set-validation" in rendered.supplement_sections
    assert "sh-alrt-support-summary" in rendered.supplement_sections
    assert "conflicting-support-branches" in rendered.supplement_sections
    assert "weak-backbone" in rendered.supplement_sections


def test_bootstrap_support_rejects_replicate_counts_below_iqtree_minimum(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "bootstrap-iqtree")

    with pytest.raises(
        EngineWorkflowError,
        match="ultrafast bootstrap requires at least 1000 replicates",
    ):
        run_bootstrap_support_estimation(
            fixture("alignments/example_alignment.fasta"),
            out_dir=tmp_path / "bootstrap",
            model="GTR+G",
            executable=executable,
            prefix="bootstrap",
            replicates=999,
        )


def test_alignment_trimming_workflow_report_includes_trimming_summary_section(
    tmp_path: Path,
) -> None:
    executable = _fake_trimal(tmp_path / "trimal-fixture")
    workflow = run_alignment_trimming(
        fixture("alignments/example_alignment_trim.fasta"),
        tmp_path / "trimmed.fasta",
        executable=executable,
        mode="strictplus",
    )

    rendered = render_inference_workflow_report(
        manifest_path=workflow.manifest_path,
        out_path=tmp_path / "trim-report.html",
    )

    assert "alignment-trimming-summary" in rendered.supplement_sections
    assert rendered.output_path.exists()


def test_bundle_inference_workflow_evidence_copies_inputs_outputs_and_manifests(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")
    model = run_model_selection(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "model",
        executable=executable,
        prefix="example",
    )
    ml = run_maximum_likelihood_tree_inference(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "ml",
        model="GTR+G",
        executable=executable,
        prefix="example",
    )

    result = bundle_inference_workflow_evidence(
        [model.manifest_path, ml.manifest_path],
        bundle_root=tmp_path / "inference-bundle",
    )

    assert result.bundle_root.exists()
    assert result.workflow_count == 2
    assert result.manifest_file_count == 2
    assert result.validation.mismatches == []
    assert (result.bundle_root / "outputs").exists()
