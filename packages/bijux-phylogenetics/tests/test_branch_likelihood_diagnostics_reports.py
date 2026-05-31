from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    evaluate_nucleotide_branch_likelihood_diagnostics_from_alignment,
    fit_local_clock_likelihood_from_alignment,
    fit_strict_clock_likelihood_from_alignment,
    write_branch_likelihood_diagnostic_table,
    write_local_clock_likelihood_artifacts,
    write_strict_clock_likelihood_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_write_branch_likelihood_diagnostic_table_writes_expected_header_and_flags(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "zero-branch-tree.nwk"
    tree_path.write_text(
        "(((A:0.0,B:0.5):0.4,C:0.6):0.3,D:0.8);",
        encoding="utf-8",
    )
    report = evaluate_nucleotide_branch_likelihood_diagnostics_from_alignment(
        tree_path,
        fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta"),
        model_name="jc69",
    )
    output_path = tmp_path / "branch_likelihood_diagnostics.tsv"

    write_branch_likelihood_diagnostic_table(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert output_path.read_text(encoding="utf-8").startswith(
        "model_name\tbranch_id\tchild_name\tdescendant_taxa\tbranch_length\tbaseline_log_likelihood\tcollapsed_branch_log_likelihood\tcontribution_proxy\twarning_flags\n"
    )
    row_by_descendants = {row["descendant_taxa"]: row for row in rows}
    assert "zero-branch-length" in row_by_descendants["A"]["warning_flags"].split("|")
    assert row_by_descendants["A"]["model_name"] == "JC69"


def test_write_strict_clock_likelihood_artifacts_materialize_branch_likelihood_diagnostics(
    tmp_path: Path,
) -> None:
    report = fit_strict_clock_likelihood_from_alignment(
        fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
        fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta"),
    )

    outputs = write_strict_clock_likelihood_artifacts(tmp_path, report)

    assert outputs["branch_likelihood_diagnostic_path"].name == (
        "branch_likelihood_diagnostics.tsv"
    )
    assert (
        outputs["branch_likelihood_diagnostic_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "model_name\tbranch_id\tchild_name\tdescendant_taxa\tbranch_length\tbaseline_log_likelihood\tcollapsed_branch_log_likelihood\tcontribution_proxy\twarning_flags\n"
        )
    )


def test_write_local_clock_likelihood_artifacts_materialize_branch_likelihood_diagnostics(
    tmp_path: Path,
) -> None:
    report = fit_local_clock_likelihood_from_alignment(
        fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
        fixture("alignments", "local_clock_likelihood_alignment_4_taxa.fasta"),
        fixture("metadata", "local_clock_regimes_4_taxa.tsv"),
    )

    outputs = write_local_clock_likelihood_artifacts(tmp_path, report)

    assert outputs["branch_likelihood_diagnostic_path"].name == (
        "branch_likelihood_diagnostics.tsv"
    )
    assert (
        outputs["branch_likelihood_diagnostic_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "model_name\tbranch_id\tchild_name\tdescendant_taxa\tbranch_length\tbaseline_log_likelihood\tcollapsed_branch_log_likelihood\tcontribution_proxy\twarning_flags\n"
        )
    )
