from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    fit_strict_clock_likelihood_from_alignment,
    write_strict_clock_branch_table,
    write_strict_clock_likelihood_artifacts,
    write_strict_clock_run_json,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_write_strict_clock_branch_table_writes_expected_rows(tmp_path: Path) -> None:
    report = fit_strict_clock_likelihood_from_alignment(
        fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
        fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta"),
    )
    output_path = tmp_path / "branch_rates.tsv"

    write_strict_clock_branch_table(output_path, report)

    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "branch_id\tchild_name\tdescendant_taxa\ttime_duration\toptimized_branch_length\toptimized_clock_rate",
        "root:clade:A|B|C|D/clade:A|B|C\t\tA|B|C\t1\t0.521738366943483\t0.521738366943483",
        "root:clade:A|B|C|D/clade:A|B|C/clade:A|B\t\tA|B\t1\t0.521738366943483\t0.521738366943483",
        "root:clade:A|B|C|D/clade:A|B|C/clade:A|B/taxon:A\tA\tA\t1\t0.521738366943483\t0.521738366943483",
        "root:clade:A|B|C|D/clade:A|B|C/clade:A|B/taxon:B\tB\tB\t1\t0.521738366943483\t0.521738366943483",
        "root:clade:A|B|C|D/clade:A|B|C/taxon:C\tC\tC\t2\t1.04347673388697\t0.521738366943483",
        "root:clade:A|B|C|D/taxon:D\tD\tD\t3\t1.56521510083045\t0.521738366943483",
    ]


def test_write_strict_clock_run_json_serializes_report_fields(tmp_path: Path) -> None:
    report = fit_strict_clock_likelihood_from_alignment(
        fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
        fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta"),
    )
    output_path = tmp_path / "run.json"

    write_strict_clock_run_json(output_path, report)

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["model_name"] == "JC69"
    assert payload["site_count"] == 12
    assert payload["pattern_count"] == 2
    assert payload["branch_count"] == 6
    assert payload["optimized_clock_rate"] == 0.5217383669434831
    assert payload["optimized_log_likelihood"] == -64.38812241070909
    assert len(payload["branch_rows"]) == 6


def test_write_strict_clock_likelihood_artifacts_materializes_governed_outputs(
    tmp_path: Path,
) -> None:
    report = fit_strict_clock_likelihood_from_alignment(
        fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
        fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta"),
    )

    outputs = write_strict_clock_likelihood_artifacts(tmp_path, report)

    assert sorted(outputs) == [
        "branch_likelihood_diagnostic_path",
        "branch_table_path",
        "run_json_path",
        "scaled_tree_path",
        "site_log_likelihood_path",
    ]
    assert outputs["scaled_tree_path"].read_text(encoding="utf-8").strip() == (
        "(((A:0.521738366943483,B:0.521738366943483):0.521738366943483,"
        "C:1.04347673388697):0.521738366943483,D:1.56521510083045);"
    )
    assert (
        outputs["branch_table_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "branch_id\tchild_name\tdescendant_taxa\ttime_duration\toptimized_branch_length\toptimized_clock_rate\n"
        )
    )
    assert (
        outputs["site_log_likelihood_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "model_name\ttaxon_order\tpattern_id\tpattern_weight\tsite_position\tsite_states\tlog_likelihood\n"
        )
    )
    assert (
        outputs["branch_likelihood_diagnostic_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "model_name\tbranch_id\tchild_name\tdescendant_taxa\tbranch_length\tbaseline_log_likelihood\tcollapsed_branch_log_likelihood\tcontribution_proxy\twarning_flags\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["model_name"] == "JC69"
    assert payload["optimized_clock_rate"] == 0.5217383669434831
