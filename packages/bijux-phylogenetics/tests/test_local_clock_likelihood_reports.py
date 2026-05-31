from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    fit_local_clock_likelihood_from_alignment,
    write_local_clock_branch_table,
    write_local_clock_likelihood_artifacts,
    write_local_clock_regime_table,
    write_local_clock_run_json,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def build_report():
    return fit_local_clock_likelihood_from_alignment(
        fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
        fixture("alignments", "local_clock_likelihood_alignment_4_taxa.fasta"),
        fixture("metadata", "local_clock_regimes_4_taxa.tsv"),
    )


def test_write_local_clock_branch_table_writes_expected_assignments(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "branch_rates.tsv"

    write_local_clock_branch_table(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 6
    rows_by_taxa = {row["descendant_taxa"]: row for row in rows}
    assert rows_by_taxa["A|B|C"]["regime_id"] == "abc_stem"
    assert rows_by_taxa["A|B"]["regime_id"] == "ab_clade"
    assert rows_by_taxa["C"]["regime_id"] == "background"
    assert rows_by_taxa["D"]["target_kind"] == "background"


def test_write_local_clock_regime_table_writes_expected_rows(tmp_path: Path) -> None:
    report = build_report()
    output_path = tmp_path / "regimes.tsv"

    write_local_clock_regime_table(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 3
    rows_by_id = {row["regime_id"]: row for row in rows}
    assert rows_by_id["background"]["branch_count"] == "2"
    assert rows_by_id["ab_clade"]["target_kind"] == "clade"
    assert rows_by_id["abc_stem"]["target_kind"] == "branch"
    assert float(rows_by_id["abc_stem"]["optimized_clock_rate"]) > float(
        rows_by_id["ab_clade"]["optimized_clock_rate"]
    )


def test_write_local_clock_run_json_serializes_report_fields(tmp_path: Path) -> None:
    report = build_report()
    output_path = tmp_path / "run.json"

    write_local_clock_run_json(output_path, report)

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["model_name"] == "JC69 local-clock"
    assert payload["regime_count"] == 3
    assert payload["preferred_model_by_aic"] == "local-clock"
    assert payload["strict_clock_aic"] > payload["aic"]
    assert len(payload["branch_rows"]) == 6
    assert len(payload["regime_rows"]) == 3


def test_write_local_clock_likelihood_artifacts_materializes_governed_outputs(
    tmp_path: Path,
) -> None:
    report = build_report()

    outputs = write_local_clock_likelihood_artifacts(tmp_path, report)

    assert sorted(outputs) == [
        "branch_likelihood_diagnostic_path",
        "branch_table_path",
        "regime_table_path",
        "run_json_path",
        "scaled_tree_path",
        "site_log_likelihood_path",
    ]
    assert outputs["scaled_tree_path"].read_text(encoding="utf-8").strip() == (
        report.scaled_tree_newick
    )
    assert (
        outputs["branch_table_path"]
        .read_text(encoding="utf-8")
        .startswith("branch_id\tchild_name\tdescendant_taxa\tregime_id\t")
    )
    assert (
        outputs["regime_table_path"]
        .read_text(encoding="utf-8")
        .startswith("regime_id\ttarget_kind\ttarget_label\tdescendant_taxa\t")
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
    assert payload["preferred_model_by_aic"] == "local-clock"
