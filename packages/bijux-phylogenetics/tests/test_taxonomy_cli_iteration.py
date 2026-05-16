from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_cli_taxonomy_synonyms_json_output(capsys) -> None:
    exit_code = main(
        [
            "taxonomy",
            "synonyms",
            str(fixture("example_taxonomy_tree.nwk")),
            "--synonym-table",
            str(fixture("example_taxon_synonyms.tsv")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["candidate_count"] == 1
    assert payload["metrics"]["ambiguous_mapping_count"] == 1


def test_cli_taxonomy_resolve_synonyms_writes_outputs(tmp_path: Path, capsys) -> None:
    out_path = tmp_path / "resolved.nwk"
    mapping_path = tmp_path / "resolved.tsv"
    exit_code = main(
        [
            "taxonomy",
            "resolve-synonyms",
            str(fixture("example_taxonomy_tree.nwk")),
            "--synonym-table",
            str(fixture("example_taxon_synonyms.tsv")),
            "--out",
            str(out_path),
            "--mapping-out",
            str(mapping_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert out_path.exists()
    assert mapping_path.exists()
    assert payload["metrics"]["renamed_taxa"] == 1


def test_cli_taxonomy_stability_json_output(capsys) -> None:
    exit_code = main(
        [
            "taxonomy",
            "stability",
            "--run",
            f"tree={fixture('example_taxon_workflow_tree.nwk')}",
            "--run",
            f"alignment={fixture('example_taxon_workflow_alignment.fasta')}",
            "--run",
            f"filtered={fixture('example_taxon_workflow_filtered_alignment.fasta')}",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["source_count"] == 3
    assert payload["metrics"]["unstable_taxa"] == 2


def test_cli_taxonomy_rank_consistency_json_output(capsys) -> None:
    exit_code = main(
        [
            "taxonomy",
            "rank-consistency",
            str(fixture("example_taxonomy_rank_mixed.nwk")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["mixed_ranks"] is True
    assert payload["metrics"]["rank_count"] >= 4


def test_cli_taxonomy_accepted_names_writes_mapping(tmp_path: Path, capsys) -> None:
    out_path = tmp_path / "accepted.tsv"
    exit_code = main(
        [
            "taxonomy",
            "accepted-names",
            str(fixture("example_taxonomy_tree.nwk")),
            "--synonym-table",
            str(fixture("example_taxon_synonyms.tsv")),
            "--out",
            str(out_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert out_path.exists()
    assert payload["metrics"]["resolved_count"] == 1
    assert payload["metrics"]["ambiguous_count"] == 1


def test_cli_taxonomy_audit_json_output(capsys) -> None:
    exit_code = main(
        [
            "taxonomy",
            "audit",
            str(fixture("example_taxonomy_rank_mixed.nwk")),
            "--synonym-table",
            str(fixture("example_taxon_synonyms.tsv")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["status"] == "needs_review"
    assert payload["metrics"]["mapping_conflict_count"] >= 1


def test_cli_report_taxonomy_json_output_includes_workflow_metrics(
    tmp_path: Path, capsys
) -> None:
    out_path = tmp_path / "taxonomy-report.html"
    exit_code = main(
        [
            "report",
            "taxonomy",
            "--tree",
            str(fixture("example_taxon_workflow_tree.nwk")),
            "--metadata",
            str(fixture("example_taxon_workflow_metadata.csv")),
            "--traits",
            str(fixture("example_taxon_workflow_traits.csv")),
            "--alignment",
            str(fixture("example_taxon_workflow_alignment.fasta")),
            "--filtered-alignment",
            str(fixture("example_taxon_workflow_filtered_alignment.fasta")),
            "--inference-tree",
            str(fixture("example_taxon_workflow_inference.nwk")),
            "--reported-taxa",
            str(fixture("example_taxon_workflow_reported.csv")),
            "--out",
            str(out_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert out_path.exists()
    assert payload["metrics"]["crosswalk_rows"] == 4
    assert payload["metrics"]["excluded_taxa"] == 2
    assert payload["metrics"]["loss_stage_count"] == 3
    assert payload["metrics"]["unstable_taxa"] == 3
