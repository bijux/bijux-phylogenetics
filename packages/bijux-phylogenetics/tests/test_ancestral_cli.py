from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.cli import main

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


def test_ancestral_continuous_cli_can_export_table(tmp_path: Path, capsys) -> None:
    table_path = tmp_path / "ancestral.tsv"
    summary_path = tmp_path / "ancestral-summary.tsv"
    uncertainty_path = tmp_path / "ancestral-uncertainty.tsv"
    exclusions_path = tmp_path / "ancestral-excluded.tsv"
    exit_code = main(
        [
            "ancestral",
            "continuous",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--table-out",
            str(table_path),
            "--summary-out",
            str(summary_path),
            "--uncertainty-out",
            str(uncertainty_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["model"] == "brownian"
    assert payload["metrics"]["internal_node_count"] == 3
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["unstable_node_count"] >= 0
    assert "estimate\tstandard_error" in table_path.read_text(encoding="utf-8")
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\talpha"
    )
    assert uncertainty_path.read_text(encoding="utf-8").startswith(
        "node\tnode_name\tdescendant_taxa\testimate\tstandard_error"
    )
    assert exclusions_path.read_text(encoding="utf-8") == "taxon\treason\n"


def test_ancestral_discrete_cli_reports_sparse_state_warning(capsys) -> None:
    exit_code = main(
        [
            "ancestral",
            "discrete",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_ancestral_sparse.tsv")),
            "--trait",
            "habitat",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert (
        "one or more discrete states are represented by fewer than two taxa and should be interpreted cautiously"
        in payload["warnings"]
    )
    assert (
        "one or more discrete ancestral nodes remain unstable across candidate states"
        in payload["warnings"]
    )


def test_ancestral_render_cli_writes_svg_with_internal_annotations(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "ancestral.svg"
    exit_code = main(
        [
            "ancestral",
            "render",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["rendered_internal_annotation_count"] == 3
    assert 'class="internal-annotation-label"' in output.read_text(encoding="utf-8")


def test_ancestral_report_cli_writes_html_and_svg(tmp_path: Path, capsys) -> None:
    output = tmp_path / "ancestral-report.html"
    exit_code = main(
        [
            "ancestral",
            "report",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--compare-model",
            "ou",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["report_kind"] == "ancestral-state"
    assert output.exists()
    assert output.with_suffix(".svg").exists()


def test_ancestral_sensitivity_cli_reports_available_comparisons(capsys) -> None:
    exit_code = main(
        [
            "ancestral",
            "sensitivity",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--compare-model",
            "ou",
            "--compare-tree",
            str(fixture("example_tree_topology_diff.nwk")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["has_model_sensitivity"] is True
    assert payload["metrics"]["has_tree_sensitivity"] is True


def test_ancestral_package_cli_writes_publication_bundle(
    tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "ancestral-package"
    exit_code = main(
        [
            "ancestral",
            "package",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--out-dir",
            str(output_dir),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["artifact_count"] == 7
    assert (output_dir / "figure-manifest.json").exists()


def test_ancestral_discrete_cli_rejects_ordered_fitch_model(capsys) -> None:
    try:
        main(
            [
                "ancestral",
                "report",
                str(fixture("example_tree.nwk")),
                str(fixture("example_traits_geography.tsv")),
                "--trait",
                "region",
                "--kind",
                "discrete",
                "--state-ordering",
                "ordered",
                "--out",
                "ignored.html",
            ]
        )
    except SystemExit as error:
        assert error.code == 2
    captured = capsys.readouterr()
    assert (
        "ordered ancestral discrete reconstruction requires a likelihood model"
        in captured.err
    )
