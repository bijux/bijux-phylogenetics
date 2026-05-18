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


def test_comparative_dtt_cli_writes_ledgers_and_svg(
    tmp_path: Path,
    capsys,
) -> None:
    summary_out = tmp_path / "dtt-summary.tsv"
    curve_out = tmp_path / "dtt-curve.tsv"
    clades_out = tmp_path / "dtt-clades.tsv"
    bins_out = tmp_path / "dtt-bins.tsv"
    excluded_out = tmp_path / "dtt-excluded.tsv"
    svg_out = tmp_path / "dtt.svg"

    exit_code = main(
        [
            "comparative",
            "dtt",
            str(fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")),
            str(
                fixture(
                    "example_traits_geiger_continuous_model_panel_twenty_four_taxa.tsv"
                )
            ),
            "--traits",
            "ou_truth,early_burst_truth",
            "--time-bin-count",
            "4",
            "--summary-out",
            str(summary_out),
            "--curve-out",
            str(curve_out),
            "--clades-out",
            str(clades_out),
            "--bins-out",
            str(bins_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--svg-out",
            str(svg_out),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["analyzed_taxon_count"] == 24
    assert payload["metrics"]["trait_column_count"] == 2
    assert payload["metrics"]["curve_point_count"] == 24
    assert payload["metrics"]["time_bin_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["rendered_point_count"] == 24
    assert payload["data"]["relative_scaling_applied"] is True
    assert summary_out.exists()
    assert curve_out.exists()
    assert clades_out.exists()
    assert bins_out.exists()
    assert excluded_out.exists()
    assert "<svg" in svg_out.read_text(encoding="utf-8")
