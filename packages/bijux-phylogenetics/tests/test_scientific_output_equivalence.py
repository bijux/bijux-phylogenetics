from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.validation import compare_scientific_output

FIXTURES = Path(__file__).parent / "fixtures"


def test_tree_support_and_branch_length_formatting_can_drift_without_failing(
    tmp_path: Path,
) -> None:
    expected = FIXTURES / "trees" / "example_tree_support_iqtree_composite.nwk"
    observed = tmp_path / "observed.nwk"
    observed.write_text(
        "((A:0.100000000000,B:0.100000000000)82/97:0.200000000000,(C:0.100000000000,D:0.100000000000)79/96:0.200000000000)91/99;",
        encoding="utf-8",
    )

    report = compare_scientific_output(expected, observed)

    assert report.equivalent is True
    assert report.compared_file_count == 1


def test_tree_support_change_reports_scientific_mismatch() -> None:
    expected = FIXTURES / "trees" / "example_tree_support_conflict_left.nwk"
    observed = FIXTURES / "trees" / "example_tree_support_conflict_right.nwk"

    report = compare_scientific_output(expected, observed)

    assert report.equivalent is False
    assert any(issue.kind == "support_value_mismatch" for issue in report.issues)


def test_tabular_numeric_tolerance_accepts_harmless_rounding(tmp_path: Path) -> None:
    expected = tmp_path / "expected.tsv"
    observed = tmp_path / "observed.tsv"
    expected.write_text("name\tvalue\nalpha\t0.1\nbeta\t2\n", encoding="utf-8")
    observed.write_text(
        "name\tvalue\nalpha\t0.1000000005\nbeta\t2.0000000001\n",
        encoding="utf-8",
    )

    report = compare_scientific_output(expected, observed)

    assert report.equivalent is True


def test_tabular_identity_fields_match_rows_by_stable_node_id(tmp_path: Path) -> None:
    expected = tmp_path / "expected.tsv"
    observed = tmp_path / "observed.tsv"
    expected.write_text(
        "node\tvalue\nA|B\t1\nC|D\t2\n",
        encoding="utf-8",
    )
    observed.write_text(
        "node\tvalue\nC|D\t2\nA|B\t1\n",
        encoding="utf-8",
    )

    report = compare_scientific_output(expected, observed)

    assert report.equivalent is True


def test_probability_tables_allow_governed_confidence_tolerance(tmp_path: Path) -> None:
    expected = tmp_path / "host-state-nodes.tsv"
    observed = tmp_path / "observed.tsv"
    expected.write_text(
        (
            "node\thost_probabilities\tparent_confidence\n"
            'A|B\t"{""bat"": 0.01, ""canid"": 0.79, ""livestock"": 0.20}"\t0.79\n'
        ),
        encoding="utf-8",
    )
    observed.write_text(
        (
            "node\thost_probabilities\tparent_confidence\n"
            'A|B\t"{""bat"": 0.03, ""canid"": 0.77, ""livestock"": 0.20}"\t0.77\n'
        ),
        encoding="utf-8",
    )

    report = compare_scientific_output(expected, observed)

    assert report.equivalent is True


def test_html_report_contract_ignores_whitespace_but_keeps_links_and_manifest(
    tmp_path: Path,
) -> None:
    manifest = {
        "report_kind": "tree",
        "sections": ["reviewer-summary", "tree-validation"],
    }
    expected = tmp_path / "expected.html"
    observed = tmp_path / "observed.html"
    expected.write_text(
        """
<!doctype html>
<html>
  <head><title>Bijux Tree Report</title></head>
  <body>
    <h1>Bijux Tree Report</h1>
    <h2>reviewer-summary</h2>
    <h2>tree-validation</h2>
    <a href="workflow-summary.tsv">workflow summary</a>
    <script id="bijux-report-manifest" type="application/json">{"report_kind":"tree","sections":["reviewer-summary","tree-validation"]}</script>
  </body>
</html>
""".strip(),
        encoding="utf-8",
    )
    observed.write_text(
        f"""
<html><head>
<title>Bijux Tree Report</title>
</head><body>
<h1>Bijux Tree Report</h1>
<h2>
reviewer-summary
</h2>
<h2>tree-validation</h2>
<script id="bijux-report-manifest" type="application/json">
{json.dumps(manifest, indent=2)}
</script>
<a href="workflow-summary.tsv">workflow summary</a>
</body></html>
""".strip(),
        encoding="utf-8",
    )

    report = compare_scientific_output(expected, observed)

    assert report.equivalent is True


def test_fasta_sequence_equivalence_ignores_wrapping(tmp_path: Path) -> None:
    expected = tmp_path / "expected.fasta"
    observed = tmp_path / "observed.fasta"
    expected.write_text(">alpha\nACGTACGT\n>beta\nTTAA\n", encoding="utf-8")
    observed.write_text(">alpha\nACGT\nACGT\n>beta\nTTAA\n", encoding="utf-8")

    report = compare_scientific_output(expected, observed)

    assert report.equivalent is True
