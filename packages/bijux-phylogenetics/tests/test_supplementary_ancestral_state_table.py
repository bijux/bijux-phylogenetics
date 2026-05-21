from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.reports import write_supplementary_ancestral_state_table

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


def test_write_supplementary_ancestral_state_table_writes_continuous_internal_nodes(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-ancestral-continuous.tsv"

    result = write_supplementary_ancestral_state_table(
        output_path,
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_comparative.tsv"),
        trait="response",
        reconstruction_kind="continuous",
        model="brownian",
    )

    assert result.output_path == output_path
    assert result.reconstruction_kind == "continuous"
    assert result.model == "brownian"
    assert result.row_count == 3
    assert result.analysis_taxon_count == 4
    assert result.excluded_taxon_count == 0
    assert result.unstable_node_count >= 0
    assert all(row.reconstruction_kind == "continuous" for row in result.rows)
    assert all(row.estimate_value is not None for row in result.rows)
    assert all(row.standard_error is not None for row in result.rows)
    assert all(row.most_likely_state is None for row in result.rows)
    root_rows = [row for row in result.rows if row.descendant_taxon_count == 4]
    assert len(root_rows) == 1
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3
    assert rows[0]["reconstruction_kind"] == "continuous"
    assert rows[0]["model"] == "brownian"
    assert rows[0]["state_probabilities"] == "{}"


def test_write_supplementary_ancestral_state_table_writes_discrete_uncertainty_and_warnings(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-ancestral-discrete.tsv"

    result = write_supplementary_ancestral_state_table(
        output_path,
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_ancestral_sparse.tsv"),
        trait="habitat",
        reconstruction_kind="discrete",
        model="equal-rates",
    )

    assert result.output_path == output_path
    assert result.reconstruction_kind == "discrete"
    assert result.model == "equal-rates"
    assert result.row_count == 3
    assert result.analysis_taxon_count == 4
    assert result.excluded_taxon_count == 0
    assert result.unstable_node_count >= 0
    assert all(row.reconstruction_kind == "discrete" for row in result.rows)
    assert all(row.most_likely_state is not None for row in result.rows)
    assert all(row.estimate_value is None for row in result.rows)
    assert any(row.warning_count > 0 for row in result.rows)
    assert any(
        "one or more discrete states are represented by fewer than two taxa" in warning
        for row in result.rows
        for warning in row.warnings
    )
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3
    first_probabilities = json.loads(rows[0]["state_probabilities"])
    assert set(first_probabilities) == {"desert", "forest"}
    assert rows[0]["root_prior_mode"] == "equal"
    assert rows[0]["estimate_value"] == ""
