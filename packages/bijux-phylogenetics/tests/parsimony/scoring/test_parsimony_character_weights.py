from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.parsimony import (
    reconstruct_acctran,
    reconstruct_deltran,
    resolve_parsimony_character_weights,
    score_camin_sokal,
    score_dollo,
    score_fitch,
    score_sankoff,
    score_wagner,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


@pytest.mark.parametrize(
    ("factory", "args", "weights_name", "expected_total"),
    [
        (
            score_fitch,
            (fixture("fitch_tree.nwk"), fixture("fitch_binary_matrix.tsv")),
            "fitch_character_weights.tsv",
            3.5,
        ),
        (
            score_wagner,
            (fixture("fitch_tree.nwk"), fixture("wagner_ordinal_matrix.tsv")),
            "wagner_character_weights.tsv",
            7.5,
        ),
        (
            score_sankoff,
            (
                fixture("sankoff_tree_5_taxa.nwk"),
                fixture("sankoff_character_matrix.tsv"),
                fixture("sankoff_cost_matrix.tsv"),
            ),
            "sankoff_character_weights.tsv",
            7.5,
        ),
        (
            score_dollo,
            (fixture("dollo_tree_5_taxa.nwk"), fixture("dollo_binary_matrix.tsv")),
            "dollo_character_weights.tsv",
            9.5,
        ),
        (
            score_camin_sokal,
            (
                fixture("camin_sokal_tree_5_taxa.nwk"),
                fixture("camin_sokal_binary_matrix.tsv"),
            ),
            "camin_sokal_character_weights.tsv",
            5.5,
        ),
        (
            reconstruct_acctran,
            (
                fixture("acctran_tree_5_taxa.nwk"),
                fixture("acctran_ambiguous_matrix.tsv"),
            ),
            "acctran_character_weights.tsv",
            5.0,
        ),
        (
            reconstruct_deltran,
            (
                fixture("acctran_tree_5_taxa.nwk"),
                fixture("acctran_ambiguous_matrix.tsv"),
            ),
            "acctran_character_weights.tsv",
            5.0,
        ),
    ],
)
def test_native_parsimony_surfaces_honor_explicit_character_weights(
    factory,
    args,
    weights_name: str,
    expected_total: float,
) -> None:
    report = factory(*args, character_weights=fixture(weights_name))

    assert report.weights_path == fixture(weights_name)
    assert report.total_weighted_score == expected_total
    assert all(row.character_weight >= 0.0 for row in report.step_rows)
    assert sum(row.weighted_score for row in report.step_rows) == expected_total


def test_resolve_parsimony_character_weights_rejects_negative_values() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        resolve_parsimony_character_weights(
            ["char01_split", "char02_left_conflict"],
            fixture("parsimony_negative_character_weight.tsv"),
        )

    assert error_info.value.code == "parsimony_character_weight_negative_value"
    assert error_info.value.details["character_id"] == "char01_split"
