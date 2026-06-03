from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.ancestral.common import load_discrete_dataset
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


def test_load_discrete_dataset_rejects_explicit_missing_state_tokens(
    tmp_path: Path,
) -> None:
    traits_path = tmp_path / "missing_states.tsv"
    traits_path.write_text(
        "taxon\tstate\nA\t?\nB\t2\n",
        encoding="utf-8",
    )

    with pytest.raises(
        AncestralReconstructionError,
        match="uses missing token '\\?'",
    ):
        load_discrete_dataset(
            fixture("felsenstein_two_tip_tree.nwk"),
            traits_path,
            trait="state",
            taxon_column="taxon",
        )


def test_load_discrete_dataset_rejects_ambiguous_state_tokens(
    tmp_path: Path,
) -> None:
    traits_path = tmp_path / "ambiguous_states.tsv"
    traits_path.write_text(
        "taxon\tstate\nA\t0|1\nB\t2\n",
        encoding="utf-8",
    )

    with pytest.raises(
        AncestralReconstructionError,
        match="uses ambiguous token '0\\|1'",
    ):
        load_discrete_dataset(
            fixture("felsenstein_two_tip_tree.nwk"),
            traits_path,
            trait="state",
            taxon_column="taxon",
        )
