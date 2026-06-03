from __future__ import annotations

import pytest

from bijux_phylogenetics.phylo.likelihood.ctmc import (
    verify_ctmc_stationary_distribution,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_verify_ctmc_stationary_distribution_rejects_nonunique_fixture() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        verify_ctmc_stationary_distribution(
            [
                [-1.0, 1.0, 0.0],
                [1.0, -1.0, 0.0],
                [0.0, 0.0, 0.0],
            ],
            {"A": 0.5, "B": 0.5, "C": 0.0},
            state_labels=("A", "B", "C"),
        )

    assert error_info.value.code == "ctmc_stationary_distribution_not_unique"
