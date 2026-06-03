from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralEstimate,
    DiscreteTransitionRateRow,
)
from bijux_phylogenetics.biogeography import (
    summarize_constrained_geographic_model,
    summarize_constrained_geographic_report,
    write_constrained_geographic_exclusion_table,
    write_constrained_geographic_fit_table,
    write_constrained_geographic_summary_table,
    write_constrained_geographic_transition_table,
    write_unsupported_geographic_transition_claim_table,
)
import bijux_phylogenetics.biogeography.state_models.constrained as constrained_geography_module
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

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


@pytest.mark.slow
def test_summarize_constrained_geographic_model_supports_aliases() -> None:
    for model in ("er", "sym", "ard"):
        report = summarize_constrained_geographic_model(
            fixture("example_tree.nwk"),
            fixture("example_traits_geography.tsv"),
            fixture("example_geographic_adjacency.tsv"),
            trait="region",
            model=model,
        )

        assert report.model == model
        assert report.fit_rows
        assert report.transition_rows


def test_summarize_constrained_geographic_model_reports_unsupported_transition_claims(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _FakeNode:
        def __init__(
            self, name: str | None, children: list[_FakeNode] | None = None
        ) -> None:
            self.name = name
            self.children = children or []

        def is_leaf(self) -> bool:
            return not self.children

    def _estimate(
        node: str,
        *,
        descendant_taxa: list[str],
        most_likely_state: str,
        state_probabilities: dict[str, float],
        is_tip: bool = False,
        node_name: str | None = None,
    ) -> DiscreteAncestralEstimate:
        return DiscreteAncestralEstimate(
            node=node,
            node_name=node_name,
            is_tip=is_tip,
            descendant_taxa=descendant_taxa,
            state_set=sorted(state_probabilities),
            most_likely_state=most_likely_state,
            state_probabilities=state_probabilities,
            ambiguous=sum(
                probability > 0.0 for probability in state_probabilities.values()
            )
            > 1,
            confidence=max(state_probabilities.values()),
            interpretation="fixture",
            unstable=False,
            downstream_risks=[],
        )

    def _report(
        *,
        estimates: list[DiscreteAncestralEstimate],
        transition_rate_rows: list[DiscreteTransitionRateRow],
        warnings: list[str],
    ) -> SimpleNamespace:
        return SimpleNamespace(
            estimates=estimates,
            transition_rate_rows=transition_rate_rows,
            warnings=warnings,
            model="all-rates-different",
            taxon_count=4,
            log_likelihood=-4.2,
            parameter_count=6,
            aic=20.4,
        )

    ab = _FakeNode(None, [_FakeNode("A"), _FakeNode("B")])
    cd = _FakeNode(None, [_FakeNode("C"), _FakeNode("D")])
    root = _FakeNode(None, [ab, cd])

    unconstrained_report = _report(
        estimates=[
            _estimate(
                "A|B|C|D",
                descendant_taxa=["A", "B", "C", "D"],
                most_likely_state="south",
                state_probabilities={
                    "island": 0.49,
                    "north": 0.0,
                    "south": 0.51,
                },
            ),
            _estimate(
                "A|B",
                descendant_taxa=["A", "B"],
                most_likely_state="north",
                state_probabilities={
                    "island": 0.0,
                    "north": 1.0,
                    "south": 0.0,
                },
            ),
            _estimate(
                "C|D",
                descendant_taxa=["C", "D"],
                most_likely_state="south",
                state_probabilities={
                    "island": 0.0,
                    "north": 0.0,
                    "south": 1.0,
                },
            ),
            _estimate(
                "A",
                descendant_taxa=["A"],
                most_likely_state="north",
                state_probabilities={"north": 1.0},
                is_tip=True,
                node_name="A",
            ),
            _estimate(
                "B",
                descendant_taxa=["B"],
                most_likely_state="north",
                state_probabilities={"north": 1.0},
                is_tip=True,
                node_name="B",
            ),
            _estimate(
                "C",
                descendant_taxa=["C"],
                most_likely_state="island",
                state_probabilities={"island": 1.0},
                is_tip=True,
                node_name="C",
            ),
            _estimate(
                "D",
                descendant_taxa=["D"],
                most_likely_state="south",
                state_probabilities={"south": 1.0},
                is_tip=True,
                node_name="D",
            ),
        ],
        transition_rate_rows=[
            DiscreteTransitionRateRow(
                source_state="south",
                target_state="north",
                transition_allowed=True,
                step_distance=1,
                rate=0.51,
            ),
            DiscreteTransitionRateRow(
                source_state="island",
                target_state="north",
                transition_allowed=False,
                step_distance=1,
                rate=0.49,
            ),
        ],
        warnings=[
            "one or more discrete ancestral nodes remain unstable across candidate states",
            "low-confidence ancestral state assignments should not be overinterpreted as definitive transitions",
        ],
    )
    constrained_report = _report(
        estimates=[
            _estimate(
                "A|B|C|D",
                descendant_taxa=["A", "B", "C", "D"],
                most_likely_state="south",
                state_probabilities={"south": 1.0},
            ),
            _estimate(
                "A|B",
                descendant_taxa=["A", "B"],
                most_likely_state="north",
                state_probabilities={"north": 1.0},
            ),
            _estimate(
                "C|D",
                descendant_taxa=["C", "D"],
                most_likely_state="south",
                state_probabilities={"south": 1.0},
            ),
            _estimate(
                "A",
                descendant_taxa=["A"],
                most_likely_state="north",
                state_probabilities={"north": 1.0},
                is_tip=True,
                node_name="A",
            ),
            _estimate(
                "B",
                descendant_taxa=["B"],
                most_likely_state="north",
                state_probabilities={"north": 1.0},
                is_tip=True,
                node_name="B",
            ),
            _estimate(
                "C",
                descendant_taxa=["C"],
                most_likely_state="island",
                state_probabilities={"island": 1.0},
                is_tip=True,
                node_name="C",
            ),
            _estimate(
                "D",
                descendant_taxa=["D"],
                most_likely_state="south",
                state_probabilities={"south": 1.0},
                is_tip=True,
                node_name="D",
            ),
        ],
        transition_rate_rows=[
            DiscreteTransitionRateRow(
                source_state="south",
                target_state="north",
                transition_allowed=True,
                step_distance=1,
                rate=0.51,
            ),
            DiscreteTransitionRateRow(
                source_state="island",
                target_state="north",
                transition_allowed=False,
                step_distance=1,
                rate=0.0,
            ),
        ],
        warnings=[],
    )

    monkeypatch.setattr(
        constrained_geography_module,
        "audit_discrete_state_coding",
        lambda *args, **kwargs: SimpleNamespace(rows=[]),
    )
    monkeypatch.setattr(
        constrained_geography_module,
        "load_discrete_dataset",
        lambda *args, **kwargs: SimpleNamespace(
            tree=SimpleNamespace(root=root),
            observed_states=["north", "south", "island"],
            taxon_column="taxon",
        ),
    )
    monkeypatch.setattr(
        constrained_geography_module,
        "reconstruct_discrete_ancestral_states",
        lambda *args, allowed_transition_pairs=None, **kwargs: (
            constrained_report
            if allowed_transition_pairs is not None
            else unconstrained_report
        ),
    )

    traits_path = tmp_path / "simulated-geography.tsv"
    traits_path.write_text(
        "taxon\tregion\nA\tnorth\nB\tnorth\nC\tisland\nD\tsouth\n",
        encoding="utf-8",
    )

    report = summarize_constrained_geographic_model(
        fixture("example_tree.nwk"),
        traits_path,
        fixture("example_geographic_adjacency.tsv"),
        trait="region",
        model="ard",
    )
    summary = summarize_constrained_geographic_report(report)

    assert summary.unsupported_transition_claim_count == len(
        report.unsupported_claim_rows
    )
    assert report.unsupported_claim_rows
    assert any(
        {
            row.unconstrained_source_region,
            row.unconstrained_target_region,
        }
        == {"north", "island"}
        for row in report.unsupported_claim_rows
    )


def test_summarize_constrained_geographic_model_rejects_invalid_adjacency_matrix(
    tmp_path: Path,
) -> None:
    adjacency_path = tmp_path / "bad-adjacency.tsv"
    adjacency_path.write_text(
        "region\tnorth\tsouth\nnorth\t0\t1\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="missing rows"):
        summarize_constrained_geographic_model(
            fixture("example_tree.nwk"),
            fixture("example_traits_geography.tsv"),
            adjacency_path,
            trait="region",
            model="ard",
        )


@pytest.mark.slow
def test_write_constrained_geographic_tables_emit_expected_ledgers(
    tmp_path: Path,
) -> None:
    report = summarize_constrained_geographic_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        fixture("example_geographic_adjacency.tsv"),
        trait="region",
        model="ard",
    )

    summary_path = write_constrained_geographic_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    fits_path = write_constrained_geographic_fit_table(
        tmp_path / "fits.tsv",
        report,
    )
    transitions_path = write_constrained_geographic_transition_table(
        tmp_path / "transitions.tsv",
        report,
    )
    unsupported_path = write_unsupported_geographic_transition_claim_table(
        tmp_path / "unsupported.tsv",
        report,
    )
    exclusions_path = write_constrained_geographic_exclusion_table(
        tmp_path / "exclusions.tsv",
        report,
    )

    assert "preferred_constraint" in summary_path.read_text(encoding="utf-8")
    assert "constraint_mode" in fits_path.read_text(encoding="utf-8")
    assert "transition_allowed" in transitions_path.read_text(encoding="utf-8")
    assert "claim_resolved" in unsupported_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")


def test_summarize_constrained_geographic_model_rejects_in_tree_unknown_regions(
    tmp_path: Path,
) -> None:
    traits_path = tmp_path / "bad-geography.tsv"
    traits_path.write_text(
        "taxon\tregion\nA\tnorth\nB\tlagoon\nC\tsouth\nD\tisland\n",
        encoding="utf-8",
    )

    with pytest.raises(
        AncestralReconstructionError,
        match="does not define one or more analyzed region labels",
    ):
        summarize_constrained_geographic_model(
            fixture("example_tree.nwk"),
            traits_path,
            fixture("example_geographic_adjacency.tsv"),
            trait="region",
            model="ard",
        )
