from __future__ import annotations

from pathlib import Path
import random

import pytest

from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.ancestral.discrete.review import (
    summarize_irreversible_discrete_reconstruction,
    summarize_irreversible_discrete_report,
    write_irreversible_discrete_fit_table,
    write_irreversible_discrete_node_table,
    write_irreversible_discrete_summary_table,
    write_irreversible_discrete_transition_table,
)
from bijux_phylogenetics.io.trees import load_tree

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


def test_discrete_reconstruction_rejects_asymmetric_constraints_under_symmetric_model() -> (
    None
):
    with pytest.raises(ValueError):
        reconstruct_discrete_ancestral_states(
            fixture("example_tree.nwk"),
            fixture("example_traits_geography.tsv"),
            trait="region",
            model="symmetric",
            allowed_transition_pairs=[("north", "south")],
        )


def test_summarize_irreversible_discrete_reconstruction_prefers_constraint_on_simulated_losses(
    tmp_path: Path,
) -> None:
    traits_path = _write_simulated_irreversible_traits(
        tmp_path / "simulated-irreversible.tsv"
    )
    report = summarize_irreversible_discrete_reconstruction(
        fixture("example_tree_six_taxa.nwk"),
        traits_path,
        trait="state",
        model="all-rates-different",
        allowed_transition_pairs=[("present", "absent")],
    )
    summary = summarize_irreversible_discrete_report(report)

    assert summary.constrained_parameter_count == 1
    assert summary.unconstrained_parameter_count == 2
    assert summary.forbidden_transition_count == 1
    assert summary.preferred_constraint == "constrained"
    forbidden_row = next(
        row
        for row in report.transition_rows
        if row.source_state == "absent" and row.target_state == "present"
    )
    assert forbidden_row.constrained_transition_allowed is False
    assert forbidden_row.unconstrained_transition_allowed is True
    assert forbidden_row.constrained_rate == 0.0


def test_write_irreversible_discrete_tables_emit_expected_ledgers(
    tmp_path: Path,
) -> None:
    traits_path = _write_simulated_irreversible_traits(
        tmp_path / "simulated-irreversible.tsv"
    )
    report = summarize_irreversible_discrete_reconstruction(
        fixture("example_tree_six_taxa.nwk"),
        traits_path,
        trait="state",
        model="all-rates-different",
        allowed_transition_pairs=[("present", "absent")],
    )
    summary_path = tmp_path / "irreversible-summary.tsv"
    fit_path = tmp_path / "irreversible-fits.tsv"
    node_path = tmp_path / "irreversible-nodes.tsv"
    transition_path = tmp_path / "irreversible-transitions.tsv"

    write_irreversible_discrete_summary_table(summary_path, report)
    write_irreversible_discrete_fit_table(fit_path, report)
    write_irreversible_discrete_node_table(node_path, report)
    write_irreversible_discrete_transition_table(transition_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    fit_rows = fit_path.read_text(encoding="utf-8").splitlines()
    node_rows = node_path.read_text(encoding="utf-8").splitlines()
    transition_rows = transition_path.read_text(encoding="utf-8").splitlines()

    assert summary_rows[0].startswith(
        "trait\ttaxon_column\tmodel\tanalyzed_taxon_count\tconstrained_log_likelihood"
    )
    assert fit_rows[0].startswith(
        "constraint_mode\tmodel\tanalyzed_taxon_count\tlog_likelihood"
    )
    assert node_rows[0].startswith(
        "node\tdescendant_taxa\tconstrained_state\tunconstrained_state"
    )
    assert transition_rows[0].startswith(
        "source_state\ttarget_state\tconstrained_transition_allowed"
    )
    assert len(fit_rows) == 3


def _write_simulated_irreversible_traits(path: Path) -> Path:
    tree = load_tree(fixture("example_tree_six_taxa.nwk"))
    for seed in range(1, 50):
        rng = random.Random(seed)  # nosec B311
        states = _simulate_irreversible_node_states(tree.root, rng)
        tip_states = dict(states)
        if len(set(tip_states.values())) == 2:
            path.write_text(
                "taxon\tstate\n"
                + "".join(
                    f"{taxon}\t{tip_states[taxon]}\n" for taxon in sorted(tip_states)
                ),
                encoding="utf-8",
            )
            return path
    raise AssertionError("failed to simulate both irreversible states")


def _simulate_irreversible_node_states(
    node, rng: random.Random, state: str = "present"
):
    states: dict[str, str] = {}
    current_state = state
    if node.is_leaf():
        states[node.name] = current_state
        return states
    for child in node.children:
        next_state = current_state
        if current_state == "present":
            branch_length = child.branch_length or 0.0
            if rng.random() < min(0.45 * branch_length, 0.95):
                next_state = "absent"
        states.update(_simulate_irreversible_node_states(child, rng, next_state))
    return states
