from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FixedTopologyDnaRunReport,
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
    run_fixed_topology_dna_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    build_fixed_topology_dna_run_manifest,
    list_metropolis_hastings_retained_sample_ids,
    load_bayesian_run_manifest,
    replay_fixed_topology_dna_run_manifest,
    write_bayesian_run_manifest,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_fixed_topology_dna_run_manifest_records_reproducibility_fields(
    tmp_path: Path,
) -> None:
    run_report = _build_fixed_topology_dna_run_report(
        tree_path=fixture("trees", "hky85_likelihood_tree_2_taxa.nwk"),
        alignment_path=fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta"),
        seed=3,
    )
    sample_id_path = _write_sample_id_output(
        tmp_path / "retained-sample-ids.json",
        run_report,
    )

    manifest = build_fixed_topology_dna_run_manifest(
        run_report=run_report,
        tree_path=fixture("trees", "hky85_likelihood_tree_2_taxa.nwk"),
        alignment_path=fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta"),
        output_paths=[sample_id_path],
    )

    assert manifest.manifest_kind == "bayesian-run"
    assert manifest.run_kind == "fixed-topology-dna"
    assert manifest.model_name == "HKY85"
    assert manifest.seed == 3
    assert manifest.chain_count == 1
    assert manifest.retained_sample_count == len(run_report.chain_report.sampled_states)
    assert manifest.burnin_policy.policy_name == "none"
    assert manifest.execution_configuration == {
        "iteration_count": 9,
        "sample_every": 1,
        "observation_policy": "reject",
    }
    assert manifest.input_checksums[
        str(fixture("trees", "hky85_likelihood_tree_2_taxa.nwk"))
    ]
    assert manifest.input_checksums[
        str(fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta"))
    ]
    assert manifest.output_checksums[str(sample_id_path)]
    assert manifest.retained_sample_ids == list_metropolis_hastings_retained_sample_ids(
        chain_report=run_report.chain_report
    )
    prior_names = [row.prior_name for row in manifest.prior_rows]
    assert prior_names == [
        "branch-lengths",
        "substitution:kappa",
        "substitution:base-frequencies",
    ]
    assert manifest.proposal_schedule["kappa_move_weight"] == 1.0
    assert manifest.model_configuration["substitution_model_name"] == "HKY85"


def test_bayesian_run_manifest_round_trips_json(tmp_path: Path) -> None:
    run_report = _build_fixed_topology_dna_run_report(
        tree_path=fixture("trees", "hky85_likelihood_tree_2_taxa.nwk"),
        alignment_path=fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta"),
        seed=5,
    )
    sample_id_path = _write_sample_id_output(
        tmp_path / "retained-sample-ids.json",
        run_report,
    )
    manifest = build_fixed_topology_dna_run_manifest(
        run_report=run_report,
        tree_path=fixture("trees", "hky85_likelihood_tree_2_taxa.nwk"),
        alignment_path=fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta"),
        output_paths=[sample_id_path],
    )

    manifest_path = write_bayesian_run_manifest(
        tmp_path / "fixed-topology-dna.manifest.json",
        manifest,
    )
    loaded_manifest = load_bayesian_run_manifest(manifest_path)

    assert loaded_manifest == manifest


def test_replay_fixed_topology_dna_run_manifest_reproduces_retained_sample_ids(
    tmp_path: Path,
) -> None:
    tree_path = fixture("trees", "hky85_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta")
    run_report = _build_fixed_topology_dna_run_report(
        tree_path=tree_path,
        alignment_path=alignment_path,
        seed=7,
    )
    sample_id_path = _write_sample_id_output(
        tmp_path / "retained-sample-ids.json",
        run_report,
    )
    manifest = build_fixed_topology_dna_run_manifest(
        run_report=run_report,
        tree_path=tree_path,
        alignment_path=alignment_path,
        output_paths=[sample_id_path],
    )

    replay_report = replay_fixed_topology_dna_run_manifest(manifest)

    assert replay_report.run_kind == "fixed-topology-dna"
    assert replay_report.model_name == "HKY85"
    assert replay_report.matches_retained_sample_ids is True
    assert replay_report.expected_retained_sample_ids == manifest.retained_sample_ids
    assert replay_report.replayed_retained_sample_ids == manifest.retained_sample_ids


def test_replay_fixed_topology_dna_run_manifest_rejects_input_checksum_drift(
    tmp_path: Path,
) -> None:
    tree_copy = tmp_path / "tree.nwk"
    alignment_copy = tmp_path / "alignment.fasta"
    tree_copy.write_text(
        fixture("trees", "hky85_likelihood_tree_2_taxa.nwk").read_text(
            encoding="utf-8"
        ),
        encoding="utf-8",
    )
    alignment_copy.write_text(
        fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta").read_text(
            encoding="utf-8"
        ),
        encoding="utf-8",
    )
    run_report = _build_fixed_topology_dna_run_report(
        tree_path=tree_copy,
        alignment_path=alignment_copy,
        seed=11,
    )
    sample_id_path = _write_sample_id_output(
        tmp_path / "retained-sample-ids.json",
        run_report,
    )
    manifest = build_fixed_topology_dna_run_manifest(
        run_report=run_report,
        tree_path=tree_copy,
        alignment_path=alignment_copy,
        output_paths=[sample_id_path],
    )
    alignment_copy.write_text(
        ">A\nAAAA\n>B\nAAAA\n",
        encoding="utf-8",
    )

    with pytest.raises(PhylogeneticsError) as error_info:
        replay_fixed_topology_dna_run_manifest(manifest)

    assert (
        error_info.value.code == "bayesian_run_manifest_replay_input_checksum_mismatch"
    )


def _build_fixed_topology_dna_run_report(
    *,
    tree_path: Path,
    alignment_path: Path,
    seed: int,
) -> FixedTopologyDnaRunReport:
    model_definition = build_fixed_topology_dna_model_definition(
        substitution_model_name="HKY85",
        branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
        substitution_parameter_prior_bundle=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=1.5
            ),
            base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
                expected_component_names=("A", "C", "G", "T"),
                concentration_parameters={"A": 2.0, "C": 2.0, "G": 2.0, "T": 2.0},
            ),
        ),
    )
    proposal_schedule = build_fixed_topology_dna_proposal_schedule(
        model_definition=model_definition,
        branch_length_move_weight=1.0,
        branch_length_log_scale_standard_deviation=0.25,
        kappa_move_weight=1.0,
        kappa_log_scale_standard_deviation=0.35,
        base_frequency_move_weight=1.0,
        base_frequency_coordinate_standard_deviation=0.45,
    )
    return run_fixed_topology_dna_metropolis_hastings(
        tree=load_tree(tree_path),
        records=load_fasta_alignment(alignment_path),
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=9,
        sample_every=1,
        seed=seed,
    )


def _write_sample_id_output(path: Path, run_report: FixedTopologyDnaRunReport) -> Path:
    path.write_text(
        json.dumps(
            list_metropolis_hastings_retained_sample_ids(
                chain_report=run_report.chain_report
            ),
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return path
