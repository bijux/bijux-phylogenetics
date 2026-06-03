from __future__ import annotations

from pathlib import Path
import random

from bijux_phylogenetics.io.fasta import write_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment import AlignmentRecord


def _simulate_alignment_records(
    tree_path: Path,
    *,
    alphabet: tuple[str, ...],
    model: str,
    sequence_length: int,
    substitution_rate: float,
    seed: int,
):
    from .._state_propagation import _iter_tip_trait_values
    from .._stochastic import _poisson_count
    from ..contracts import AlignmentSimulationReport

    if sequence_length < 1:
        raise ValueError(f"sequence_length must be at least 1, got {sequence_length}")
    if substitution_rate < 0.0:
        raise ValueError(
            f"substitution_rate must be nonnegative, got {substitution_rate}"
        )
    tree = load_tree(tree_path)
    rng = random.Random(seed)  # nosec B311
    root_sequence = "".join(rng.choice(alphabet) for _ in range(sequence_length))

    def mutate_sequence(sequence: str, branch_length: float) -> str:
        residues: list[str] = []
        for residue in sequence:
            next_residue = residue
            for _ in range(_poisson_count(substitution_rate * branch_length, rng)):
                alternatives = [
                    candidate for candidate in alphabet if candidate != next_residue
                ]
                next_residue = rng.choice(alternatives)
            residues.append(next_residue)
        return "".join(residues)

    values = _iter_tip_trait_values(
        tree,
        root_state=root_sequence,
        propagate=lambda state, branch_length: mutate_sequence(state, branch_length),
    )
    return AlignmentSimulationReport(
        model=model,
        tree_path=tree_path,
        tip_count=tree.tip_count,
        seed=seed,
        sequence_length=sequence_length,
        substitution_rate=substitution_rate,
        inferred_alphabet="dna" if alphabet == ("A", "C", "G", "T") else "protein",
        records=[
            AlignmentRecord(identifier=taxon, sequence=sequence)
            for taxon, sequence in sorted(values.items())
        ],
    )


def simulate_dna_alignment(
    tree_path: Path,
    *,
    sequence_length: int,
    substitution_rate: float = 1.0,
    seed: int = 1,
):
    return _simulate_alignment_records(
        tree_path,
        alphabet=("A", "C", "G", "T"),
        model="jukes-cantor-like",
        sequence_length=sequence_length,
        substitution_rate=substitution_rate,
        seed=seed,
    )


def simulate_protein_alignment(
    tree_path: Path,
    *,
    sequence_length: int,
    substitution_rate: float = 1.0,
    seed: int = 1,
):
    return _simulate_alignment_records(
        tree_path,
        alphabet=tuple("ACDEFGHIKLMNPQRSTVWY"),
        model="symmetric-protein",
        sequence_length=sequence_length,
        substitution_rate=substitution_rate,
        seed=seed,
    )


def write_simulated_alignment(path: Path, report) -> Path:
    return write_fasta_alignment(path, report.records)
