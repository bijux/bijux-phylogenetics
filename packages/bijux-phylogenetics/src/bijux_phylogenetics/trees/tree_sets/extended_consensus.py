from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .contracts import MajorityRuleExtendedConsensusReport


def compute_majority_rule_extended_consensus(
    path: Path,
) -> tuple[PhyloTree, MajorityRuleExtendedConsensusReport]:
    raise NotImplementedError


def write_majority_rule_extended_consensus_inclusion_table(
    path: Path,
    report: MajorityRuleExtendedConsensusReport,
) -> Path:
    raise NotImplementedError


def write_majority_rule_extended_consensus_rejected_conflict_table(
    path: Path,
    report: MajorityRuleExtendedConsensusReport,
) -> Path:
    raise NotImplementedError


def write_majority_rule_extended_consensus_artifacts(
    out_dir: Path,
    report: MajorityRuleExtendedConsensusReport,
) -> dict[str, Path]:
    raise NotImplementedError
