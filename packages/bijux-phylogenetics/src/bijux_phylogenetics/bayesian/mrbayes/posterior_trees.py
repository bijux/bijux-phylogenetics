from __future__ import annotations

from pathlib import Path
import re

from bijux_phylogenetics.io.biopython import loads_biophylo
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.nexus_translate import (
    parse_nexus_translate_map,
    translate_nexus_tip_labels,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import EngineWorkflowError, TreeParseError
from bijux_phylogenetics.trees import (
    compute_clade_frequency_table,
    compute_consensus_tree,
    load_tree_set,
)

from .artifacts import _mrbayes_artifact_error
from .models import (
    MrBayesConsensusTreeReport,
    MrBayesPosteriorSummaryReport,
    MrBayesPosteriorTreeSample,
    MrBayesPosteriorTreeSetReport,
)

_MRBAYES_TREE_PATTERN = re.compile(
    r"tree\s+([^\s=]+)\s*=\s*(.+?);", flags=re.IGNORECASE | re.DOTALL
)
_MRBAYES_TREE_GENERATION_PATTERN = re.compile(r"(\d+)$")
_MRBAYES_PROBABILITY_PATTERN = re.compile(
    r'prob\s*=\s*"?(?P<value>[0-9.eE+-]+)"?',
    flags=re.IGNORECASE,
)
_MRBAYES_PROBABILITY_PERCENT_PATTERN = re.compile(
    r'prob\s*\(\s*percent\s*\)\s*=\s*"?(?P<value>[0-9.eE+-]+)"?',
    flags=re.IGNORECASE,
)


def _extract_mrbayes_tree_entries(text: str) -> list[tuple[str, str]]:
    entries = [
        (match.group(1), match.group(2).strip())
        for match in _MRBAYES_TREE_PATTERN.finditer(text)
    ]
    if not entries:
        raise EngineWorkflowError("MrBayes tree file contains no tree entries")
    return entries


def _strip_square_bracket_comments(text: str) -> str:
    stripped: list[str] = []
    depth = 0
    for character in text:
        if character == "[":
            depth += 1
            continue
        if character == "]" and depth:
            depth -= 1
            continue
        if depth == 0:
            stripped.append(character)
    return "".join(stripped)


def _detect_mrbayes_rooted_flag(tree_text: str) -> bool | None:
    prefix = tree_text.lstrip()
    if prefix.startswith("[&R]"):
        return True
    if prefix.startswith("[&U]"):
        return False
    return None


def _parse_mrbayes_tree_generation(tree_name: str) -> int | None:
    match = _MRBAYES_TREE_GENERATION_PATTERN.search(tree_name)
    return None if match is None else int(match.group(1))


def _parse_mrbayes_tree_text(
    tree_text: str, *, translation: dict[str, str]
) -> tuple[str, PhyloTree, bool | None]:
    rooted = _detect_mrbayes_rooted_flag(tree_text)
    stripped = _strip_square_bracket_comments(tree_text).strip()
    translated = translate_nexus_tip_labels(stripped, translation)
    tree = loads_biophylo(f"{translated};", source_format="newick")
    return dumps_newick(tree), tree, rooted


def parse_mrbayes_posterior_tree_samples(path: Path) -> MrBayesPosteriorTreeSetReport:
    """Parse a MrBayes posterior tree set into generation-tagged samples."""
    if not path.exists():
        raise _mrbayes_artifact_error(
            f"MrBayes posterior tree file was not found: {path}",
            code="mrbayes_tree_missing_file",
            path=path,
            artifact_kind="mrbayes-posterior-trees",
            details={"expected_section": "trees block"},
        )
    text = path.read_text(encoding="utf-8")
    try:
        entries = _extract_mrbayes_tree_entries(text)
    except EngineWorkflowError as error:
        raise _mrbayes_artifact_error(
            f"MrBayes posterior tree file contains no tree entries: {path}",
            code="mrbayes_tree_missing_entries",
            path=path,
            artifact_kind="mrbayes-posterior-trees",
            details={"expected_section": "trees block"},
        ) from error
    try:
        translation = parse_nexus_translate_map(text)
    except ValueError as error:
        raise _mrbayes_artifact_error(
            f"MrBayes posterior tree file contains an invalid translate block: {path}",
            code="mrbayes_tree_invalid_translate_block",
            path=path,
            artifact_kind="mrbayes-posterior-trees",
            details={
                "cause": str(error),
                "expected_section": "translate block",
            },
        ) from error
    samples: list[MrBayesPosteriorTreeSample] = []
    for tree_name, tree_text in entries:
        try:
            newick, tree, rooted = _parse_mrbayes_tree_text(
                tree_text, translation=translation
            )
        except TreeParseError as error:
            raise _mrbayes_artifact_error(
                f"MrBayes posterior tree entry '{tree_name}' could not be parsed: {path}",
                code="mrbayes_tree_parse_error",
                path=path,
                artifact_kind="mrbayes-posterior-trees",
                details={
                    "tree_name": tree_name,
                    "cause": error.message,
                    "expected_section": "trees block entry",
                },
            ) from error
        samples.append(
            MrBayesPosteriorTreeSample(
                tree_name=tree_name,
                generation=_parse_mrbayes_tree_generation(tree_name),
                rooted=rooted if rooted is not None else tree.rooted,
                tip_names=tree.tip_names,
                newick=newick,
            )
        )
    rooted_tree_count = sum(1 for sample in samples if sample.rooted)
    sampled_generations = [
        generation
        for generation in (sample.generation for sample in samples)
        if generation is not None
    ]
    return MrBayesPosteriorTreeSetReport(
        path=path,
        tree_count=len(samples),
        rooted_tree_count=rooted_tree_count,
        sampled_generations=sampled_generations,
        tip_names=samples[0].tip_names,
        trees=samples,
    )


def parse_mrbayes_consensus_tree(
    path: Path,
) -> tuple[PhyloTree, MrBayesConsensusTreeReport]:
    """Parse a MrBayes consensus tree with posterior-probability annotations."""
    if not path.exists():
        raise _mrbayes_artifact_error(
            f"MrBayes consensus tree file was not found: {path}",
            code="mrbayes_consensus_missing_file",
            path=path,
            artifact_kind="mrbayes-consensus-tree",
            details={"expected_section": "consensus tree"},
        )
    text = path.read_text(encoding="utf-8")
    try:
        translation = parse_nexus_translate_map(text)
    except ValueError as error:
        raise _mrbayes_artifact_error(
            f"MrBayes consensus tree file contains an invalid translate block: {path}",
            code="mrbayes_consensus_invalid_translate_block",
            path=path,
            artifact_kind="mrbayes-consensus-tree",
            details={
                "cause": str(error),
                "expected_section": "translate block",
            },
        ) from error
    try:
        entries = _extract_mrbayes_tree_entries(text)
    except EngineWorkflowError as error:
        raise _mrbayes_artifact_error(
            f"MrBayes consensus tree file contains no tree entries: {path}",
            code="mrbayes_consensus_missing_entries",
            path=path,
            artifact_kind="mrbayes-consensus-tree",
            details={"expected_section": "consensus tree"},
        ) from error
    if len(entries) != 1:
        raise _mrbayes_artifact_error(
            f"MrBayes consensus tree file must contain exactly one tree: {path}",
            code="mrbayes_consensus_invalid_tree_count",
            path=path,
            artifact_kind="mrbayes-consensus-tree",
            details={"observed_tree_count": len(entries), "expected_section": "tree"},
        )
    tree_name, tree_text = entries[0]
    try:
        consensus_newick, tree, rooted = _parse_mrbayes_tree_text(
            tree_text, translation=translation
        )
    except TreeParseError as error:
        raise _mrbayes_artifact_error(
            f"MrBayes consensus tree entry '{tree_name}' could not be parsed: {path}",
            code="mrbayes_consensus_parse_error",
            path=path,
            artifact_kind="mrbayes-consensus-tree",
            details={
                "tree_name": tree_name,
                "cause": error.message,
                "expected_section": "tree",
            },
        ) from error
    posterior_probabilities = [
        float(match.group("value"))
        for match in _MRBAYES_PROBABILITY_PATTERN.finditer(tree_text)
    ]
    posterior_probability_percents = [
        float(match.group("value"))
        for match in _MRBAYES_PROBABILITY_PERCENT_PATTERN.finditer(tree_text)
    ]
    report = MrBayesConsensusTreeReport(
        path=path,
        tree_name=tree_name,
        rooted=rooted if rooted is not None else tree.rooted,
        tip_names=tree.tip_names,
        consensus_newick=consensus_newick,
        annotated_node_count=len(posterior_probabilities),
        minimum_posterior_probability=(
            None if not posterior_probabilities else min(posterior_probabilities)
        ),
        maximum_posterior_probability=(
            None if not posterior_probabilities else max(posterior_probabilities)
        ),
        minimum_posterior_probability_percent=(
            None
            if not posterior_probability_percents
            else min(posterior_probability_percents)
        ),
        maximum_posterior_probability_percent=(
            None
            if not posterior_probability_percents
            else max(posterior_probability_percents)
        ),
    )
    return tree, report


def _write_mrbayes_posterior_tree_set(
    path: Path,
    *,
    trees: list[MrBayesPosteriorTreeSample],
) -> Path:
    path.write_text("".join(f"{sample.newick}\n" for sample in trees), encoding="utf-8")
    return path


def summarize_mrbayes_posterior_trees(
    tree_set_path: Path,
    *,
    burnin_fraction: float = 0.25,
) -> tuple[PhyloTree, MrBayesPosteriorSummaryReport]:
    """Summarize MrBayes posterior trees after discarding a burn-in fraction."""
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    tree_set_report = parse_mrbayes_posterior_tree_samples(tree_set_path)
    burnin_tree_count = int(tree_set_report.tree_count * burnin_fraction)
    kept_trees = tree_set_report.trees[burnin_tree_count:]
    if not kept_trees:
        raise EngineWorkflowError(
            f"MrBayes posterior tree file is empty after burn-in filtering: {tree_set_path}"
        )
    filtered_tree_set_path = tree_set_path.with_suffix(".postburnin.nwk")
    _write_mrbayes_posterior_tree_set(filtered_tree_set_path, trees=kept_trees)
    summary = load_tree_set(filtered_tree_set_path)
    consensus_tree, consensus = compute_consensus_tree(filtered_tree_set_path)
    clade_frequencies = compute_clade_frequency_table(filtered_tree_set_path)
    return consensus_tree, MrBayesPosteriorSummaryReport(
        source_path=tree_set_path,
        filtered_tree_set_path=filtered_tree_set_path,
        total_tree_count=tree_set_report.tree_count,
        burnin_tree_count=burnin_tree_count,
        kept_tree_count=len(kept_trees),
        rooted_topology_count=summary.rooted_topology_count,
        shared_taxa=summary.shared_taxa,
        consensus_newick=consensus.consensus_newick,
        clade_frequency_count=len(clade_frequencies.clade_frequencies),
    )
