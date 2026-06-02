from __future__ import annotations

import tempfile

from bijux_phylogenetics.io.nexus_translate import (
    parse_nexus_translate_map,
    translate_nexus_tip_labels,
)

from ._shared import (
    _BEAST_TREE_PATTERN,
    _BEAST_TREE_STATE_PATTERN,
    EngineWorkflowError,
    InvalidAlignmentError,
    Path,
    PhyloTree,
    TreeParseError,
    _beast_artifact_error,
    compute_clade_frequency_table,
    compute_consensus_tree,
    descendant_taxa,
    dumps_newick,
    load_tree_set,
    loads_biophylo,
    re,
    summarize_posterior_topology_diversity,
)
from .models import (
    BeastPosteriorClade,
    BeastPosteriorConsensusReport,
    BeastPosteriorTopologyDiversityReport,
    BeastPosteriorTreeSample,
    BeastPosteriorTreeSetReport,
)


def parse_beast_posterior_tree_samples(
    path: Path,
    *,
    burnin_fraction: float = 0.0,
) -> BeastPosteriorTreeSetReport:
    """Parse a BEAST posterior tree set into state-tagged normalized trees."""
    if not path.exists():
        raise _beast_artifact_error(
            f"BEAST posterior tree file was not found: {path}",
            code="beast_tree_missing_file",
            path=path,
            artifact_kind="beast-posterior-trees",
            details={"expected_section": "trees block"},
        )
    text = path.read_text(encoding="utf-8")
    entries = _extract_beast_tree_entries(text)
    if not entries:
        raise _beast_artifact_error(
            f"BEAST posterior tree file contains no tree entries: {path}",
            code="beast_tree_missing_entries",
            path=path,
            artifact_kind="beast-posterior-trees",
            details={"expected_section": "trees block"},
        )
    burnin_tree_count, kept_entries = _split_beast_tree_entries(
        entries, burnin_fraction=burnin_fraction, path=path
    )
    try:
        translation = parse_nexus_translate_map(text)
    except ValueError as error:
        raise _beast_artifact_error(
            f"BEAST posterior tree file contains an invalid translate block: {path}",
            code="beast_tree_invalid_translate_block",
            path=path,
            artifact_kind="beast-posterior-trees",
            details={
                "cause": str(error),
                "expected_section": "translate block",
            },
        ) from error
    samples: list[BeastPosteriorTreeSample] = []
    trees: list[PhyloTree] = []
    for tree_name, tree_text in kept_entries:
        try:
            (
                newick,
                tree,
                rooted,
                annotation_values,
                annotation_keys,
                annotation_record_count,
            ) = _parse_beast_tree_text(tree_text, translation=translation)
        except TreeParseError as error:
            raise _beast_artifact_error(
                f"BEAST posterior tree entry '{tree_name}' could not be parsed: {path}",
                code="beast_tree_parse_error",
                path=path,
                artifact_kind="beast-posterior-trees",
                details={
                    "tree_name": tree_name,
                    "cause": error.message,
                    "expected_section": "trees block entry",
                },
            ) from error
        samples.append(
            BeastPosteriorTreeSample(
                tree_name=tree_name,
                state=_parse_beast_tree_state(tree_name),
                rooted=rooted if rooted is not None else True,
                tip_names=tree.tip_names,
                newick=newick,
                annotation_key_count=len(annotation_keys),
                annotation_record_count=annotation_record_count,
                annotation_keys=annotation_keys,
                annotation_values=annotation_values,
            )
        )
        trees.append(tree)
    clades = _summarize_beast_clades(trees)
    sampled_states = [
        state for state in (sample.state for sample in samples) if state is not None
    ]
    rooted_tree_count = sum(1 for sample in samples if sample.rooted)
    return BeastPosteriorTreeSetReport(
        path=path,
        burnin_fraction=burnin_fraction,
        total_tree_count=len(entries),
        burnin_tree_count=burnin_tree_count,
        kept_tree_count=len(samples),
        rooted_tree_count=rooted_tree_count,
        sampled_states=sampled_states,
        tip_names=sorted(samples[0].tip_names),
        clades=clades,
        trees=samples,
    )


def write_beast_posterior_tree_set(
    path: Path, report: BeastPosteriorTreeSetReport
) -> Path:
    """Write a normalized Newick tree set from parsed BEAST posterior samples."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(f"{sample.newick}\n" for sample in report.trees),
        encoding="utf-8",
    )
    return path


def _annotate_consensus_tree_with_posterior_probabilities(tree: PhyloTree) -> None:
    for node in tree.iter_nodes():
        if node is tree.root or node.is_leaf() or node.name is None:
            continue
        node.name = format(float(node.name) / 100.0, ".15g")


def summarize_beast_posterior_trees(
    tree_set_path: Path,
    *,
    burnin_fraction: float = 0.25,
) -> tuple[PhyloTree, BeastPosteriorConsensusReport]:
    """Summarize BEAST posterior trees into a majority-rule consensus tree."""
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    tree_set_report = parse_beast_posterior_tree_samples(
        tree_set_path,
        burnin_fraction=burnin_fraction,
    )
    if not tree_set_report.trees:
        raise EngineWorkflowError(
            f"BEAST posterior tree file is empty after burn-in filtering: {tree_set_path}"
        )
    retained_tree_set_path = _retained_tree_set_path(
        tree_set_path,
        burnin_fraction=tree_set_report.burnin_fraction,
    )
    write_beast_posterior_tree_set(retained_tree_set_path, tree_set_report)
    summary = load_tree_set(retained_tree_set_path)
    consensus_tree, consensus = compute_consensus_tree(retained_tree_set_path)
    _annotate_consensus_tree_with_posterior_probabilities(consensus_tree)
    consensus_newick = dumps_newick(consensus_tree)
    clade_frequencies = compute_clade_frequency_table(retained_tree_set_path)
    posterior_probabilities = sorted(
        float(node.name)
        for node in consensus_tree.iter_nodes()
        if node is not consensus_tree.root
        and not node.is_leaf()
        and node.name is not None
    )
    return consensus_tree, BeastPosteriorConsensusReport(
        source_path=tree_set_path,
        retained_tree_set_path=retained_tree_set_path,
        burnin_fraction=tree_set_report.burnin_fraction,
        total_tree_count=tree_set_report.total_tree_count,
        burnin_tree_count=tree_set_report.burnin_tree_count,
        kept_tree_count=tree_set_report.kept_tree_count,
        rooted_topology_count=summary.rooted_topology_count,
        shared_taxa=summary.shared_taxa,
        consensus_newick=consensus_newick,
        clade_frequency_count=len(clade_frequencies.clade_frequencies),
        annotated_node_count=len(posterior_probabilities),
        minimum_posterior_probability=(
            None if not posterior_probabilities else posterior_probabilities[0]
        ),
        maximum_posterior_probability=(
            None if not posterior_probabilities else posterior_probabilities[-1]
        ),
    )


def summarize_beast_posterior_topology_diversity(
    tree_set_path: Path,
    *,
    burnin_fraction: float = 0.25,
) -> BeastPosteriorTopologyDiversityReport:
    """Summarize BEAST posterior topology diversity after burn-in filtering."""
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    tree_set_report = parse_beast_posterior_tree_samples(
        tree_set_path,
        burnin_fraction=burnin_fraction,
    )
    if not tree_set_report.trees:
        raise EngineWorkflowError(
            f"BEAST posterior tree file is empty after burn-in filtering: {tree_set_path}"
        )
    retained_tree_set_path = _retained_tree_set_path(
        tree_set_path,
        burnin_fraction=tree_set_report.burnin_fraction,
    )
    write_beast_posterior_tree_set(retained_tree_set_path, tree_set_report)
    diversity = summarize_posterior_topology_diversity(retained_tree_set_path)
    return BeastPosteriorTopologyDiversityReport(
        source_path=tree_set_path,
        retained_tree_set_path=retained_tree_set_path,
        burnin_fraction=tree_set_report.burnin_fraction,
        total_tree_count=tree_set_report.total_tree_count,
        burnin_tree_count=tree_set_report.burnin_tree_count,
        kept_tree_count=tree_set_report.kept_tree_count,
        rooted_topology_count=diversity.rooted_topology_count,
        dominant_topology_frequency=diversity.dominant_topology_frequency,
        effective_topology_count=diversity.effective_topology_count,
        pair_count=diversity.pair_count,
        mean_robinson_foulds_distance=diversity.mean_robinson_foulds_distance,
        mean_normalized_robinson_foulds_distance=(
            diversity.mean_normalized_robinson_foulds_distance
        ),
        maximum_robinson_foulds_distance=diversity.maximum_robinson_foulds_distance,
        maximum_normalized_robinson_foulds_distance=(
            diversity.maximum_normalized_robinson_foulds_distance
        ),
        unstable_clade_count=diversity.unstable_clade_count,
    )


def _extract_beast_tree_entries(text: str) -> list[tuple[str, str]]:
    return [
        (match.group(1), match.group(2).strip())
        for match in _BEAST_TREE_PATTERN.finditer(text)
    ]


def _retained_tree_set_path(tree_set_path: Path, *, burnin_fraction: float) -> Path:
    fraction_token = (
        format(burnin_fraction, ".6f").rstrip("0").rstrip(".").replace(".", "p") or "0"
    )
    return Path(
        tempfile.mkstemp(
            prefix=f"{tree_set_path.stem}.burnin-{fraction_token}-",
            suffix=".nwk",
        )[1]
    )


def _split_beast_tree_entries(
    entries: list[tuple[str, str]],
    *,
    burnin_fraction: float,
    path: Path,
) -> tuple[int, list[tuple[str, str]]]:
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    burnin_tree_count = int(len(entries) * burnin_fraction)
    kept_entries = entries[burnin_tree_count:]
    if not kept_entries:
        raise _beast_artifact_error(
            f"BEAST posterior tree file is empty after burn-in filtering: {path}",
            code="beast_tree_empty_after_burnin",
            path=path,
            artifact_kind="beast-posterior-trees",
            details={
                "burnin_fraction": burnin_fraction,
                "total_tree_count": len(entries),
            },
        )
    return burnin_tree_count, kept_entries


def _parse_beast_tree_text(
    tree_text: str, *, translation: dict[str, str]
) -> tuple[str, PhyloTree, bool | None, dict[str, str], list[str], int]:
    rooted = _detect_nexus_rooted_flag(tree_text)
    annotation_values, annotation_keys, annotation_record_count = (
        _extract_beast_tree_annotations(tree_text)
    )
    stripped = _strip_square_bracket_comments(tree_text).strip()
    translated = translate_nexus_tip_labels(stripped, translation)
    tree = loads_biophylo(f"{translated};", source_format="newick")
    return (
        dumps_newick(tree),
        tree,
        rooted,
        annotation_values,
        annotation_keys,
        annotation_record_count,
    )


def _detect_nexus_rooted_flag(tree_text: str) -> bool | None:
    prefix = tree_text.lstrip()
    if prefix.startswith("[&R]"):
        return True
    if prefix.startswith("[&U]"):
        return False
    return None


def _parse_beast_tree_state(tree_name: str) -> int | None:
    match = _BEAST_TREE_STATE_PATTERN.search(tree_name)
    return None if match is None else int(match.group(1))


def _summarize_beast_clades(trees: list[PhyloTree]) -> list[BeastPosteriorClade]:
    if not trees:
        return []
    taxa_sets = {frozenset(tree.tip_names) for tree in trees}
    if len(taxa_sets) != 1:
        raise InvalidAlignmentError(
            "BEAST posterior tree summaries require all trees to share the exact same taxon set"
        )
    shared_taxa = set(next(iter(taxa_sets)))
    counts: dict[frozenset[str], int] = {}
    for tree in trees:
        for clade in _informative_tree_clades(tree, shared_taxa):
            counts[clade] = counts.get(clade, 0) + 1
    return [
        BeastPosteriorClade(
            clade="|".join(sorted(clade)),
            tree_count=tree_count,
            frequency=round(tree_count / len(trees), 15),
        )
        for clade, tree_count in sorted(
            counts.items(),
            key=lambda item: (-item[1], sorted(item[0])),
        )
    ]


def _informative_tree_clades(
    tree: PhyloTree, shared_taxa: set[str]
) -> list[frozenset[str]]:
    clades: list[frozenset[str]] = []
    for node in tree.iter_nodes():
        taxa = frozenset(descendant_taxa(node))
        if 1 < len(taxa) < len(shared_taxa):
            clades.append(taxa)
    return clades


def _extract_beast_tree_annotations(
    tree_text: str,
) -> tuple[dict[str, str], list[str], int]:
    annotation_values: dict[str, str] = {}
    annotation_keys: list[str] = []
    annotation_record_count = 0
    for match in re.finditer(r"\[(.*?)\]", tree_text):
        raw = match.group(1).strip()
        if not raw.startswith("&"):
            continue
        directive = raw[1:].strip()
        if directive in {"R", "U"}:
            continue
        for token in _split_beast_annotation_tokens(directive):
            if not token:
                continue
            if "=" in token:
                key, value = token.split("=", 1)
                normalized_value = value.strip()
            else:
                key = token
                normalized_value = "true"
            normalized_key = key.strip()
            if not normalized_key:
                continue
            if normalized_key not in annotation_values:
                annotation_keys.append(normalized_key)
            annotation_values[normalized_key] = normalized_value
            annotation_record_count += 1
    return annotation_values, sorted(annotation_keys), annotation_record_count


def _split_beast_annotation_tokens(text: str) -> list[str]:
    tokens: list[str] = []
    current: list[str] = []
    brace_depth = 0
    quote_character: str | None = None
    index = 0
    while index < len(text):
        char = text[index]
        if char in {"'", '"'}:
            if quote_character is None:
                quote_character = char
            elif char == quote_character:
                if (
                    quote_character == "'"
                    and index + 1 < len(text)
                    and text[index + 1] == "'"
                ):
                    current.append("''")
                    index += 2
                    continue
                quote_character = None
        if char == "{":
            brace_depth += 1
        elif char == "}" and brace_depth:
            brace_depth -= 1
        if char == "," and brace_depth == 0 and quote_character is None:
            token = "".join(current).strip()
            if token:
                tokens.append(token)
            current = []
            index += 1
            continue
        current.append(char)
        index += 1
    token = "".join(current).strip()
    if token:
        tokens.append(token)
    return tokens


def _strip_square_bracket_comments(text: str) -> str:
    result: list[str] = []
    depth = 0
    for char in text:
        if char == "[":
            depth += 1
            continue
        if char == "]" and depth:
            depth -= 1
            continue
        if depth == 0:
            result.append(char)
    return "".join(result)
