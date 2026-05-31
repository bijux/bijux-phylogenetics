from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import random

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import loads_newick, write_newick, write_newick_tree_set
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodBootstrapCladeSupportRow,
    NucleotideLikelihoodBootstrapReplicateRow,
    NucleotideLikelihoodBootstrapTreeInferenceReport,
    NucleotideLikelihoodTreeInferenceReport,
)
from bijux_phylogenetics.phylo.likelihood.tree_inference import (
    infer_nucleotide_likelihood_tree,
    validate_nucleotide_likelihood_tree_inference_model_name,
    validate_nucleotide_likelihood_tree_inference_model_selection_criterion,
    validate_nucleotide_likelihood_tree_inference_search_method,
    validate_nucleotide_likelihood_tree_inference_start_tree_count,
)
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_clade_id,
    informative_rooted_clade_nodes,
    split_sort_key,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.trees import (
    compute_clade_frequency_table,
    compute_consensus_tree,
    write_clade_frequency_table,
)


def validate_nucleotide_likelihood_bootstrap_replicate_count(replicate_count: int) -> int:
    """Validate the native ML bootstrap replicate count."""
    if replicate_count < 1:
        raise ValueError("replicate_count must be at least one")
    return replicate_count


def bootstrap_nucleotide_likelihood_tree_inference(
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str = "auto",
    model_selection_criterion: str = "aic",
    search_method: str = "nni",
    start_tree_count: int = 4,
    start_tree_seed: int = 1,
    branch_reoptimization_policy: str = "coordinate-branch-lengths",
    replicate_count: int = 100,
    bootstrap_seed: int = 1,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodBootstrapTreeInferenceReport:
    """Infer one native ML tree and bootstrap its clade support by site resampling."""
    resolved_records, resolved_alignment_path = _resolve_bootstrap_alignment(records)
    normalized_model_name = validate_nucleotide_likelihood_tree_inference_model_name(
        model_name
    )
    normalized_model_selection_criterion = (
        validate_nucleotide_likelihood_tree_inference_model_selection_criterion(
            model_selection_criterion
        )
    )
    normalized_search_method = (
        validate_nucleotide_likelihood_tree_inference_search_method(search_method)
    )
    validated_start_tree_count = (
        validate_nucleotide_likelihood_tree_inference_start_tree_count(start_tree_count)
    )
    validated_replicate_count = (
        validate_nucleotide_likelihood_bootstrap_replicate_count(replicate_count)
    )

    reference_report = infer_nucleotide_likelihood_tree(
        resolved_records,
        model_name=normalized_model_name,
        model_selection_criterion=normalized_model_selection_criterion,
        search_method=normalized_search_method,
        start_tree_count=validated_start_tree_count,
        start_tree_seed=start_tree_seed,
        branch_reoptimization_policy=branch_reoptimization_policy,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    replicate_model_name = reference_report.selected_model_name.lower()
    rng = random.Random(bootstrap_seed)  # nosec B311
    replicate_rows: list[NucleotideLikelihoodBootstrapReplicateRow] = []
    replicate_trees: list[PhyloTree] = []
    for replicate_index in range(1, validated_replicate_count + 1):
        sampled_site_indices, replicate_records = _bootstrap_resample_alignment_columns(
            resolved_records,
            rng=rng,
        )
        replicate_start_tree_seed = (
            start_tree_seed + (replicate_index - 1) * validated_start_tree_count
        )
        replicate_report = infer_nucleotide_likelihood_tree(
            replicate_records,
            model_name=replicate_model_name,
            model_selection_criterion=normalized_model_selection_criterion,
            search_method=normalized_search_method,
            start_tree_count=validated_start_tree_count,
            start_tree_seed=replicate_start_tree_seed,
            branch_reoptimization_policy=branch_reoptimization_policy,
            lower_branch_length_bound=lower_branch_length_bound,
            upper_branch_length_bound=upper_branch_length_bound,
            improvement_tolerance=improvement_tolerance,
            max_coordinate_passes=max_coordinate_passes,
        )
        replicate_trees.append(loads_newick(replicate_report.best_final_tree_newick))
        replicate_rows.append(
            NucleotideLikelihoodBootstrapReplicateRow(
                replicate_index=replicate_index,
                sampled_site_indices=sampled_site_indices,
                replicate_start_tree_seed=replicate_start_tree_seed,
                selected_model_name=replicate_report.selected_model_name,
                best_run_source_label=replicate_report.best_run_source_label,
                final_tree_newick=replicate_report.best_final_tree_newick,
                final_log_likelihood=replicate_report.best_final_log_likelihood,
                final_topology_fingerprint=replicate_report.best_final_topology_fingerprint,
                accepted_move_count=replicate_report.best_search_report.accepted_move_count,
                search_iteration_count=_resolve_bootstrap_search_iteration_count(
                    replicate_report
                ),
            )
        )
    clade_support_rows = _build_bootstrap_clade_support_rows(
        reference_tree=loads_newick(reference_report.best_final_tree_newick),
        replicate_trees=replicate_trees,
    )
    return NucleotideLikelihoodBootstrapTreeInferenceReport(
        algorithm="nucleotide-likelihood-bootstrap-tree-inference",
        alignment_path=None
        if resolved_alignment_path is None
        else str(resolved_alignment_path),
        requested_model_name=normalized_model_name,
        model_selection_strategy=reference_report.model_selection_strategy,
        model_selection_criterion=reference_report.model_selection_criterion,
        selected_reference_model_name=reference_report.selected_model_name,
        search_method=normalized_search_method,
        branch_reoptimization_policy=branch_reoptimization_policy,
        taxon_count=reference_report.taxon_count,
        site_count=reference_report.site_count,
        pattern_count=reference_report.pattern_count,
        start_tree_count=validated_start_tree_count,
        start_tree_seed=start_tree_seed,
        replicate_count=validated_replicate_count,
        bootstrap_seed=bootstrap_seed,
        reference_tree_newick=reference_report.best_final_tree_newick,
        reference_log_likelihood=reference_report.best_final_log_likelihood,
        reference_topology_fingerprint=reference_report.best_final_topology_fingerprint,
        reference_best_run_source_label=reference_report.best_run_source_label,
        replicate_rows=replicate_rows,
        clade_support_rows=clade_support_rows,
    )


def bootstrap_nucleotide_likelihood_tree_inference_from_alignment(
    alignment_path: Path,
    *,
    model_name: str = "auto",
    model_selection_criterion: str = "aic",
    search_method: str = "nni",
    start_tree_count: int = 4,
    start_tree_seed: int = 1,
    branch_reoptimization_policy: str = "coordinate-branch-lengths",
    replicate_count: int = 100,
    bootstrap_seed: int = 1,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodBootstrapTreeInferenceReport:
    """Infer one native ML tree and bootstrap it from one FASTA alignment path."""
    return bootstrap_nucleotide_likelihood_tree_inference(
        alignment_path,
        model_name=model_name,
        model_selection_criterion=model_selection_criterion,
        search_method=search_method,
        start_tree_count=start_tree_count,
        start_tree_seed=start_tree_seed,
        branch_reoptimization_policy=branch_reoptimization_policy,
        replicate_count=replicate_count,
        bootstrap_seed=bootstrap_seed,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def write_nucleotide_likelihood_bootstrap_replicate_draws_table(
    path: Path,
    report: NucleotideLikelihoodBootstrapTreeInferenceReport,
) -> Path:
    """Write one deterministic native likelihood bootstrap replicate ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "replicate_index\tsampled_site_indices\treplicate_start_tree_seed\tselected_model_name\tbest_run_source_label\tfinal_log_likelihood\tfinal_topology_fingerprint\taccepted_move_count\tsearch_iteration_count\tfinal_tree_newick"
    ]
    for row in report.replicate_rows:
        lines.append(
            (
                f"{row.replicate_index}\t"
                + ",".join(str(index) for index in row.sampled_site_indices)
                + f"\t{row.replicate_start_tree_seed}"
                + f"\t{row.selected_model_name}"
                + f"\t{row.best_run_source_label}"
                + f"\t{format(row.final_log_likelihood, '.15g')}"
                + f"\t{row.final_topology_fingerprint}"
                + f"\t{row.accepted_move_count}"
                + f"\t{row.search_iteration_count}"
                + f"\t{row.final_tree_newick}"
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_bootstrap_clade_support_table(
    path: Path,
    report: NucleotideLikelihoodBootstrapTreeInferenceReport,
) -> Path:
    """Write native bootstrap clade support mapped by descendant-tip set."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "branch_id\tnode_label\tdescendant_taxa\tsupporting_tree_count\tclade_frequency\tsupport_percent"
    ]
    for row in report.clade_support_rows:
        lines.append(
            (
                f"{row.branch_id}\t"
                f"{'' if row.node_label is None else row.node_label}\t"
                + "|".join(row.descendant_taxa)
                + f"\t{row.supporting_tree_count}"
                + f"\t{format(row.clade_frequency, '.15g')}"
                + f"\t{format(row.support_percent, '.15g')}"
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_bootstrap_run_json(
    path: Path,
    report: NucleotideLikelihoodBootstrapTreeInferenceReport,
) -> Path:
    """Write one governed native likelihood bootstrap run payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_nucleotide_likelihood_bootstrap_artifacts(
    out_dir: Path,
    report: NucleotideLikelihoodBootstrapTreeInferenceReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one native ML bootstrap run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    reference_tree_path = write_newick(
        out_dir / "final_tree.nwk",
        loads_newick(report.reference_tree_newick),
    )
    replicate_trees_path = write_newick_tree_set(
        out_dir / "replicate_trees.nwk",
        [loads_newick(row.final_tree_newick) for row in report.replicate_rows],
    )
    replicate_draws_path = write_nucleotide_likelihood_bootstrap_replicate_draws_table(
        out_dir / "replicate_draws.tsv",
        report,
    )
    clade_support_path = write_nucleotide_likelihood_bootstrap_clade_support_table(
        out_dir / "clade_support.tsv",
        report,
    )
    consensus_tree, _consensus_report = compute_consensus_tree(replicate_trees_path)
    consensus_tree_path = write_newick(out_dir / "consensus_tree.nwk", consensus_tree)
    clade_frequency_report = compute_clade_frequency_table(replicate_trees_path)
    clade_frequencies_path = write_clade_frequency_table(
        out_dir / "clade_frequencies.tsv",
        clade_frequency_report,
    )
    run_json_path = write_nucleotide_likelihood_bootstrap_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "reference_tree_path": reference_tree_path,
        "replicate_trees_path": replicate_trees_path,
        "replicate_draws_path": replicate_draws_path,
        "clade_support_path": clade_support_path,
        "consensus_tree_path": consensus_tree_path,
        "clade_frequencies_path": clade_frequencies_path,
        "run_json_path": run_json_path,
    }


def _resolve_bootstrap_alignment(
    records: list[AlignmentRecord] | Path,
) -> tuple[list[AlignmentRecord], Path | None]:
    if isinstance(records, Path):
        return load_fasta_alignment(records), records
    return records, None


def _bootstrap_resample_alignment_columns(
    records: list[AlignmentRecord],
    *,
    rng: random.Random,
) -> tuple[list[int], list[AlignmentRecord]]:
    sampled_site_indices = [
        rng.randrange(len(records[0].sequence)) for _ in range(len(records[0].sequence))
    ]
    return sampled_site_indices, [
        AlignmentRecord(
            identifier=record.identifier,
            sequence="".join(
                record.sequence[position] for position in sampled_site_indices
            ),
        )
        for record in records
    ]


def _resolve_bootstrap_search_iteration_count(
    report: NucleotideLikelihoodTreeInferenceReport,
) -> int:
    best_search_report = report.best_search_report
    if hasattr(best_search_report, "iteration_count"):
        return int(best_search_report.iteration_count)
    return max((row.iteration for row in best_search_report.trace_rows), default=0)


def _build_bootstrap_clade_support_rows(
    *,
    reference_tree: PhyloTree,
    replicate_trees: list[PhyloTree],
) -> list[NucleotideLikelihoodBootstrapCladeSupportRow]:
    shared_taxa = {
        node.name for node in reference_tree.iter_leaves() if node.name is not None
    }
    replicate_clade_sets = [
        set(informative_rooted_clade_nodes(tree, shared_taxa).keys())
        for tree in replicate_trees
    ]
    rows: list[NucleotideLikelihoodBootstrapCladeSupportRow] = []
    reference_clades = informative_rooted_clade_nodes(reference_tree, shared_taxa)
    for clade_signature, node in sorted(
        reference_clades.items(),
        key=lambda item: split_sort_key(item[0]),
    ):
        supporting_tree_count = sum(
            1
            for replicate_clades in replicate_clade_sets
            if clade_signature in replicate_clades
        )
        clade_frequency = round(supporting_tree_count / len(replicate_trees), 15)
        rows.append(
            NucleotideLikelihoodBootstrapCladeSupportRow(
                branch_id=canonical_clade_id(clade_signature),
                node_label=node.name,
                descendant_taxa=sorted(clade_signature),
                supporting_tree_count=supporting_tree_count,
                clade_frequency=clade_frequency,
                support_percent=round(clade_frequency * 100.0, 15),
            )
        )
    return rows
