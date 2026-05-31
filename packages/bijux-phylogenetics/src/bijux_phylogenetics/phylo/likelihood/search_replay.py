from __future__ import annotations

from bijux_phylogenetics.io.newick import dumps_newick, loads_newick
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodNniSearchReport,
    NucleotideLikelihoodSearchTraceReplayReport,
    NucleotideLikelihoodSearchTraceReplayStep,
)
from bijux_phylogenetics.phylo.topology import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.rooted_nni import (
    apply_rooted_nni_move,
    iter_rooted_nni_move_candidates,
)


def replay_nucleotide_likelihood_nni_search_trace(
    report: NucleotideLikelihoodNniSearchReport,
) -> NucleotideLikelihoodSearchTraceReplayReport:
    """Replay accepted rooted NNI moves from one native likelihood search trace."""
    current_tree = loads_newick(report.start_tree_newick).refresh()
    step_rows: list[NucleotideLikelihoodSearchTraceReplayStep] = []
    failure_reason: str | None = None
    for row in report.trace_rows:
        if row.event_kind != "accepted-move":
            continue
        matches = []
        for candidate in iter_rooted_nni_move_candidates(current_tree):
            if (
                candidate.pivot_branch_id != row.pivot_branch_id
                or candidate.sibling_clade_id != row.sibling_clade_id
                or candidate.exchanged_clade_id != row.exchanged_clade_id
            ):
                continue
            replayed_tree = apply_rooted_nni_move(current_tree, candidate).refresh()
            replayed_topology_fingerprint = rooted_topology_fingerprint(replayed_tree)
            if replayed_topology_fingerprint == row.candidate_topology_fingerprint:
                matches.append((candidate, replayed_tree, replayed_topology_fingerprint))
        step_rows.append(
            NucleotideLikelihoodSearchTraceReplayStep(
                trace_event_index=row.event_index,
                iteration=row.iteration,
                move_type=row.move_type,
                trace_reason=row.trace_reason,
                candidate_topology_fingerprint=row.candidate_topology_fingerprint,
                replayed_topology_fingerprint=(
                    None if len(matches) != 1 else matches[0][2]
                ),
                matched_candidate_count=len(matches),
                step_replayed=len(matches) == 1,
            )
        )
        if len(matches) != 1:
            failure_reason = (
                "accepted trace row did not resolve to exactly one rooted NNI move"
            )
            break
        current_tree = matches[0][1]
    return _build_search_trace_replay_report(
        source_search_algorithm=report.algorithm,
        start_tree_newick=report.start_tree_newick,
        stored_final_tree_newick=report.final_tree_newick,
        current_tree=current_tree,
        failure_reason=failure_reason,
        step_rows=step_rows,
    )

def _build_search_trace_replay_report(
    *,
    source_search_algorithm: str,
    start_tree_newick: str,
    stored_final_tree_newick: str,
    current_tree,
    failure_reason: str | None,
    step_rows: list[NucleotideLikelihoodSearchTraceReplayStep],
) -> NucleotideLikelihoodSearchTraceReplayReport:
    stored_final_topology_fingerprint = rooted_topology_fingerprint(
        loads_newick(stored_final_tree_newick)
    )
    replayed_final_tree_newick = dumps_newick(current_tree)
    replayed_final_topology_fingerprint = rooted_topology_fingerprint(current_tree)
    final_topology_matches = (
        failure_reason is None
        and replayed_final_topology_fingerprint == stored_final_topology_fingerprint
    )
    return NucleotideLikelihoodSearchTraceReplayReport(
        algorithm="nucleotide-likelihood-search-trace-replay",
        source_search_algorithm=source_search_algorithm,
        accepted_trace_event_count=len(step_rows),
        replayed_step_count=sum(1 for row in step_rows if row.step_replayed),
        start_tree_newick=start_tree_newick,
        stored_final_tree_newick=stored_final_tree_newick,
        stored_final_topology_fingerprint=stored_final_topology_fingerprint,
        replayed_final_tree_newick=replayed_final_tree_newick,
        replayed_final_topology_fingerprint=replayed_final_topology_fingerprint,
        final_topology_matches=final_topology_matches,
        replay_failed=failure_reason is not None,
        failure_reason=failure_reason,
        step_rows=step_rows,
    )
