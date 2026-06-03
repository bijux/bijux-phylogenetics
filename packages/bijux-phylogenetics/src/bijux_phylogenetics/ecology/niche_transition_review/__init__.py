from __future__ import annotations

from .artifact_outputs import (
    write_niche_state_node_table,
    write_niche_transition_branch_table,
    write_niche_transition_clade_table,
    write_niche_transition_count_table,
    write_niche_transition_exclusion_table,
    write_niche_transition_rate_table,
    write_niche_transition_summary_table,
)
from .builder import summarize_niche_transitions
from .contracts import (
    NicheStateNodeRow,
    NicheTransitionBranchRow,
    NicheTransitionCladeRow,
    NicheTransitionCountRow,
    NicheTransitionExclusionRow,
    NicheTransitionRateRow,
    NicheTransitionReport,
    NicheTransitionSummary,
)
