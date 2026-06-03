from __future__ import annotations

from .artifact_outputs import (
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_count_table,
    write_host_switch_exclusion_table,
    write_host_switch_fit_table,
    write_host_switch_summary_table,
    write_unsupported_host_switch_claim_table,
)
from .builder import summarize_host_switching
from .contracts import (
    HostStateNodeRow,
    HostSwitchBranchRow,
    HostSwitchCountRow,
    HostSwitchExclusionRow,
    HostSwitchFitRow,
    HostSwitchingReport,
    HostSwitchSummary,
    UnsupportedHostSwitchClaimRow,
)

__all__ = [
    "HostStateNodeRow",
    "HostSwitchBranchRow",
    "HostSwitchCountRow",
    "HostSwitchExclusionRow",
    "HostSwitchFitRow",
    "HostSwitchSummary",
    "HostSwitchingReport",
    "UnsupportedHostSwitchClaimRow",
    "summarize_host_switching",
    "write_host_state_node_table",
    "write_host_switch_branch_table",
    "write_host_switch_count_table",
    "write_host_switch_exclusion_table",
    "write_host_switch_fit_table",
    "write_host_switch_summary_table",
    "write_unsupported_host_switch_claim_table",
]
