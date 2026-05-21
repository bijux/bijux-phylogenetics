from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.ancestral.discrete.state_resolution import (
    resolve_clade_consensus_state,
)


def load_allowed_host_transitions(
    path: Path,
    *,
    observed_hosts: list[str],
) -> list[tuple[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"host-transition constraint file not found: {path}")
    raw_text = path.read_text(encoding="utf-8")
    sample = raw_text[:1024]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
    except csv.Error:
        dialect = csv.excel_tab if "\t" in sample else csv.excel
    reader = csv.DictReader(raw_text.splitlines(), dialect=dialect)
    if reader.fieldnames is None:
        raise ValueError("host-transition constraint file must contain a header row")
    required = {"source_host", "target_host"}
    if not required.issubset(reader.fieldnames):
        raise ValueError(
            "host-transition constraint file must contain source_host and target_host columns"
        )
    allowed_field = (
        "transition_allowed" if "transition_allowed" in reader.fieldnames else None
    )
    observed_host_set = set(observed_hosts)
    allowed_pairs: list[tuple[str, str]] = []
    for row in reader:
        source_host = (row.get("source_host") or "").strip()
        target_host = (row.get("target_host") or "").strip()
        if not source_host or not target_host:
            raise ValueError(
                "host-transition constraint rows must name both source_host and target_host"
            )
        if source_host == target_host:
            raise ValueError(
                "host-transition constraint rows must connect distinct hosts"
            )
        if source_host not in observed_host_set:
            raise ValueError(
                "host-transition source host is not present in the analyzed host vocabulary: "
                f"{source_host}"
            )
        if target_host not in observed_host_set:
            raise ValueError(
                "host-transition target host is not present in the analyzed host vocabulary: "
                f"{target_host}"
            )
        if allowed_field is not None and not parse_truthy_cell(
            row.get(allowed_field, "")
        ):
            continue
        allowed_pairs.append((source_host, target_host))
    if not allowed_pairs:
        raise ValueError(
            "host-transition constraint file must allow at least one directed host transition"
        )
    return sorted(set(allowed_pairs))


def parse_truthy_cell(raw: str) -> bool:
    normalized = raw.strip().lower()
    if normalized in {"", "0", "false", "no", "forbidden"}:
        return False
    if normalized in {"1", "true", "yes", "allowed", "x"}:
        return True
    raise ValueError(
        "host-transition constraint transition_allowed cells must be one of "
        "0,1,false,true,no,yes,forbidden,allowed,x"
    )


def transition_certainty_class(
    *,
    changed: bool,
    overlapping_hosts: list[str],
    parent_host_set: list[str],
    child_host_set: list[str],
) -> str:
    if not changed:
        return "no_switch"
    if overlapping_hosts:
        return "uncertain_switch"
    if len(parent_host_set) == 1 and len(child_host_set) == 1:
        return "certain_switch"
    return "uncertain_switch"


def node_signature(node) -> str:
    taxa = sorted(node_descendant_taxa(node))
    if taxa:
        return "|".join(taxa)
    return node.name or "<unnamed>"


def node_descendant_taxa(node) -> list[str]:
    if not node.children:
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(node_descendant_taxa(child))
    return taxa


def resolve_host_state(
    *,
    descendant_taxa: list[str],
    candidate_hosts: list[str],
    observed_hosts_by_taxon: dict[str, str],
    fallback_host: str,
) -> str:
    return resolve_clade_consensus_state(
        clade_taxa=descendant_taxa,
        candidate_states=candidate_hosts,
        observed_states_by_taxon=observed_hosts_by_taxon,
        fallback_state=fallback_host,
    )


def stable_float(value: float) -> float:
    return float(format(round(value, 15), ".15g"))


def format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(value)
