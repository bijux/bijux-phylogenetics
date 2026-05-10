from __future__ import annotations

import json
from pathlib import Path


STUDY_ID = "primate-longevity-signal"
EVIDENCE_ID = "evidence-001"
SOURCE_LOCATOR = "external:lund/pcm1-plots-signal/script"

FAMILY_DEFINITIONS = {
    "workflow-contracts": {
        "title": "Workflow contracts",
        "summary": "Reproducibility and package-loading context that stays visible without being overstated as analytical parity.",
    },
    "data-preparation": {
        "title": "Data preparation",
        "summary": "Workbook-derived preprocessing, tip alignment, and vector assembly that must match before downstream inference is credible.",
    },
    "tree-operations": {
        "title": "Tree operations",
        "summary": "Tree import, pruning, and topology-manipulation steps checked by stable taxon or node identity.",
    },
    "visual-surfaces": {
        "title": "Visual surfaces",
        "summary": "Plotting-oriented teaching surfaces tracked honestly without claiming rendered-figure equivalence.",
    },
    "simulation-inputs": {
        "title": "Simulation inputs",
        "summary": "Seeded random inputs frozen for cross-tool comparison before fitting behavior is judged.",
    },
    "comparative-signal": {
        "title": "Comparative signal",
        "summary": "Pagel lambda fitting and lambda-zero statistical checks that rely on governed tolerance rules.",
    },
    "ancestral-reconstruction": {
        "title": "Ancestral reconstruction",
        "summary": "Internal-node estimates, intervals, and derivative summaries tied to the Brownian comparative workflow.",
    },
    "artifact-provenance": {
        "title": "Artifact provenance",
        "summary": "Checked-in processed outputs and saved workspaces tracked as provenance, not overstated as direct numerical parity claims.",
    },
}

CLAIM_DEFINITIONS = {
    "pcm1-reproducibility-contract-tracked": {
        "claim_title": "PCM1 reproducibility contract is tracked explicitly",
        "summary": "The lecture setup and package context remain visible for reviewers without being misreported as a numerical parity claim.",
        "verdict": "not_comparable",
    },
    "pcm1-primate-data-preparation-parity": {
        "claim_title": "PCM1 primate data preparation matches the governed reference",
        "summary": "Workbook-derived preprocessing and grouped trait preparation agree with the checked-in R reference artifacts.",
        "verdict": "matched",
    },
    "pcm1-tree-import-pruning-parity": {
        "claim_title": "PCM1 tree import and pruning behavior matches the governed reference",
        "summary": "Tree loading, diagnostics, and pruning outcomes align with the R-derived trimmed tree and its taxon set.",
        "verdict": "matched",
    },
    "pcm1-topology-operation-parity": {
        "claim_title": "PCM1 topology teaching operations preserve the intended structure",
        "summary": "Unrooting, clade extraction, and node rotation are reproduced with matching taxon identity or tip-order outcomes.",
        "verdict": "matched",
    },
    "pcm1-tree-data-join-parity": {
        "claim_title": "PCM1 tree and trait joining logic matches across R and Bijux",
        "summary": "Tree tip order, aligned trait vectors, and representative joined node mappings stay consistent across both sides.",
        "verdict": "matched",
    },
    "pcm1-random-lambda-fit-parity": {
        "claim_title": "PCM1 random-data lambda fits agree within governed tolerance",
        "summary": "Seeded simulation inputs and downstream lambda fits stay within explicitly bounded numerical tolerance.",
        "verdict": "matched_with_tolerance",
    },
    "pcm1-primate-lambda-fit-parity": {
        "claim_title": "PCM1 primate lambda inference agrees within governed tolerance",
        "summary": "Real-data lambda estimation and lambda-zero likelihood-ratio checks remain numerically aligned within stated tolerance.",
        "verdict": "matched_with_tolerance",
    },
    "pcm1-ancestral-state-parity": {
        "claim_title": "PCM1 ancestral-state outputs agree within governed tolerance",
        "summary": "Internal-node estimates, confidence intervals, MRCA spot checks, and increase counts agree to floating-point noise or exact counts.",
        "verdict": "matched_with_tolerance",
    },
    "pcm1-visual-surface-tracking": {
        "claim_title": "PCM1 visual surfaces are tracked without false equivalence claims",
        "summary": "Plotting examples remain indexed and reviewer-visible while figure-equivalence claims stay intentionally out of scope.",
        "verdict": "not_comparable",
    },
    "pcm1-artifact-provenance-tracking": {
        "claim_title": "PCM1 saved artifacts remain governed provenance surfaces",
        "summary": "Processed files and saved workspaces are indexed as provenance outputs rather than overstated as analytical matches.",
        "verdict": "not_comparable",
    },
}

FRAGMENT_CLASSIFICATIONS = {
    "environment-and-package-contract": {
        "concept_family": "workflow-contracts",
        "claim_ids": ["pcm1-reproducibility-contract-tracked"],
        "parity_expectation": "not_comparable",
        "scope": "workflow",
    },
    "primate-data-preprocessing": {
        "concept_family": "data-preparation",
        "claim_ids": ["pcm1-primate-data-preparation-parity"],
        "parity_expectation": "exact",
        "scope": "analytical",
    },
    "tree-import-and-pruning": {
        "concept_family": "tree-operations",
        "claim_ids": ["pcm1-tree-import-pruning-parity"],
        "parity_expectation": "exact",
        "scope": "structural",
    },
    "processed-analysis-artifacts": {
        "concept_family": "artifact-provenance",
        "claim_ids": ["pcm1-artifact-provenance-tracking"],
        "parity_expectation": "not_comparable",
        "scope": "artifact",
    },
    "ape-plotting-basics": {
        "concept_family": "visual-surfaces",
        "claim_ids": ["pcm1-visual-surface-tracking"],
        "parity_expectation": "not_comparable",
        "scope": "visual",
    },
    "ape-alternate-layouts": {
        "concept_family": "visual-surfaces",
        "claim_ids": ["pcm1-visual-surface-tracking"],
        "parity_expectation": "not_comparable",
        "scope": "visual",
    },
    "unrooted-tree-demo": {
        "concept_family": "tree-operations",
        "claim_ids": ["pcm1-topology-operation-parity"],
        "parity_expectation": "exact",
        "scope": "structural",
    },
    "phytools-tree-plotting": {
        "concept_family": "visual-surfaces",
        "claim_ids": ["pcm1-visual-surface-tracking"],
        "parity_expectation": "not_comparable",
        "scope": "visual",
    },
    "extract-clade-node-77": {
        "concept_family": "tree-operations",
        "claim_ids": ["pcm1-topology-operation-parity"],
        "parity_expectation": "exact",
        "scope": "structural",
    },
    "rotate-nodes-behavior": {
        "concept_family": "tree-operations",
        "claim_ids": ["pcm1-topology-operation-parity"],
        "parity_expectation": "exact",
        "scope": "structural",
    },
    "ggtree-tree-visualization": {
        "concept_family": "visual-surfaces",
        "claim_ids": ["pcm1-visual-surface-tracking"],
        "parity_expectation": "not_comparable",
        "scope": "visual",
    },
    "tip-order-alignment": {
        "concept_family": "data-preparation",
        "claim_ids": ["pcm1-tree-data-join-parity"],
        "parity_expectation": "exact",
        "scope": "structural",
    },
    "ape-longevity-overlay": {
        "concept_family": "visual-surfaces",
        "claim_ids": ["pcm1-visual-surface-tracking"],
        "parity_expectation": "not_comparable",
        "scope": "visual",
    },
    "treeio-node-mapping-and-join": {
        "concept_family": "tree-operations",
        "claim_ids": ["pcm1-tree-data-join-parity"],
        "parity_expectation": "exact",
        "scope": "structural",
    },
    "joined-ggtree-trait-plotting": {
        "concept_family": "visual-surfaces",
        "claim_ids": ["pcm1-visual-surface-tracking"],
        "parity_expectation": "not_comparable",
        "scope": "visual",
    },
    "random-simulation-inputs": {
        "concept_family": "simulation-inputs",
        "claim_ids": ["pcm1-random-lambda-fit-parity"],
        "parity_expectation": "exact",
        "scope": "artifact",
    },
    "random-simulation-plotting": {
        "concept_family": "visual-surfaces",
        "claim_ids": ["pcm1-visual-surface-tracking"],
        "parity_expectation": "not_comparable",
        "scope": "visual",
    },
    "random-signal-lambda-fits": {
        "concept_family": "comparative-signal",
        "claim_ids": ["pcm1-random-lambda-fit-parity"],
        "parity_expectation": "statistical_tolerance",
        "scope": "analytical",
    },
    "primate-longevity-visual-inspection": {
        "concept_family": "visual-surfaces",
        "claim_ids": ["pcm1-visual-surface-tracking"],
        "parity_expectation": "not_comparable",
        "scope": "visual",
    },
    "primate-longevity-vector-assembly": {
        "concept_family": "data-preparation",
        "claim_ids": ["pcm1-tree-data-join-parity"],
        "parity_expectation": "exact",
        "scope": "analytical",
    },
    "primate-lambda-fit": {
        "concept_family": "comparative-signal",
        "claim_ids": ["pcm1-primate-lambda-fit-parity"],
        "parity_expectation": "statistical_tolerance",
        "scope": "analytical",
    },
    "lambda-zero-visual-comparison": {
        "concept_family": "visual-surfaces",
        "claim_ids": ["pcm1-visual-surface-tracking"],
        "parity_expectation": "not_comparable",
        "scope": "visual",
    },
    "lambda-zero-covariance-and-lrt": {
        "concept_family": "comparative-signal",
        "claim_ids": ["pcm1-primate-lambda-fit-parity"],
        "parity_expectation": "near_exact",
        "scope": "analytical",
    },
    "continuous-ancestral-point-estimates": {
        "concept_family": "ancestral-reconstruction",
        "claim_ids": ["pcm1-ancestral-state-parity"],
        "parity_expectation": "near_exact",
        "scope": "analytical",
    },
    "continuous-ancestral-intervals": {
        "concept_family": "ancestral-reconstruction",
        "claim_ids": ["pcm1-ancestral-state-parity"],
        "parity_expectation": "near_exact",
        "scope": "analytical",
    },
    "ancestral-table-assembly": {
        "concept_family": "ancestral-reconstruction",
        "claim_ids": ["pcm1-ancestral-state-parity"],
        "parity_expectation": "exact",
        "scope": "structural",
    },
    "ancestral-colored-tree-plot": {
        "concept_family": "visual-surfaces",
        "claim_ids": ["pcm1-visual-surface-tracking"],
        "parity_expectation": "not_comparable",
        "scope": "visual",
    },
    "bonobo-gibbon-mrca-estimate": {
        "concept_family": "ancestral-reconstruction",
        "claim_ids": ["pcm1-ancestral-state-parity"],
        "parity_expectation": "near_exact",
        "scope": "analytical",
    },
    "lifespan-increase-counts": {
        "concept_family": "ancestral-reconstruction",
        "claim_ids": ["pcm1-ancestral-state-parity"],
        "parity_expectation": "exact",
        "scope": "analytical",
    },
    "final-workspace-artifact": {
        "concept_family": "artifact-provenance",
        "claim_ids": ["pcm1-artifact-provenance-tracking"],
        "parity_expectation": "not_comparable",
        "scope": "artifact",
    },
}


def _study_root(repo_root: Path) -> Path:
    return Path(repo_root) / "evidence-book" / "studies" / STUDY_ID


def _bundle_root(repo_root: Path) -> Path:
    return _study_root(repo_root) / EVIDENCE_ID


def _comparison_payload(repo_root: Path) -> dict[str, object]:
    path = _bundle_root(repo_root) / "comparison.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _line_spans(spec: str) -> list[dict[str, int]]:
    spans: list[dict[str, int]] = []
    for raw_part in spec.split(","):
        part = raw_part.strip()
        if not part:
            continue
        if "-" in part:
            start_text, end_text = part.split("-", maxsplit=1)
            start = int(start_text)
            end = int(end_text)
        else:
            start = int(part)
            end = int(part)
        spans.append({"start_line": start, "end_line": end})
    return spans


def _line_locators(spec: str) -> list[str]:
    locators: list[str] = []
    for span in _line_spans(spec):
        if span["start_line"] == span["end_line"]:
            locators.append(f"{SOURCE_LOCATOR}#L{span['start_line']}")
        else:
            locators.append(
                f"{SOURCE_LOCATOR}#L{span['start_line']}-L{span['end_line']}"
            )
    return locators


def build_primate_source_fragment_map(repo_root: Path) -> dict[str, object]:
    comparison = _comparison_payload(repo_root)
    fragments: list[dict[str, object]] = []
    for row in comparison["comparisons"]:
        block_id = row["block_id"]
        classification = FRAGMENT_CLASSIFICATIONS[block_id]
        fragments.append(
            {
                "fragment_id": block_id,
                "fragment_title": row["title"],
                "script_line_spec": row["script_lines"],
                "script_line_spans": _line_spans(row["script_lines"]),
                "script_locators": _line_locators(row["script_lines"]),
                "evidence_id": EVIDENCE_ID,
                "concept_family": classification["concept_family"],
                "claim_ids": classification["claim_ids"],
                "parity_expectation": classification["parity_expectation"],
                "scope": classification["scope"],
                "block_status": row["status"],
                "review_note": row["note"],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": EVIDENCE_ID,
        "source_locator": SOURCE_LOCATOR,
        "fragment_count": len(fragments),
        "fragments": fragments,
    }


def build_primate_claim_registry(repo_root: Path) -> dict[str, object]:
    fragment_map = build_primate_source_fragment_map(repo_root)
    fragments = fragment_map["fragments"]
    claims: list[dict[str, object]] = []
    for claim_id, definition in CLAIM_DEFINITIONS.items():
        claim_fragments = [
            fragment["fragment_id"]
            for fragment in fragments
            if claim_id in fragment["claim_ids"]
        ]
        claims.append(
            {
                "claim_id": claim_id,
                "claim_title": definition["claim_title"],
                "summary": definition["summary"],
                "verdict": definition["verdict"],
                "evidence_ids": [EVIDENCE_ID],
                "source_fragments": claim_fragments,
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": EVIDENCE_ID,
        "claim_count": len(claims),
        "claims": claims,
    }


def _family_verdict(fragment_statuses: list[str]) -> str:
    if all(status in {"workflow_only", "artifact_only", "plot_only"} for status in fragment_statuses):
        return "not_comparable"
    if "verified_with_tolerance" in fragment_statuses:
        return "matched_with_tolerance"
    return "matched"


def build_primate_family_index(repo_root: Path) -> dict[str, object]:
    fragment_map = build_primate_source_fragment_map(repo_root)
    claim_registry = build_primate_claim_registry(repo_root)
    fragments = fragment_map["fragments"]
    claims = claim_registry["claims"]
    families: list[dict[str, object]] = []
    for family_id, definition in FAMILY_DEFINITIONS.items():
        family_fragments = [
            fragment for fragment in fragments if fragment["concept_family"] == family_id
        ]
        family_claim_ids = sorted(
            {
                claim_id
                for fragment in family_fragments
                for claim_id in fragment["claim_ids"]
            }
        )
        claim_titles = {
            claim["claim_id"]: claim["claim_title"]
            for claim in claims
            if claim["claim_id"] in family_claim_ids
        }
        families.append(
            {
                "family_id": family_id,
                "family_title": definition["title"],
                "summary": definition["summary"],
                "evidence_ids": [EVIDENCE_ID],
                "fragment_ids": [fragment["fragment_id"] for fragment in family_fragments],
                "claim_ids": family_claim_ids,
                "claim_titles": claim_titles,
                "fragment_count": len(family_fragments),
                "family_verdict": _family_verdict(
                    [fragment["block_status"] for fragment in family_fragments]
                ),
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "source_family_id": "pcm1-plots-signal",
        "source_family_title": "PCM1 plots signal evidence family",
        "evidence_id": EVIDENCE_ID,
        "family_count": len(families),
        "families": families,
    }
