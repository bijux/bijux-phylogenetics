from __future__ import annotations

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
