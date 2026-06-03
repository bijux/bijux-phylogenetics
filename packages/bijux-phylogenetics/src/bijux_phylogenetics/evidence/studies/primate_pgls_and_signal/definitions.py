from __future__ import annotations

from pathlib import Path

STUDY_ID = "primate-pgls-and-signal"
SUMMARY_EVIDENCE_ID = "evidence-001"
PCM2_SOURCE_LOCATOR = "external:lund/pcm2-modes-pgls/script"
PCM2_REFERENCE_SCRIPT_PATH = (
    "evidence-book/studies/primate-pgls-and-signal/reference/"
    "primate_pgls_and_signal_reference_r.R"
)
STUDY_ONE_REFERENCE_ROOT = (
    Path("evidence-book") / "studies" / "primate-longevity-signal" / "datasets"
)

FAMILY_DEFINITIONS = {
    "workflow-contracts": {
        "title": "Workflow contracts",
        "summary": "Lecture workspace assumptions are represented explicitly without being overstated as numerical parity.",
    },
    "transformed-tree-workflows": {
        "title": "Transformed tree workflows",
        "summary": "OU, early-burst, and late-burst branch rescaling outputs are compared before downstream model claims are judged.",
    },
    "continuous-model-fitting": {
        "title": "Continuous model fitting",
        "summary": "Lecture BM, OU, and early-burst fitContinuous-style intercept fits are checked with governed parameters and fit statistics.",
    },
    "likelihood-ratio-tests": {
        "title": "Likelihood-ratio tests",
        "summary": "BM, OU, and early-burst model-comparison statistics are checked explicitly instead of being inferred from prose.",
    },
    "ancestral-reconstruction": {
        "title": "Ancestral reconstruction",
        "summary": "Brownian and early-burst ancestral-state estimates are compared directly on governed node identities.",
    },
    "baseline-regression": {
        "title": "Baseline regression",
        "summary": "Non-phylogenetic regression outputs stay visible before phylogenetic correction is judged.",
    },
    "phylogenetic-regression": {
        "title": "Phylogenetic regression",
        "summary": "Pagel-lambda regression outputs are compared with explicit tolerance and conclusion rules.",
    },
    "phylogenetic-signal": {
        "title": "Phylogenetic signal",
        "summary": "Intercept-only signal testing and lambda-difference statistics are governed as their own trust surface.",
    },
    "diagnostics": {
        "title": "Diagnostics",
        "summary": "Residual, QQ, and heteroscedasticity checks are machine-recorded so visual diagnostics do not stay trapped in plots.",
    },
    "coverage-boundaries": {
        "title": "Coverage boundaries",
        "summary": "The remaining lecture intercept-mode likelihood sweep stays explicit until canonical runtime coverage exists for that exact boundary.",
    },
}

CLAIM_DEFINITIONS = {
    "pcm2-reload-contract-governed": {
        "claim_title": "PCM2 reload semantics are represented explicitly",
        "summary": "The lecture one-line workspace reload contract is reconstructed from governed repository artifacts instead of being left as an opaque external assumption.",
        "verdict": "matched",
    },
    "pcm2-transformed-tree-parity": {
        "claim_title": "PCM2 transformed-tree workflows match on governed branch summaries",
        "summary": "OU, early-burst, and late-burst rescaling stay reviewable through deterministic branch and total-length comparisons before downstream comparative claims are made.",
        "verdict": "matched",
    },
    "pcm2-fitcontinuous-parity": {
        "claim_title": "PCM2 fitContinuous-style evolutionary mode fits agree within governed tolerance",
        "summary": "Brownian, OU, and early-burst intercept fits preserve the same parameter and model-fit story while allowing bounded numerical drift from the R reference.",
        "verdict": "matched_with_tolerance",
    },
    "pcm2-likelihood-ratio-parity": {
        "claim_title": "PCM2 likelihood-ratio model comparisons agree within governed tolerance",
        "summary": "BM-versus-OU, BM-versus-early-burst, and OU-versus-early-burst test statistics stay aligned under explicit tolerance rules.",
        "verdict": "matched_with_tolerance",
    },
    "pcm2-ancestral-parity": {
        "claim_title": "PCM2 Brownian and early-burst ancestral reconstructions agree within governed tolerance",
        "summary": "Internal-node ancestral estimates stay comparable across the canonical runtime and governed R reference on named node identities.",
        "verdict": "matched_with_tolerance",
    },
    "pcm2-baseline-gls-parity": {
        "claim_title": "PCM2 baseline GLS outputs match before phylogenetic correction",
        "summary": "The non-phylogenetic regression slope, intercept, fit statistics, and coefficient decisions align before any phylogenetic covariance is introduced.",
        "verdict": "matched",
    },
    "pcm2-pagel-lambda-regression-parity": {
        "claim_title": "PCM2 Pagel-lambda regression agrees within governed tolerance",
        "summary": "Estimated-lambda regression outputs keep the same analytical conclusion while allowing bounded numerical drift between R and Bijux implementations.",
        "verdict": "matched_with_tolerance",
    },
    "pcm2-phylogenetic-signal-parity": {
        "claim_title": "PCM2 phylogenetic signal testing agrees within governed tolerance",
        "summary": "Lambda estimation and lambda-zero likelihood-ratio testing preserve the same signal conclusion while tolerating bounded numerical differences.",
        "verdict": "matched_with_tolerance",
    },
    "pcm2-diagnostics-parity": {
        "claim_title": "PCM2 residual diagnostics are machine-recorded and scientifically equivalent",
        "summary": "Residual diagnostics remain reviewer-visible as governed scalar summaries instead of being trapped in hand-inspected plots.",
        "verdict": "matched_with_tolerance",
    },
    "pcm2-coverage-boundary-explicit": {
        "claim_title": "PCM2 remaining intercept-mode likelihood sweep boundary stays explicit",
        "summary": "The lecture corBlomberg likelihood sweep is kept visible as a bounded trust surface instead of being silently implied by the new parity bundles.",
        "verdict": "not_comparable",
    },
}

FRAGMENT_DEFINITIONS = [
    {
        "fragment_id": "workspace-reload-contract",
        "fragment_title": "Package loading and one-line primate workspace reload",
        "family_id": "workflow-contracts",
        "claim_ids": ["pcm2-reload-contract-governed"],
        "evidence_id": "evidence-001",
        "supporting_evidence_ids": [],
        "script_line_spec": "8-16",
        "parity_expectation": "exact",
        "comparison_kind": "exact_answer",
        "block_status": "verified",
        "review_note": "The lecture `load()` contract is reconstructed from governed repository artifacts with explicit object and path semantics.",
        "scope": "workflow",
    },
    {
        "fragment_id": "transformed-tree-workflows",
        "fragment_title": "OU and EB tree rescaling exploration",
        "family_id": "transformed-tree-workflows",
        "claim_ids": ["pcm2-transformed-tree-parity"],
        "evidence_id": "evidence-006",
        "supporting_evidence_ids": [],
        "script_line_spec": "18-30",
        "parity_expectation": "exact",
        "comparison_kind": "exact_answer",
        "block_status": "verified",
        "review_note": "The lecture tree-rescaling surface is checked through deterministic transformed-branch summaries before fit statistics are compared.",
        "scope": "analytical",
    },
    {
        "fragment_id": "continuous-model-comparison",
        "fragment_title": "BM, OU, and EB fitContinuous model comparison",
        "family_id": "continuous-model-fitting",
        "claim_ids": ["pcm2-fitcontinuous-parity", "pcm2-likelihood-ratio-parity"],
        "evidence_id": "evidence-007",
        "supporting_evidence_ids": ["evidence-008"],
        "script_line_spec": "36-87",
        "parity_expectation": "statistical_tolerance",
        "comparison_kind": "tolerance_or_equivalence",
        "block_status": "verified",
        "review_note": "The lecture fitContinuous surfaces are checked through governed Brownian, OU, and early-burst intercept fits plus their explicit likelihood-ratio rows.",
        "scope": "analytical",
    },
    {
        "fragment_id": "ancestral-mode-comparison",
        "fragment_title": "Ancestral-state comparison under BM and EB",
        "family_id": "ancestral-reconstruction",
        "claim_ids": ["pcm2-ancestral-parity"],
        "evidence_id": "evidence-009",
        "supporting_evidence_ids": [],
        "script_line_spec": "89-111",
        "parity_expectation": "statistical_tolerance",
        "comparison_kind": "tolerance_or_equivalence",
        "block_status": "verified",
        "review_note": "The lecture ancestral comparison is matched on governed node identities using Brownian and early-burst reconstructions.",
        "scope": "analytical",
    },
    {
        "fragment_id": "baseline-gls-fit",
        "fragment_title": "Non-phylogenetic GLS fit",
        "family_id": "baseline-regression",
        "claim_ids": ["pcm2-baseline-gls-parity"],
        "evidence_id": "evidence-002",
        "supporting_evidence_ids": ["evidence-005"],
        "script_line_spec": "122-136",
        "parity_expectation": "exact",
        "comparison_kind": "exact_answer",
        "block_status": "verified",
        "review_note": "The baseline regression is checked before phylogenetic covariance enters the workflow.",
        "scope": "analytical",
    },
    {
        "fragment_id": "baseline-gls-diagnostics",
        "fragment_title": "Baseline GLS heteroscedasticity and QQ diagnostics",
        "family_id": "diagnostics",
        "claim_ids": ["pcm2-diagnostics-parity"],
        "evidence_id": "evidence-005",
        "supporting_evidence_ids": ["evidence-002"],
        "script_line_spec": "126-133",
        "parity_expectation": "scientific_equivalence",
        "comparison_kind": "scientific_equivalence",
        "block_status": "verified",
        "review_note": "Plot-only diagnostics are converted into machine-recorded scalar summaries with explicit equivalence rules.",
        "scope": "diagnostic",
    },
    {
        "fragment_id": "pagel-lambda-regression",
        "fragment_title": "Fixed-lambda and estimated-lambda PGLS fits",
        "family_id": "phylogenetic-regression",
        "claim_ids": ["pcm2-pagel-lambda-regression-parity"],
        "evidence_id": "evidence-003",
        "supporting_evidence_ids": ["evidence-005"],
        "script_line_spec": "138-179",
        "parity_expectation": "statistical_tolerance",
        "comparison_kind": "tolerance_or_equivalence",
        "block_status": "verified",
        "review_note": "The lecture fixed-lambda equivalence and estimated-lambda fit are compared with separate tolerance and conclusion rules.",
        "scope": "analytical",
    },
    {
        "fragment_id": "estimated-lambda-diagnostics",
        "fragment_title": "Estimated-lambda residual and fitted diagnostics",
        "family_id": "diagnostics",
        "claim_ids": ["pcm2-diagnostics-parity"],
        "evidence_id": "evidence-005",
        "supporting_evidence_ids": ["evidence-003"],
        "script_line_spec": "168-179",
        "parity_expectation": "scientific_equivalence",
        "comparison_kind": "scientific_equivalence",
        "block_status": "verified",
        "review_note": "The estimated-lambda residual checks are reduced to governed scalar diagnostics instead of staying as plot-only anecdotes.",
        "scope": "diagnostic",
    },
    {
        "fragment_id": "phylogenetic-signal-test",
        "fragment_title": "Intercept-only PGLS and lambda-zero likelihood-ratio testing",
        "family_id": "phylogenetic-signal",
        "claim_ids": ["pcm2-phylogenetic-signal-parity"],
        "evidence_id": "evidence-004",
        "supporting_evidence_ids": [],
        "script_line_spec": "181-192",
        "parity_expectation": "statistical_tolerance",
        "comparison_kind": "tolerance_or_equivalence",
        "block_status": "verified",
        "review_note": "Signal is judged by intercept-only likelihood surfaces and the lambda-zero model comparison, not by prose.",
        "scope": "analytical",
    },
    {
        "fragment_id": "evolutionary-mode-likelihood-ratios",
        "fragment_title": "BM, OU, and EB likelihood-ratio tests",
        "family_id": "likelihood-ratio-tests",
        "claim_ids": ["pcm2-likelihood-ratio-parity"],
        "evidence_id": "evidence-008",
        "supporting_evidence_ids": ["evidence-007"],
        "script_line_spec": "59-74",
        "parity_expectation": "statistical_tolerance",
        "comparison_kind": "tolerance_or_equivalence",
        "block_status": "verified",
        "review_note": "The model-comparison rows are judged separately from fitted parameter rows so the evidence-book can show where fit statistics and test statistics agree.",
        "scope": "analytical",
    },
    {
        "fragment_id": "mode-linked-intercept-models",
        "fragment_title": "Mode-linked intercept-only GLS surrogates for BM, OU, and EB",
        "family_id": "coverage-boundaries",
        "claim_ids": ["pcm2-coverage-boundary-explicit"],
        "evidence_id": "evidence-010",
        "supporting_evidence_ids": [],
        "script_line_spec": "194-227",
        "parity_expectation": "not_comparable",
        "comparison_kind": "not_comparable",
        "block_status": "coverage_gap",
        "review_note": "The lecture corBlomberg intercept sweep remains visible, but the canonical runtime does not yet expose parity for that exact likelihood-profile surface.",
        "scope": "analytical",
    },
]

BUNDLE_DEFINITIONS = [
    {
        "evidence_id": "evidence-001",
        "report_filename": "rdata-reload-semantics.json",
        "title": "Primate reload semantics bundle",
        "summary": "Governed representation of the lecture one-line primate workspace reload contract.",
        "claim_id": "pcm2-reload-contract-governed",
        "claim_tags": ["teaching", "parity", "workflow", "reload-semantics"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["workflow-contracts"],
        "source_fragments": ["workspace-reload-contract"],
        "limitations": [
            "The raw lecture `primate.RData` file remains external; this bundle governs the object contract using repository-owned reference artifacts."
        ],
    },
    {
        "evidence_id": "evidence-002",
        "report_filename": "baseline-gls-parity.json",
        "title": "Primate baseline GLS parity bundle",
        "summary": "Governed parity for the non-phylogenetic baseline regression before phylogenetic covariance enters the workflow.",
        "claim_id": "pcm2-baseline-gls-parity",
        "claim_tags": ["teaching", "parity", "baseline-regression", "gls"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["baseline-regression"],
        "source_fragments": ["baseline-gls-fit"],
        "limitations": [
            "This bundle covers the baseline regression only; phylogenetic-covariance behavior is isolated in a separate bundle."
        ],
    },
    {
        "evidence_id": "evidence-003",
        "report_filename": "pagel-lambda-regression-parity.json",
        "title": "Primate Pagel-lambda regression parity bundle",
        "summary": "Governed parity for fixed-lambda and estimated-lambda regression surfaces derived from the lecture workflow.",
        "claim_id": "pcm2-pagel-lambda-regression-parity",
        "claim_tags": ["teaching", "parity", "phylogenetic-regression", "pagel-lambda"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["phylogenetic-regression"],
        "source_fragments": ["pagel-lambda-regression"],
        "limitations": [
            "Likelihood and coefficient comparisons allow bounded numerical drift as long as the same scientific conclusion is preserved."
        ],
    },
    {
        "evidence_id": "evidence-004",
        "report_filename": "phylogenetic-signal-parity.json",
        "title": "Primate phylogenetic signal parity bundle",
        "summary": "Governed parity for intercept-only signal testing and lambda-zero likelihood-ratio logic.",
        "claim_id": "pcm2-phylogenetic-signal-parity",
        "claim_tags": ["teaching", "parity", "phylogenetic-signal", "lambda-zero"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["phylogenetic-signal"],
        "source_fragments": ["phylogenetic-signal-test"],
        "limitations": [
            "This bundle records the lecture’s intercept-only signal workflow, not the later BM/OU/EB surrogate-model section."
        ],
    },
    {
        "evidence_id": "evidence-005",
        "report_filename": "residual-diagnostics-parity.json",
        "title": "Primate residual diagnostics parity bundle",
        "summary": "Governed machine-readable residual diagnostics for the baseline and estimated-lambda regression surfaces.",
        "claim_id": "pcm2-diagnostics-parity",
        "claim_tags": ["teaching", "parity", "diagnostics", "residuals"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["diagnostics"],
        "source_fragments": [
            "baseline-gls-diagnostics",
            "estimated-lambda-diagnostics",
        ],
        "limitations": [
            "The lecture plots are summarized as scalar diagnostics; this bundle does not claim rendered-figure equivalence."
        ],
    },
    {
        "evidence_id": "evidence-006",
        "report_filename": "transformed-tree-parity.json",
        "title": "Primate transformed tree parity bundle",
        "summary": "Governed parity for the lecture OU, early-burst, and late-burst transformed-tree workflows.",
        "claim_id": "pcm2-transformed-tree-parity",
        "claim_tags": ["teaching", "parity", "transformed-tree", "evolutionary-modes"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["transformed-tree-workflows"],
        "source_fragments": ["transformed-tree-workflows"],
        "limitations": [
            "This bundle checks deterministic branch and total-length parity, not rendered figure equivalence."
        ],
    },
    {
        "evidence_id": "evidence-007",
        "report_filename": "continuous-mode-fit-parity.json",
        "title": "Primate evolutionary mode fit parity bundle",
        "summary": "Governed parity for the lecture Brownian, OU, and early-burst fitContinuous-style intercept fits.",
        "claim_id": "pcm2-fitcontinuous-parity",
        "claim_tags": ["teaching", "parity", "fitcontinuous", "evolutionary-modes"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["continuous-model-fitting"],
        "source_fragments": ["continuous-model-comparison"],
        "limitations": [
            "This bundle tracks intercept-only Brownian, OU, and early-burst fits; the later corBlomberg likelihood sweep remains separate."
        ],
    },
    {
        "evidence_id": "evidence-008",
        "report_filename": "likelihood-ratio-parity.json",
        "title": "Primate likelihood-ratio parity bundle",
        "summary": "Governed parity for the lecture Brownian, OU, and early-burst likelihood-ratio test logic.",
        "claim_id": "pcm2-likelihood-ratio-parity",
        "claim_tags": ["teaching", "parity", "likelihood-ratio", "evolutionary-modes"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["likelihood-ratio-tests"],
        "source_fragments": ["evolutionary-mode-likelihood-ratios"],
        "limitations": [
            "These rows judge the test statistics and p-values directly; they do not replace the underlying fitted-parameter bundle."
        ],
    },
    {
        "evidence_id": "evidence-009",
        "report_filename": "ancestral-mode-parity.json",
        "title": "Primate ancestral mode parity bundle",
        "summary": "Governed parity for the lecture Brownian and early-burst ancestral-state reconstruction comparison.",
        "claim_id": "pcm2-ancestral-parity",
        "claim_tags": ["teaching", "parity", "ancestral-reconstruction", "early-burst"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["ancestral-reconstruction"],
        "source_fragments": ["ancestral-mode-comparison"],
        "limitations": [
            "This bundle compares node-level ancestral estimates and does not claim rendered-node-size figure equivalence."
        ],
    },
    {
        "evidence_id": "evidence-010",
        "report_filename": "coverage-boundaries.json",
        "title": "Primate intercept sweep coverage boundary bundle",
        "summary": "Governed record of the remaining lecture corBlomberg intercept-mode likelihood sweep boundary.",
        "claim_id": "pcm2-coverage-boundary-explicit",
        "claim_tags": ["teaching", "coverage-gap", "not-comparable"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["coverage-boundaries"],
        "source_fragments": ["mode-linked-intercept-models"],
        "limitations": [
            "This bundle is intentionally a boundary register, not a parity claim."
        ],
    },
]
