from __future__ import annotations

from pathlib import Path

from ..primate_longevity_signal import STUDY_ID as PRIMATE_LONGEVITY_SIGNAL_STUDY_ID

STUDY_ID = PRIMATE_LONGEVITY_SIGNAL_STUDY_ID

REFERENCE_SCRIPT_PATH = (
    "evidence-book/studies/primate-longevity-signal/reference/"
    "primate_lifespan_signal_reference_r.R"
)
REFERENCE_BUNDLE_ID = "evidence-001"
REFERENCE_BUNDLE_ROOT = (
    Path("evidence-book") / "studies" / "primate-longevity-signal" / REFERENCE_BUNDLE_ID
)
COMPONENT_BUNDLE_DEFINITIONS = [
    {
        "evidence_id": "evidence-002",
        "report_filename": "workbook-loading-parity.json",
        "title": "Primate workbook loading parity bundle",
        "summary": "Governed evidence for the raw workbook intake boundary and the resulting processed table contract.",
        "claim_id": "pcm1-workbook-loading-parity",
        "claim_title": "PCM1 workbook loading contract is governed explicitly",
        "claim_summary": "The lecture workbook boundary is represented explicitly and the resulting processed table contract matches the checked-in governed reference artifact.",
        "verdict": "matched",
        "claim_tags": ["teaching", "parity", "data-preparation", "workbook-loading"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["data-preparation"],
        "limitations": [
            "The raw workbook remains an external teaching source; the Python-side review starts from the governed processed reference artifact."
        ],
        "source_fragments": ["primate-data-preprocessing"],
        "reference_line_specs": ["33-34", "77"],
        "source_basis": [
            {
                "kind": "external-course-material",
                "label": "lund pcm1 raw primate workbook",
                "locator": "external:lund/pcm1-plots-signal/data/primate_raw.xlsx",
            },
            {
                "kind": "repository-reference",
                "label": "governed processed primate table",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
            },
        ],
    },
    {
        "evidence_id": "evidence-003",
        "report_filename": "type-repair-parity.json",
        "title": "Primate type repair parity bundle",
        "summary": "Governed evidence for the factor and numeric coercion contract before grouped downstream analysis.",
        "claim_id": "pcm1-type-repair-parity",
        "claim_title": "PCM1 type repair contract is explicit and reviewable",
        "claim_summary": "The reference factor and numeric coercion rules are made explicit and the governed processed artifact satisfies the repaired column contract.",
        "verdict": "matched",
        "claim_tags": ["teaching", "parity", "data-preparation", "type-repair"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["data-preparation"],
        "limitations": [
            "The Python-side review validates the repaired processed artifact rather than replaying the external workbook import directly."
        ],
        "source_fragments": ["primate-data-preprocessing"],
        "reference_line_specs": ["39-40", "59-71"],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": "governed processed primate table",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
            }
        ],
    },
    {
        "evidence_id": "evidence-004",
        "report_filename": "missing-data-accounting-parity.json",
        "title": "Primate missing-data accounting parity bundle",
        "summary": "Governed evidence for how missing values are handled before the processed primate table is used downstream.",
        "claim_id": "pcm1-missing-data-accounting-parity",
        "claim_title": "PCM1 missing-data accounting is visible before downstream inference",
        "claim_summary": "Missing-value handling is surfaced explicitly at the preprocessing boundary instead of being hidden beneath later comparative analyses.",
        "verdict": "matched",
        "claim_tags": ["teaching", "parity", "data-preparation", "missing-data"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["data-preparation"],
        "limitations": [
            "The raw workbook remains external, so the bundle records the governed missing-value contract at the processed artifact boundary."
        ],
        "source_fragments": ["primate-data-preprocessing"],
        "reference_line_specs": ["41-47"],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": "governed processed primate table",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
            }
        ],
    },
    {
        "evidence_id": "evidence-005",
        "report_filename": "duplicate-species-aggregation-parity.json",
        "title": "Primate duplicate-species aggregation parity bundle",
        "summary": "Governed evidence for duplicate-species consolidation before processed trait analysis.",
        "claim_id": "pcm1-duplicate-species-aggregation-parity",
        "claim_title": "PCM1 duplicate-species aggregation is explicit and quantified",
        "claim_summary": "The raw-to-processed row contraction and the resulting unique-species contract are broken out as standalone evidence instead of being buried inside later model checks.",
        "verdict": "matched",
        "claim_tags": [
            "teaching",
            "parity",
            "data-preparation",
            "duplicate-aggregation",
        ],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["data-preparation"],
        "limitations": [
            "The raw workbook is external, so duplicate consolidation is demonstrated through the governed row-count and species-count contract."
        ],
        "source_fragments": ["primate-data-preprocessing"],
        "reference_line_specs": ["41-47"],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": "governed processed primate table",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
            }
        ],
    },
    {
        "evidence_id": "evidence-006",
        "report_filename": "tree-import-parity.json",
        "title": "Primate tree import parity bundle",
        "summary": "Governed evidence for reading the primate tree and preserving the expected original and trimmed tree contract.",
        "claim_id": "pcm1-tree-import-parity",
        "claim_title": "PCM1 tree import behavior is governed explicitly",
        "claim_summary": "Original and trimmed tree imports are surfaced as their own evidence boundary with explicit tip-count and branch-length expectations.",
        "verdict": "matched",
        "claim_tags": ["teaching", "parity", "tree-operations", "tree-import"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["tree-preprocessing"],
        "limitations": [
            "This bundle focuses on import and preserved tree shape; downstream correspondence and diagnostics are broken out separately."
        ],
        "source_fragments": ["tree-import-and-pruning"],
        "reference_line_specs": ["79-85"],
        "source_basis": [
            {
                "kind": "external-course-material",
                "label": "lund pcm1 original primate tree",
                "locator": "external:lund/pcm1-plots-signal/data/primatetree.nex",
            },
            {
                "kind": "repository-reference",
                "label": "governed trimmed tree export",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_trimmed_primatetree.nwk",
            },
        ],
    },
    {
        "evidence_id": "evidence-007",
        "report_filename": "tree-diagnostics-parity.json",
        "title": "Primate tree diagnostics parity bundle",
        "summary": "Governed evidence for rootedness, binary-tree status, ultrametricity, and node-label interpretation.",
        "claim_id": "pcm1-tree-diagnostics-parity",
        "claim_title": "PCM1 tree diagnostics are separated from later comparative results",
        "claim_summary": "Tree rootedness, binary structure, ultrametricity, and node-label interpretation are recorded as standalone trust surfaces.",
        "verdict": "matched",
        "claim_tags": ["teaching", "parity", "tree-operations", "tree-diagnostics"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["tree-validation"],
        "limitations": [
            "This bundle records the checked-in diagnostic contract; visual tree surfaces remain intentionally separate."
        ],
        "source_fragments": ["tree-import-and-pruning", "extract-clade-node-77"],
        "reference_line_specs": ["80-83", "87-95", "109"],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": "governed tree-processing comparison payload",
                "locator": "evidence-book/studies/primate-longevity-signal/evidence-001/results/block-payloads/tree-import-and-pruning.json",
            }
        ],
    },
    {
        "evidence_id": "evidence-008",
        "report_filename": "tree-data-correspondence-parity.json",
        "title": "Primate tree-data correspondence parity bundle",
        "summary": "Governed evidence for missing-tip trimming, tip-order alignment, and node-aware tree-data joining.",
        "claim_id": "pcm1-tree-data-correspondence-parity",
        "claim_title": "PCM1 tree-data correspondence is explicit before model fitting",
        "claim_summary": "Missing-tip handling, aligned species ordering, and joined node identity are isolated as prerequisite trust surfaces before downstream inference.",
        "verdict": "matched",
        "claim_tags": [
            "teaching",
            "parity",
            "data-preparation",
            "tree-data-correspondence",
        ],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["tree-data-alignment"],
        "limitations": [
            "This bundle covers correspondence and join readiness only; later lambda and ancestral-state results are separate evidence."
        ],
        "source_fragments": [
            "tree-import-and-pruning",
            "tip-order-alignment",
            "treeio-node-mapping-and-join",
        ],
        "reference_line_specs": ["81-83", "99-107"],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": "governed processed primate table",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
            },
            {
                "kind": "repository-reference",
                "label": "governed trimmed tree export",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_trimmed_primatetree.nwk",
            },
        ],
    },
    {
        "evidence_id": "evidence-009",
        "report_filename": "processed-export-parity.json",
        "title": "Primate processed export parity bundle",
        "summary": "Governed evidence for the exported processed CSV and trimmed tree artifacts consumed by later teaching steps.",
        "claim_id": "pcm1-processed-export-parity",
        "claim_title": "PCM1 processed exports are governed as explicit downstream inputs",
        "claim_summary": "The processed CSV and trimmed tree exports are treated as first-class governed artifacts rather than hidden scaffolding for later parity claims.",
        "verdict": "matched",
        "claim_tags": ["teaching", "parity", "artifact-provenance", "processed-export"],
        "comparison_mode": "direct_parity",
        "analytical_surfaces": ["artifact-provenance"],
        "limitations": [
            "The governed export contract records the downstream artifacts and their parity expectations, not the rendered plot outputs that consume them."
        ],
        "source_fragments": ["processed-analysis-artifacts"],
        "reference_line_specs": ["77", "85"],
        "source_basis": [
            {
                "kind": "repository-reference",
                "label": "governed processed primate table",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
            },
            {
                "kind": "repository-reference",
                "label": "governed trimmed tree export",
                "locator": "evidence-book/studies/primate-longevity-signal/datasets/reference_trimmed_primatetree.nwk",
            },
        ],
    },
]

NUMERIC_COLUMNS = [
    "body_mass",
    "gestation",
    "home_range",
    "longevity",
    "social_group_size",
]
CATEGORICAL_COLUMNS = ["sex_dimorphism", "mating_system"]
TEXT_COLUMNS = ["family", "species"]
DATA_PREPARATION_BUNDLE_IDS = [
    "evidence-002",
    "evidence-003",
    "evidence-004",
    "evidence-005",
    "evidence-008",
]
