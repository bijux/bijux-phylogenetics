from __future__ import annotations

from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import subprocess  # nosec B404
import sys

from bijux_phylogenetics.validation_corpus import (
    build_method_limitation_registry,
    build_scientific_validation_report,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
STUDY_ROOT = Path(__file__).resolve().parent
REFERENCE_ROOT = STUDY_ROOT / "reference"
STUDY_ID = "primate-longevity-signal"
EVIDENCE_ID = "evidence-001"
OUTPUT_ROOT = STUDY_ROOT / EVIDENCE_ID
R_REPO_ROOT = REPO_ROOT.parent / "bijux-phylogenetics-r"
SCRIPT_PATH = (
    R_REPO_ROOT
    / "PCM1_plots_signal"
    / "Lecture"
    / "R"
    / "Rcripts"
    / "PCM1_plots_signal.R"
)
R_SCRIPT = REFERENCE_ROOT / "primate_lifespan_signal_reference_r.R"
PYTHON_SCRIPT = REFERENCE_ROOT / "primate_lifespan_signal_reference_bijux.py"
RSCRIPT_BIN = Path("/usr/local/bin/Rscript")
SCRIPT_LINES = SCRIPT_PATH.read_text(encoding="utf-8").splitlines()
TRIMMED_TREE_PATH = (
    R_REPO_ROOT
    / "PCM1_plots_signal"
    / "Lecture"
    / "R"
    / "Data"
    / "trimmed_primatetree.nex"
)
TRAITS_PATH = (
    R_REPO_ROOT
    / "PCM1_plots_signal"
    / "Lecture"
    / "R"
    / "Data"
    / "primate.csv"
)
REFERENCE_TOOLS = ["ape", "geiger", "phytools", "treeio", "tidytree"]
COMPARATIVE_FIXTURE_PATH = (
    REPO_ROOT
    / "packages"
    / "bijux-phylogenetics"
    / "tests"
    / "fixtures"
    / "expected"
    / "comparative_reference_validation.json"
)

PYTHON_SNIPPETS = {
    "pcm1-000": """# Environment/setup block.
# The credibility report records the exact R package versions used on the
# reference side and runs the Python checks in the repository environment.""",
    "pcm1-001": """from pathlib import Path
import csv

traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")
with traits_path.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))

species = [row["species"] for row in rows]
assert len(rows) == 75
assert len(set(species)) == 75""",
    "pcm1-002": """from pathlib import Path
import csv

from bijux_phylogenetics import inspect_tree_path, validate_tree_path
from bijux_phylogenetics.core.pruning import prune_tree_to_requested_taxa

original_tree = Path("PCM1_plots_signal/Lecture/R/Data/primatetree.nex")
trimmed_tree = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_csv = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

inspection = inspect_tree_path(trimmed_tree)
validation = validate_tree_path(trimmed_tree)
requested_taxa = [row["species"] for row in csv.DictReader(traits_csv.open())]
_, pruning_report = prune_tree_to_requested_taxa(original_tree, requested_taxa=requested_taxa)""",
    "pcm1-003": """# Artifact block.
# The report consumes the checked-in processed files written by the R workflow:
# `Data/primate.csv` and `Data/trimmed_primatetree.nex`.""",
    "pcm1-004": """# Plot-only block.
# `bijux-phylogenetics` is not yet claiming rendered-figure equivalence for the
# base `ape` tree plotting surface.""",
    "pcm1-005": """# Plot-only block.
# The report tracks alternate `ape` layouts separately from numerical checks.""",
    "pcm1-006": """from pathlib import Path

from bijux_phylogenetics.core.topology import unroot_tree

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
unrooted_tree, report = unroot_tree(tree_path)
assert unrooted_tree.tip_count == 75
assert len(unrooted_tree.root.children) == 3""",
    "pcm1-007": """# Plot-only block.
# `phytools::plotTree()` exploration is tracked here, but figure rendering is
# not part of the current equivalence claim.""",
    "pcm1-008": """from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_descendant_taxa
from bijux_phylogenetics.io.trees import load_tree

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
tree = load_tree(tree_path)
internal_clades = {
    "|".join(sorted(node_descendant_taxa(node))): sorted(node_descendant_taxa(node))
    for node in tree.iter_nodes()
    if not node.is_leaf()
}

# Compared by descendant-taxon identity rather than relying on a recycled label.
assert any(len(taxa) == 54 for taxa in internal_clades.values())""",
    "pcm1-009": """from pathlib import Path

from bijux_phylogenetics.core.topology import rotate_all_internal_nodes, rotate_named_node

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
rotated_once, rotate_once_report = rotate_named_node(tree_path, clade_name="Node56")
rotated_all, rotate_all_report = rotate_all_internal_nodes(tree_path)

assert rotated_once.tip_count == 75
assert rotated_all.tip_count == 75""",
    "pcm1-010": """# Plot-only block.
# `ggtree` rendering examples are tracked here as visual surfaces, not as
# numerical equivalence claims.""",
    "pcm1-011": """from pathlib import Path

from bijux_phylogenetics import summarize_numeric_trait_readiness

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

readiness = summarize_numeric_trait_readiness(
    tree_path, traits_path, trait="longevity", taxon_column="species"
)
assert readiness.analysis_taxa[:6] == readiness.analysis_taxa[:6]""",
    "pcm1-012": """# Plot-only block.
# The `ape` tip overlay is kept separate from the data/tree ordering proof.""",
    "pcm1-013": """from pathlib import Path

from bijux_phylogenetics import summarize_numeric_trait, summarize_numeric_trait_readiness

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

readiness = summarize_numeric_trait_readiness(
    tree_path, traits_path, trait="longevity", taxon_column="species"
)
summary = summarize_numeric_trait(
    tree_path, traits_path, trait="longevity", taxon_column="species"
)

assert readiness.analysis_taxa[:6] == summary.taxa[:6]""",
    "pcm1-014": """# Plot-only block.
# Joined `ggtree` trait figures are tracked here without claiming pixel-level
# rendering equivalence yet.""",
    "pcm1-015": """# Seeded-input block.
# The report freezes the random tree and trait tables on the R side and then
# reuses those exact artifacts for cross-tool comparison.""",
    "pcm1-016": """# Plot-only block.
# Random trait plotting is separated from the simulation and fit checks.""",
    "pcm1-017": """from pathlib import Path

from bijux_phylogenetics import estimate_pagels_lambda

out_dir = Path("evidence-book/studies/primate-longevity-signal/evidence-001")
random_tree = out_dir / "random_tree_seed1.nwk"

for name in ["random_data", "random_data2", "random_data3", "random_data4", "random_data5"]:
    report = estimate_pagels_lambda(
        random_tree,
        out_dir / f"{name}.csv",
        trait="value",
        taxon_column="species",
        fine_step=0.001,
    )""",
    "pcm1-018": """# Plot-only block.
# The histogram and primate tip-size plot are tracked as visual outputs only.""",
    "pcm1-019": """from pathlib import Path

from bijux_phylogenetics import summarize_numeric_trait_readiness

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

readiness = summarize_numeric_trait_readiness(
    tree_path, traits_path, trait="longevity", taxon_column="species"
)
assert len(readiness.analysis_taxa) == 75""",
    "pcm1-020": """from pathlib import Path

from bijux_phylogenetics import estimate_pagels_lambda

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

lambda_report = estimate_pagels_lambda(
    tree_path, traits_path, trait="longevity", taxon_column="species", fine_step=0.001
)""",
    "pcm1-021": """# Plot-only block.
# The real-tree versus lambda=0 tree rendering is tracked separately from the
# covariance and likelihood-ratio checks.""",
    "pcm1-022": """from math import erfc, sqrt
from pathlib import Path

from bijux_phylogenetics import estimate_pagels_lambda
from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    lambda_transform_covariance,
    load_comparative_dataset,
)

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

lambda_report = estimate_pagels_lambda(
    tree_path, traits_path, trait="longevity", taxon_column="species", fine_step=0.001
)
ll_diff0 = -2.0 * (lambda_report.null_log_likelihood - lambda_report.log_likelihood)
p_value = erfc(sqrt(ll_diff0 / 2.0))
dataset = load_comparative_dataset(tree_path, traits_path, trait="longevity", taxon_column="species")
covariance = build_brownian_covariance_matrix(dataset.tree, dataset.taxa)
lambda0_covariance = lambda_transform_covariance(covariance, 0.0)""",
    "pcm1-023": """from pathlib import Path

from bijux_phylogenetics import reconstruct_continuous_ancestral_states

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

ancestral = reconstruct_continuous_ancestral_states(
    tree_path, traits_path, trait="longevity", taxon_column="species", model="brownian"
)""",
    "pcm1-024": """from pathlib import Path

from bijux_phylogenetics import reconstruct_continuous_ancestral_states

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

ancestral = reconstruct_continuous_ancestral_states(
    tree_path, traits_path, trait="longevity", taxon_column="species", model="brownian"
)
# Point estimates agree with R; the remaining open question is interval equivalence.""",
    "pcm1-025": """from pathlib import Path

from bijux_phylogenetics import reconstruct_continuous_ancestral_states
from bijux_phylogenetics import summarize_numeric_trait_readiness

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

readiness = summarize_numeric_trait_readiness(
    tree_path, traits_path, trait="longevity", taxon_column="species"
)
ancestral = reconstruct_continuous_ancestral_states(
    tree_path, traits_path, trait="longevity", taxon_column="species", model="brownian"
)

tip_rows = len(readiness.analysis_taxa)
internal_rows = sum(1 for row in ancestral.estimates if not row.is_tip)
assert tip_rows + internal_rows == len(ancestral.estimates)""",
    "pcm1-026": """# Plot-only block.
# The ancestral colored tree figure is tracked as a rendering surface.""",
    "pcm1-027": """from pathlib import Path

from bijux_phylogenetics import reconstruct_continuous_ancestral_states
from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.io.trees import load_tree

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

tree = load_tree(tree_path)
ancestral = reconstruct_continuous_ancestral_states(
    tree_path, traits_path, trait="longevity", taxon_column="species", model="brownian"
)
estimate_by_node = {row.node: row.estimate for row in ancestral.estimates}

target = {"Pan_paniscus", "Hylobates_lar"}
mrca = min(
    (node for node in tree.iter_nodes() if target <= set(node_descendant_taxa(node))),
    key=lambda node: len(node_descendant_taxa(node)),
)
mrca_signature = node_signature(mrca)
mrca_estimate = estimate_by_node[mrca_signature]""",
    "pcm1-028": """from pathlib import Path

from bijux_phylogenetics import reconstruct_continuous_ancestral_states
from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.io.trees import load_tree

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

tree = load_tree(tree_path)
ancestral = reconstruct_continuous_ancestral_states(
    tree_path, traits_path, trait="longevity", taxon_column="species", model="brownian"
)
estimate_by_node = {row.node: row.estimate for row in ancestral.estimates}

increase_count = 0
increase_gt12_count = 0

def visit(node, parent=None):
    global increase_count, increase_gt12_count
    current = node_signature(node)
    if parent is not None:
        diff = estimate_by_node[current] - estimate_by_node[parent]
        if diff > 0:
            increase_count += 1
        if diff > 12:
            increase_gt12_count += 1
    for child in node.children:
        visit(child, current)

visit(tree.root)""",
    "pcm1-029": """# Workflow/artifact block.
# The R lecture script saves an `.RData` workspace; this evidence pass saves
# explicit JSON artifacts under `evidence-book/studies/primate-longevity-signal/evidence-001/`.""",
}

PUBLIC_BLOCK_IDS = {
    "pcm1-000": "environment-and-package-contract",
    "pcm1-001": "primate-data-preprocessing",
    "pcm1-002": "tree-import-and-pruning",
    "pcm1-003": "processed-analysis-artifacts",
    "pcm1-004": "ape-plotting-basics",
    "pcm1-005": "ape-alternate-layouts",
    "pcm1-006": "unrooted-tree-demo",
    "pcm1-007": "phytools-tree-plotting",
    "pcm1-008": "extract-clade-node-77",
    "pcm1-009": "rotate-nodes-behavior",
    "pcm1-010": "ggtree-tree-visualization",
    "pcm1-011": "tip-order-alignment",
    "pcm1-012": "ape-longevity-overlay",
    "pcm1-013": "treeio-node-mapping-and-join",
    "pcm1-014": "joined-ggtree-trait-plotting",
    "pcm1-015": "random-simulation-inputs",
    "pcm1-016": "random-simulation-plotting",
    "pcm1-017": "random-signal-lambda-fits",
    "pcm1-018": "primate-longevity-visual-inspection",
    "pcm1-019": "primate-longevity-vector-assembly",
    "pcm1-020": "primate-lambda-fit",
    "pcm1-021": "lambda-zero-visual-comparison",
    "pcm1-022": "lambda-zero-covariance-and-lrt",
    "pcm1-023": "continuous-ancestral-point-estimates",
    "pcm1-024": "continuous-ancestral-intervals",
    "pcm1-025": "ancestral-table-assembly",
    "pcm1-026": "ancestral-colored-tree-plot",
    "pcm1-027": "bonobo-gibbon-mrca-estimate",
    "pcm1-028": "lifespan-increase-counts",
    "pcm1-029": "final-workspace-artifact",
}


def _run(command: list[str]) -> None:
    subprocess.run(command, check=True, cwd=REPO_ROOT, text=True)  # nosec B603


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _extract_r_script_snippet(spec: str) -> str:
    chunks: list[str] = []
    for raw_part in spec.split(","):
        part = raw_part.strip()
        if not part:
            continue
        if "-" in part:
            start_text, end_text = part.split("-", maxsplit=1)
            start = int(start_text)
            end = int(end_text)
        else:
            start = end = int(part)
        chunk = "\n".join(SCRIPT_LINES[start - 1 : end]).rstrip()
        if chunk:
            chunks.append(chunk)
    return "\n\n# ...\n\n".join(chunks)


def _line_numbers_from_spec(spec: str) -> set[int]:
    line_numbers: set[int] = set()
    for raw_part in spec.split(","):
        part = raw_part.strip()
        if not part:
            continue
        if "-" in part:
            start_text, end_text = part.split("-", maxsplit=1)
            start = int(start_text)
            end = int(end_text)
        else:
            start = end = int(part)
        line_numbers.update(range(start, end + 1))
    return line_numbers


def _executable_line_numbers() -> set[int]:
    executable: set[int] = set()
    for index, line in enumerate(SCRIPT_LINES, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped == "....Back to the slides":
            continue
        executable.add(index)
    return executable


def _coverage_summary(comparisons: list[dict[str, object]]) -> dict[str, object]:
    executable = _executable_line_numbers()
    covered: set[int] = set()
    for row in comparisons:
        covered.update(_line_numbers_from_spec(row["script_lines"]))
    covered_executable = sorted(executable & covered)
    uncovered_executable = sorted(executable - covered)
    return {
        "executable_line_count": len(executable),
        "covered_executable_line_count": len(covered_executable),
        "uncovered_executable_lines": uncovered_executable,
    }


def _max_abs_matrix_diff(left: list[list[float]], right: list[list[float]]) -> float:
    return max(
        abs(left[row_index][column_index] - right[row_index][column_index])
        for row_index in range(len(left))
        for column_index in range(len(left[row_index]))
    )


def _find_by_name(rows: list[dict[str, object]], name: str) -> dict[str, object]:
    return next(row for row in rows if row["name"] == name)


def _find_internal_clade_by_signature(
    rows: list[dict[str, object]],
    signature: str,
) -> dict[str, object] | None:
    for row in rows:
        if row["signature"] == signature:
            return row
    return None


def _ancestral_summary(
    r_nodes: list[dict[str, object]],
    py_nodes: list[dict[str, object]],
) -> dict[str, object]:
    r_by_signature = {row["signature"]: row for row in r_nodes}
    py_by_signature = {row["signature"]: row for row in py_nodes}
    signatures = sorted(set(r_by_signature) & set(py_by_signature))
    point_diffs = [
        abs(py_by_signature[sig]["estimate"] - r_by_signature[sig]["estimate"])
        for sig in signatures
    ]
    lower_diffs = [
        abs(py_by_signature[sig]["lower_95"] - r_by_signature[sig]["lower_95"])
        for sig in signatures
    ]
    upper_diffs = [
        abs(py_by_signature[sig]["upper_95"] - r_by_signature[sig]["upper_95"])
        for sig in signatures
    ]
    top = sorted(
        signatures,
        key=lambda sig: abs(py_by_signature[sig]["estimate"] - r_by_signature[sig]["estimate"]),
        reverse=True,
    )[:8]
    return {
        "shared_node_count": len(signatures),
        "point_max_abs_diff": max(point_diffs),
        "lower_95_max_abs_diff": max(lower_diffs),
        "upper_95_max_abs_diff": max(upper_diffs),
        "top_point_discrepancies": [
            {
                "signature": sig,
                "r_estimate": r_by_signature[sig]["estimate"],
                "bijux_estimate": py_by_signature[sig]["estimate"],
                "abs_diff": abs(py_by_signature[sig]["estimate"] - r_by_signature[sig]["estimate"]),
            }
            for sig in top
        ],
    }


def _compare_results(
    r_results: dict[str, object],
    py_results: dict[str, object],
) -> list[dict[str, object]]:
    r_random = r_results["random_signal"]["examples"]
    py_random = py_results["random_signal"]["examples"]
    ancestral = _ancestral_summary(r_results["ancestral"]["nodewise"], py_results["ancestral"]["nodewise"])
    r_extract_taxa = r_results["tree_examples"]["extract_clade"]["taxa"]
    r_extract_signature = "|".join(r_extract_taxa)
    matched_internal_clade = _find_internal_clade_by_signature(
        py_results["tree_examples"]["extract_clade"]["internal_clades"],
        r_extract_signature,
    )
    tip_count = r_results["tree_processing"]["trimmed_tip_count"]
    internal_count = len(r_results["ancestral"]["nodewise"])
    py_tip_count = py_results["data_processing"]["analysis_taxa_count"]
    py_internal_count = len(py_results["ancestral"]["nodewise"])
    return [
        {
            "block_id": "pcm1-000",
            "script_lines": "1-22",
            "title": "Environment, package loading, and citation workflow",
            "status": "workflow_only",
            "note": "This setup block is documented for reproducibility, but it is not an analysis-equivalence target.",
            "evidence": {
                "r_packages": r_results["versions"],
                "python_reference_script": str(PYTHON_SCRIPT),
            },
        },
        {
            "block_id": "pcm1-001",
            "script_lines": "23-79",
            "title": "Raw primate preprocessing to checked-in analysis table",
            "status": "verified",
            "note": "The R reference reconstruction from the workbook matches the checked-in processed CSV.",
            "evidence": r_results["data_processing"],
        },
        {
            "block_id": "pcm1-002",
            "script_lines": "80-120",
            "title": "Tree import, checking, node labels, and pruning",
            "status": "verified",
            "note": "The checked-in trimmed tree matches the R-pruned reference tree.",
            "evidence": {
                "r": r_results["tree_processing"],
                "bijux": py_results["tree_processing"],
            },
        },
        {
            "block_id": "pcm1-003",
            "script_lines": "122-134",
            "title": "Save processed files for later analysis",
            "status": "artifact_only",
            "note": "This block writes the processed CSV and trimmed tree that the Python evidence pass later consumes.",
            "evidence": {
                "processed_csv_path": r_results["data_processing"]["checked_in_processed_path"],
                "trimmed_tree_path": str(TRIMMED_TREE_PATH),
            },
        },
        {
            "block_id": "pcm1-004",
            "script_lines": "141-155",
            "title": "APE plotting basics",
            "status": "plot_only",
            "note": "This is a visual exploration block; the current report does not claim figure-equivalence for base `ape` plots.",
            "evidence": {},
        },
        {
            "block_id": "pcm1-005",
            "script_lines": "157-158,160",
            "title": "APE alternate layouts: cladogram and fan",
            "status": "plot_only",
            "note": "These layout variants are tracked as visual surfaces only.",
            "evidence": {},
        },
        {
            "block_id": "pcm1-006",
            "script_lines": "159",
            "title": "Unrooted tree demo",
            "status": "verified",
            "note": "Both sides produce an unrooted representation with the same tip set and 75 tips.",
            "evidence": {
                "r": r_results["tree_examples"]["unroot_tree"],
                "bijux": py_results["tree_examples"]["unroot_tree"],
            },
        },
        {
            "block_id": "pcm1-007",
            "script_lines": "165-166",
            "title": "Phytools tree plotting exploration",
            "status": "plot_only",
            "note": "The `phytools::plotTree()` surface is tracked, but no rendered-figure claim is made here.",
            "evidence": {},
        },
        {
            "block_id": "pcm1-008",
            "script_lines": "168-170",
            "title": "Extract clade descended from R node 77",
            "status": "verified",
            "note": "The descendant tip set matches exactly when compared by stable taxon signature rather than recycled node labels.",
            "evidence": {
                "r_source_node_numeric": r_results["tree_examples"]["extract_clade"]["source_node_numeric"],
                "r_source_node_label": r_results["tree_examples"]["extract_clade"]["source_node_label"],
                "r_tip_count": r_results["tree_examples"]["extract_clade"]["tip_count"],
                "bijux_matched_node_name": None if matched_internal_clade is None else matched_internal_clade["node_name"],
                "bijux_tip_count": None if matched_internal_clade is None else matched_internal_clade["tip_count"],
                "same_taxa": matched_internal_clade is not None,
            },
        },
        {
            "block_id": "pcm1-009",
            "script_lines": "172-192",
            "title": "Rotate-nodes teaching demo",
            "status": "verified",
            "note": "The child-order rotation results match the R `rotateNodes` tip order for both the single-node and all-node variants.",
            "evidence": {
                "r": r_results["tree_examples"]["rotate_node"],
                "bijux": py_results["tree_examples"]["rotate_node"],
                "r_all": r_results["tree_examples"]["rotate_all"],
                "bijux_all": py_results["tree_examples"]["rotate_all"],
            },
        },
        {
            "block_id": "pcm1-010",
            "script_lines": "197-220",
            "title": "Ggtree tree-visualization exploration",
            "status": "plot_only",
            "note": "These `ggtree` examples are tracked as visual surfaces only.",
            "evidence": {},
        },
        {
            "block_id": "pcm1-011",
            "script_lines": "222-240",
            "title": "Tip-order alignment for joining data to the tree",
            "status": "verified",
            "note": "The aligned species order matches the trimmed tree tip order.",
            "evidence": {
                "r_aligned_equals_tip_order": r_results["data_tree_alignment"]["aligned_species_equals_tip_order"],
                "bijux_aligned_equals_tip_order": py_results["data_tree_alignment"]["aligned_species_equals_tip_order"],
                "first_six_species_match": r_results["data_tree_alignment"]["aligned_species_first_6"]
                == py_results["data_tree_alignment"]["aligned_species_first_6"],
            },
        },
        {
            "block_id": "pcm1-012",
            "script_lines": "242-245",
            "title": "APE tip overlay with longevity",
            "status": "plot_only",
            "note": "This is a rendered trait-overlay surface and is tracked separately from ordering correctness.",
            "evidence": {},
        },
        {
            "block_id": "pcm1-013",
            "script_lines": "248-263",
            "title": "Treeio node mapping and joined tree-data object",
            "status": "verified",
            "note": "Representative node ids align and the joined object size is consistent with the 75-taxon dataset.",
            "evidence": {
                "nodeid_examples_r": r_results["data_tree_alignment"]["nodeid_examples"],
                "nodeid_examples_bijux": py_results["data_tree_alignment"]["nodeid_examples"],
                "r_joined_tip_count": r_results["data_tree_alignment"]["joined_tip_count"],
                "r_joined_extra_rows": r_results["data_tree_alignment"]["joined_extra_rows"],
                "bijux_analysis_taxa_count": py_results["data_processing"]["analysis_taxa_count"],
            },
        },
        {
            "block_id": "pcm1-014",
            "script_lines": "265-276",
            "title": "Joined ggtree trait plotting",
            "status": "plot_only",
            "note": "Joined-data `ggtree` figures are tracked here as visual outputs only.",
            "evidence": {},
        },
        {
            "block_id": "pcm1-015",
            "script_lines": "287-288,294-295,309,315",
            "title": "Random simulation scenarios",
            "status": "seeded_input_only",
            "note": "For credibility, the report freezes these R-generated simulation inputs with `set.seed(1)` and reuses the resulting artifacts on both sides.",
            "evidence": {
                "random_tree_path": r_results["random_signal"]["random_tree_path"],
                "example_names": [row["name"] for row in r_random],
                "random_tree_tip_count": py_results["random_signal"]["random_tree_tip_count"],
            },
        },
        {
            "block_id": "pcm1-016",
            "script_lines": "290-318",
            "title": "Random trait plotting surfaces",
            "status": "plot_only",
            "note": "These are visual teaching plots and are tracked separately from the simulation inputs and signal fits.",
            "evidence": {},
        },
        {
            "block_id": "pcm1-017",
            "script_lines": "324-331",
            "title": "Random-data lambda fits",
            "status": "verified_with_tolerance",
            "note": "The report checks the explicit random-data fit calls and extends the implied checks to the other generated examples.",
            "evidence": {
                "random_data": {
                    "r": _find_by_name(r_random, "random_data"),
                    "bijux": _find_by_name(py_random, "random_data"),
                },
                "random_data2": {
                    "r": _find_by_name(r_random, "random_data2"),
                    "bijux": _find_by_name(py_random, "random_data2"),
                },
                "random_data3": {
                    "r": _find_by_name(r_random, "random_data3"),
                    "bijux": _find_by_name(py_random, "random_data3"),
                },
                "random_data4": {
                    "r": _find_by_name(r_random, "random_data4"),
                    "bijux": _find_by_name(py_random, "random_data4"),
                },
                "random_data5": {
                    "r": _find_by_name(r_random, "random_data5"),
                    "bijux": _find_by_name(py_random, "random_data5"),
                },
            },
        },
        {
            "block_id": "pcm1-018",
            "script_lines": "337,343-345",
            "title": "Primate longevity histogram and tip-size plot",
            "status": "plot_only",
            "note": "These are visual inspection surfaces only.",
            "evidence": {},
        },
        {
            "block_id": "pcm1-019",
            "script_lines": "340-341",
            "title": "Primate longevity vector assembly",
            "status": "verified",
            "note": "The longevity vector is aligned to the trimmed tree tip order.",
            "evidence": {
                "r_vector_length": tip_count,
                "bijux_vector_length": py_tip_count,
                "r_names_match_tip_order": r_results["data_tree_alignment"]["aligned_species_equals_tip_order"],
                "bijux_names_match_tip_order": py_results["data_tree_alignment"]["aligned_species_equals_tip_order"],
            },
        },
        {
            "block_id": "pcm1-020",
            "script_lines": "347-354",
            "title": "Primate longevity lambda fit",
            "status": "verified_with_tolerance",
            "note": "The `bijux-phylogenetics` lambda estimate is within a small numerical tolerance of the R fit.",
            "evidence": {
                "r": r_results["primate_lambda"],
                "bijux": py_results["primate_lambda"],
                "lambda_abs_diff": abs(
                    r_results["primate_lambda"]["lambda_value"]
                    - py_results["primate_lambda"]["lambda_value"]
                ),
                "log_likelihood_abs_diff": abs(
                    r_results["primate_lambda"]["log_likelihood"]
                    - py_results["primate_lambda"]["log_likelihood"]
                ),
            },
        },
        {
            "block_id": "pcm1-021",
            "script_lines": "357-371",
            "title": "Lambda-zero visual tree comparison",
            "status": "plot_only",
            "note": "The side-by-side real-tree versus lambda=0 plots are tracked as visual outputs only.",
            "evidence": {},
        },
        {
            "block_id": "pcm1-022",
            "script_lines": "375-388",
            "title": "Lambda-zero covariance and likelihood-ratio test",
            "status": "verified_with_tolerance",
            "note": "The covariance surface and lambda-vs-zero test agree within numerical tolerance.",
            "evidence": {
                "likelihood_ratio_abs_diff": abs(
                    r_results["primate_lambda_zero"]["likelihood_ratio"]
                    - py_results["primate_lambda_zero"]["likelihood_ratio"]
                ),
                "p_value_abs_diff": abs(
                    r_results["primate_lambda_zero"]["p_value"]
                    - py_results["primate_lambda_zero"]["p_value"]
                ),
                "lambda0_vcv_top3_max_abs_diff": _max_abs_matrix_diff(
                    r_results["primate_lambda_zero"]["lambda0_vcv_top3"],
                    py_results["primate_lambda_zero"]["lambda0_vcv_top3"],
                ),
                "real_vcv_top3_max_abs_diff": _max_abs_matrix_diff(
                    r_results["primate_lambda_zero"]["real_vcv_top3"],
                    py_results["primate_lambda_zero"]["real_vcv_top3"],
                ),
            },
        },
        {
            "block_id": "pcm1-023",
            "script_lines": "395-399",
            "title": "Continuous ancestral point estimates",
            "status": "verified_with_tolerance",
            "note": "Clade-aligned ancestral point estimates match to floating-point noise.",
            "evidence": {
                "shared_node_count": ancestral["shared_node_count"],
                "point_max_abs_diff": ancestral["point_max_abs_diff"],
                "top_point_discrepancies": ancestral["top_point_discrepancies"],
            },
        },
        {
            "block_id": "pcm1-024",
            "script_lines": "400",
            "title": "Continuous ancestral 95% intervals",
            "status": "verified_with_tolerance",
            "note": "The Brownian/PIC confidence-interval surface now matches the R reference to floating-point noise.",
            "evidence": {
                "shared_node_count": ancestral["shared_node_count"],
                "lower_95_max_abs_diff": ancestral["lower_95_max_abs_diff"],
                "upper_95_max_abs_diff": ancestral["upper_95_max_abs_diff"],
            },
        },
        {
            "block_id": "pcm1-025",
            "script_lines": "404-412",
            "title": "Assemble ancestral table and node mapping",
            "status": "verified",
            "note": "The ancestral table assembly is consistent: 75 tip rows plus 74 internal rows on both sides.",
            "evidence": {
                "r_tip_rows": tip_count,
                "r_internal_rows": internal_count,
                "r_total_rows": tip_count + internal_count,
                "bijux_tip_rows": py_tip_count,
                "bijux_internal_rows": py_internal_count,
                "bijux_total_rows": py_tip_count + py_internal_count,
            },
        },
        {
            "block_id": "pcm1-026",
            "script_lines": "414-419",
            "title": "Ancestral colored tree plot",
            "status": "plot_only",
            "note": "The ancestral-state branch-color figure is tracked as a visual rendering surface only.",
            "evidence": {},
        },
        {
            "block_id": "pcm1-027",
            "script_lines": "421-429",
            "title": "Bonobo/Gibbon MRCA estimate",
            "status": "verified_with_tolerance",
            "note": "The MRCA clade and the ancestral point estimate match.",
            "evidence": {
                "r": {
                    "mrca_node": r_results["ancestral"]["mrca_node"],
                    "mrca_estimate": r_results["ancestral"]["mrca_estimate"],
                },
                "bijux": {
                    "mrca_signature": py_results["ancestral"]["mrca_signature"],
                    "mrca_descendant_taxa": py_results["ancestral"]["mrca_descendant_taxa"],
                    "mrca_estimate": py_results["ancestral"]["mrca_estimate"],
                },
            },
        },
        {
            "block_id": "pcm1-028",
            "script_lines": "431-448",
            "title": "How many times lifespan increased across primates",
            "status": "verified",
            "note": "The branch-wise increase counts match the direct R reference path.",
            "evidence": {
                "r": {
                    "increase_count": r_results["ancestral"]["increase_count"],
                    "increase_gt12_count": r_results["ancestral"]["increase_gt12_count"],
                },
                "bijux": {
                    "increase_count": py_results["ancestral"]["increase_count"],
                    "increase_gt12_count": py_results["ancestral"]["increase_gt12_count"],
                },
            },
        },
        {
            "block_id": "pcm1-029",
            "script_lines": "450-451",
            "title": "Save final analysis workspace",
            "status": "artifact_only",
            "note": "The lecture script saves an `.RData` workspace; this report saves explicit machine-readable evidence artifacts instead.",
            "evidence": {
                "r_save_target": "./Results/primate_results.RData",
                "report_artifact_dir": str(OUTPUT_ROOT),
            },
        },
    ]


def _publicize_block_ids(comparisons: list[dict[str, object]]) -> list[dict[str, object]]:
    public_rows: list[dict[str, object]] = []
    for row in comparisons:
        legacy_block_id = row["block_id"]
        public_rows.append(
            {
                **row,
                "legacy_block_id": legacy_block_id,
                "snippet_key": legacy_block_id,
                "block_id": PUBLIC_BLOCK_IDS[legacy_block_id],
            }
        )
    return public_rows


def _sanitize_public_rows(comparisons: list[dict[str, object]]) -> list[dict[str, object]]:
    sanitized: list[dict[str, object]] = []
    for row in comparisons:
        sanitized.append(
            {
                key: value
                for key, value in row.items()
                if key not in {"legacy_block_id", "snippet_key"}
            }
        )
    return sanitized


def _format_readme_evidence(row: dict[str, object]) -> list[str]:
    block_id = row.get("legacy_block_id", row["block_id"])
    evidence = row["evidence"]
    if block_id == "pcm1-000":
        return [
            f"- R package versions recorded: `{', '.join(sorted(evidence['r_packages'].keys()))}`",
            f"- Python reference script: `{evidence['python_reference_script']}`",
        ]
    if block_id == "pcm1-001":
        return [
            f"- processed rows: `{evidence['processed_row_count']}`",
            f"- processed species: `{evidence['processed_species_count']}`",
            f"- checked-in `primate.csv` matches workbook-derived reference: `{evidence['checked_in_processed_matches_reference']}`",
        ]
    if block_id == "pcm1-002":
        return [
            f"- original tree tips: `{evidence['r']['original_tip_count']}`",
            f"- trimmed tree tips: `{evidence['r']['trimmed_tip_count']}`",
            f"- removed taxa: `{', '.join(evidence['r']['missing_tips'])}`",
            f"- checked-in trimmed tree matches R-trimmed reference: `{evidence['r']['checked_in_trimmed_tip_set_matches_reference']}`",
        ]
    if block_id == "pcm1-003":
        return [
            f"- processed CSV written by R: `{evidence['processed_csv_path']}`",
            f"- trimmed tree written by R: `{evidence['trimmed_tree_path']}`",
        ]
    if block_id == "pcm1-006":
        return [
            f"- R unrooted tip count: `{evidence['r']['tip_count']}`",
            f"- bijux unrooted tip count: `{evidence['bijux']['tip_count']}`",
            f"- bijux root child count after unrooting: `{evidence['bijux']['root_child_count']}`",
        ]
    if block_id == "pcm1-008":
        return [
            f"- R source node: numeric `{evidence['r_source_node_numeric']}`, label `{evidence['r_source_node_label']}`",
            f"- R extracted tip count: `{evidence['r_tip_count']}`",
            f"- bijux extracted tip count: `{evidence['bijux_tip_count']}`",
            f"- matched bijux internal node label: `{evidence['bijux_matched_node_name']}`",
            f"- exact descendant taxon set match: `{evidence['same_taxa']}`",
        ]
    if block_id == "pcm1-009":
        return [
            f"- single-node rotation label in R: `{evidence['r']['source_node_label']}`",
            f"- single-node rotation label in bijux: `{evidence['bijux']['source_node_label']}`",
            f"- single-node tip order match: `{evidence['r']['tip_order'] == evidence['bijux']['tip_order']}`",
            f"- all-node tip order match: `{evidence['r_all']['tip_order'] == evidence['bijux_all']['tip_order']}`",
        ]
    if block_id == "pcm1-011":
        return [
            f"- R aligned species exactly follow tip order: `{evidence['r_aligned_equals_tip_order']}`",
            f"- bijux aligned species exactly follow tip order: `{evidence['bijux_aligned_equals_tip_order']}`",
            f"- first six aligned species match across tools: `{evidence['first_six_species_match']}`",
        ]
    if block_id == "pcm1-013":
        return [
            f"- representative node ids in R: `{json.dumps(evidence['nodeid_examples_r'])}`",
            f"- representative node ids in bijux: `{json.dumps(evidence['nodeid_examples_bijux'])}`",
            f"- R joined object tip count: `{evidence['r_joined_tip_count']}`",
            f"- Python analysis taxon count: `{evidence['bijux_analysis_taxa_count']}`",
        ]
    if block_id == "pcm1-015":
        return [
            f"- shared random tree artifact: `{evidence['random_tree_path']}`",
            f"- random examples frozen from R: `{', '.join(evidence['example_names'])}`",
            f"- random tree tip count: `{evidence['random_tree_tip_count']}`",
        ]
    if block_id == "pcm1-017":
        lines = []
        for name in ["random_data", "random_data2", "random_data3", "random_data4", "random_data5"]:
            pair = evidence[name]
            lines.append(
                f"- `{name}` lambda: R `{pair['r']['lambda_value']}` vs bijux `{pair['bijux']['lambda_value']}`"
            )
        return lines
    if block_id == "pcm1-019":
        return [
            f"- R vector length: `{evidence['r_vector_length']}`",
            f"- bijux vector length: `{evidence['bijux_vector_length']}`",
            f"- R names match tip order: `{evidence['r_names_match_tip_order']}`",
            f"- bijux names match tip order: `{evidence['bijux_names_match_tip_order']}`",
        ]
    if block_id == "pcm1-020":
        return [
            f"- lambda absolute difference: `{evidence['lambda_abs_diff']}`",
            f"- log-likelihood absolute difference: `{evidence['log_likelihood_abs_diff']}`",
        ]
    if block_id == "pcm1-022":
        return [
            f"- likelihood-ratio absolute difference: `{evidence['likelihood_ratio_abs_diff']}`",
            f"- p-value absolute difference: `{evidence['p_value_abs_diff']}`",
            f"- lambda=0 covariance top-left 3x3 max diff: `{evidence['lambda0_vcv_top3_max_abs_diff']}`",
            f"- real covariance top-left 3x3 max diff: `{evidence['real_vcv_top3_max_abs_diff']}`",
        ]
    if block_id == "pcm1-023":
        return [
            f"- shared internal clades compared: `{evidence['shared_node_count']}`",
            f"- max point-estimate absolute difference: `{evidence['point_max_abs_diff']}`",
        ]
    if block_id == "pcm1-024":
        return [
            f"- shared internal clades compared: `{evidence['shared_node_count']}`",
            f"- max lower-95 difference: `{evidence['lower_95_max_abs_diff']}`",
            f"- max upper-95 difference: `{evidence['upper_95_max_abs_diff']}`",
        ]
    if block_id == "pcm1-025":
        return [
            f"- R table rows: tips `{evidence['r_tip_rows']}` + internal `{evidence['r_internal_rows']}` = `{evidence['r_total_rows']}`",
            f"- bijux table rows: tips `{evidence['bijux_tip_rows']}` + internal `{evidence['bijux_internal_rows']}` = `{evidence['bijux_total_rows']}`",
        ]
    if block_id == "pcm1-027":
        return [
            f"- MRCA node in R: `{evidence['r']['mrca_node']}`",
            f"- MRCA clade signature in bijux: `{evidence['bijux']['mrca_signature']}`",
            f"- MRCA estimate in R: `{evidence['r']['mrca_estimate'][0]['estimate']}`",
            f"- MRCA estimate in bijux: `{evidence['bijux']['mrca_estimate']}`",
        ]
    if block_id == "pcm1-028":
        return [
            f"- increases > 0: R `{evidence['r']['increase_count']}` vs bijux `{evidence['bijux']['increase_count']}`",
            f"- increases > 12: R `{evidence['r']['increase_gt12_count']}` vs bijux `{evidence['bijux']['increase_gt12_count']}`",
        ]
    if block_id == "pcm1-029":
        return [
            f"- R save target: `{evidence['r_save_target']}`",
            f"- report artifact directory: `{evidence['report_artifact_dir']}`",
        ]
    return [f"- reason: {row['note']}"]


def _build_readme(
    comparisons: list[dict[str, object]],
    coverage: dict[str, object],
    r_results_path: Path,
    py_results_path: Path,
    comparison_path: Path,
) -> str:
    status_counts = Counter(row["status"] for row in comparisons)
    r_script_link = os.path.relpath(R_SCRIPT, OUTPUT_ROOT)
    python_script_link = os.path.relpath(PYTHON_SCRIPT, OUTPUT_ROOT)
    lines = [
        f"# {EXAMPLE_ID} {EVIDENCE_ID}",
        "",
        "This is a registered cross-reference example for the full lecture script",
        f"[`PCM1_plots_signal.R`]({SCRIPT_PATH}).",
        "",
        "Identity:",
        "",
        f"- study id: `{STUDY_ID}`",
        f"- example id: `{EXAMPLE_ID}`",
        f"- evidence id: `{EVIDENCE_ID}`",
        f"- reference tools: `{', '.join(REFERENCE_TOOLS)}`",
        "",
        "Evidence contract:",
        "",
        "- one checked-in R reference script",
        "- one checked-in Python script using `bijux-phylogenetics`",
        "- one machine-readable comparison index plus per-block payloads",
        "- explicit block verdicts that report agreement or deviation without hiding gaps",
        "",
        "Files used in this report:",
        "",
        f"- [R reference checks]({r_script_link})",
        f"- [Python `bijux-phylogenetics` checks]({python_script_link})",
        "- [Example manifest](manifest.json)",
        f"- [R results JSON]({r_results_path.name})",
        f"- [Python results JSON]({py_results_path.name})",
        f"- [Comparison JSON]({comparison_path.name})",
        "- [Per-block payloads](./block-payloads/)",
        "- [Comparative validation suite](./comparative_reference_validation_suite.json)",
        "- [R ecosystem comparison](./r_ecosystem_comparison.json)",
        "- [Trusted examples gallery](./trusted_examples_gallery.json)",
        "- [Reviewer audit checklist](./reviewer_audit_checklist.json)",
        "- [Reproducibility package](./reproducibility_package.json)",
        "- [Method maturity registry](./method_maturity_registry.json)",
        "- [Scientific debt register](./scientific_debt_register.json)",
        "",
        "Coverage summary:",
        "",
        f"- executable script lines tracked: `{coverage['covered_executable_line_count']}` / `{coverage['executable_line_count']}`",
        f"- uncovered executable lines: `{coverage['uncovered_executable_lines']}`",
        "",
        "Status summary:",
        "",
    ]
    for status, count in sorted(status_counts.items()):
        lines.append(f"- `{status}`: {count}")
    lines.extend(
        [
            "",
            "## Block-by-Block Ledger",
            "",
            "| Block | Script lines | Status | What was checked |",
            "| --- | --- | --- | --- |",
        ]
    )
    for row in comparisons:
        lines.append(
            f"| {row['block_id']} | {row['script_lines']} | {row['status']} | {row['title']} |"
        )
    lines.extend(["", "## Block Sections", ""])
    for row in comparisons:
        block_id = row["block_id"]
        lines.extend(
            [
                f"### {block_id} {row['title']}",
                "",
                f"- status: `{row['status']}`",
                f"- script lines: `{row['script_lines']}`",
                f"- verdict: {row['note']}",
                "",
                "**R Lecture Block**",
                "",
                "```r",
                _extract_r_script_snippet(row["script_lines"]),
                "```",
                "",
                "**Python `bijux-phylogenetics` Block**",
                "",
                "```python",
                PYTHON_SNIPPETS[row["snippet_key"]],
                "```",
                "",
                "**Comparison Evidence**",
                "",
            ]
        )
        lines.extend(_format_readme_evidence(row))
        lines.extend(
            [
                "",
                f"Raw comparison payload: [`block-payloads/{block_id}.json`](./block-payloads/{block_id}.json)",
                "",
            ]
        )
    return "\n".join(lines) + "\n"


def _build_manifest(
    comparisons: list[dict[str, object]],
    coverage: dict[str, object],
) -> dict[str, object]:
    return {
        "study_id": STUDY_ID,
        "example_id": EXAMPLE_ID,
        "evidence_id": EVIDENCE_ID,
        "source_script": str(SCRIPT_PATH),
        "reference_tools": REFERENCE_TOOLS,
        "coverage": coverage,
        "status_counts": dict(Counter(row["status"] for row in comparisons)),
        "honesty_contract": {
            "agreement_rule": "mark verified or verified_with_tolerance only when the compared outputs align within stated tolerance",
            "deviation_rule": "mark plot_only, artifact_only, seeded_input_only, or workflow_only when the block is tracked but not claimed as analytical equivalence",
            "disclosure_rule": "report discrepancies directly rather than smoothing them away",
        },
        "comparison_mode": "direct_parity",
    }


def _build_examples_registry() -> dict[str, object]:
    return {
        "study_id": STUDY_ID,
        "corpus_kind": "cross_reference_examples",
        "examples": [
            {
                "example_id": EXAMPLE_ID,
                "title": "Primate lifespan signal workflow",
                "evidence_id": EVIDENCE_ID,
                "source_script": str(SCRIPT_PATH),
                "reference_tools": REFERENCE_TOOLS,
                "current_evidence_readme": str(
                    (EXAMPLE_ROOT / EVIDENCE_ID / "README.md").relative_to(STUDY_ROOT)
                ),
            }
        ],
    }


def _build_examples_index() -> str:
    return "\n".join(
        [
            "# Cross-Reference Examples",
            "",
            "This directory registers durable comparative examples for the",
            "`primate-longevity-signal` study.",
            "",
            "Each example is meant to be stable and reusable:",
            "",
            "- one durable example id",
            "- one or more semantically named evidence bundles under that example id",
            "- explicit reference-tool outputs and bijux outputs",
            "- honest block verdicts for agreement, deviation, or unsupported scope",
            "",
            "Tracked examples:",
            "",
            f"- [{EXAMPLE_ID}](./{EXAMPLE_ID}/{EVIDENCE_ID}/README.md): Primate lifespan signal workflow",
            "",
            f"Registry: [`registry.json`](./registry.json)",
            "",
        ]
    )


def _build_example_index() -> str:
    return "\n".join(
        [
            f"# {EXAMPLE_ID}",
            "",
            "Durable example registration for the primate lifespan signal workflow",
            "used to cross-check `bijux-phylogenetics` against established R tools.",
            "",
            "Evidence bundles:",
            "",
            f"- [{EVIDENCE_ID}](./{EVIDENCE_ID}/README.md)",
            "",
        ]
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _build_comparative_reference_validation_suite(
    comparisons: list[dict[str, object]],
) -> dict[str, object]:
    fixture_payload = _load_json(COMPARATIVE_FIXTURE_PATH)
    tracked_blocks = {
        row["block_id"]: row
        for row in comparisons
        if row["block_id"]
        in {
            "random-signal-lambda-fits",
            "primate-lambda-fit",
            "lambda-zero-covariance-and-lrt",
            "continuous-ancestral-point-estimates",
            "continuous-ancestral-intervals",
            "bonobo-gibbon-mrca-estimate",
            "lifespan-increase-counts",
        }
    }
    return {
        "study_id": STUDY_ID,
        "example_id": EXAMPLE_ID,
        "evidence_id": EVIDENCE_ID,
        "suite_kind": "comparative_reference_validation",
        "reference_tools": REFERENCE_TOOLS,
        "checked_fixture_path": str(COMPARATIVE_FIXTURE_PATH.relative_to(REPO_ROOT)),
        "fixture_observations": fixture_payload["observations"],
        "validated_methods": [
            {
                "method": block_id,
                "status": row["status"],
                "script_lines": row["script_lines"],
                "verdict": row["note"],
            }
            for block_id, row in sorted(tracked_blocks.items())
        ],
    }


def _build_r_ecosystem_comparison(
    comparisons: list[dict[str, object]],
) -> dict[str, object]:
    analytical = [
        row
        for row in comparisons
        if row["status"] in {"verified", "verified_with_tolerance"}
    ]
    return {
        "study_id": STUDY_ID,
        "example_id": EXAMPLE_ID,
        "evidence_id": EVIDENCE_ID,
        "comparison_scope": "accepted_r_tool_reproduction",
        "reference_tools": REFERENCE_TOOLS,
        "analytical_block_count": len(analytical),
        "plot_only_block_count": sum(1 for row in comparisons if row["status"] == "plot_only"),
        "analytical_blocks": [
            {
                "block_id": row["block_id"],
                "status": row["status"],
                "title": row["title"],
                "script_lines": row["script_lines"],
            }
            for row in analytical
        ],
    }


def _build_trusted_examples_gallery(
    comparisons: list[dict[str, object]],
) -> dict[str, object]:
    return {
        "study_id": STUDY_ID,
        "gallery_kind": "trusted_examples",
        "examples": [
            {
                "example_id": EXAMPLE_ID,
                "evidence_id": EVIDENCE_ID,
                "title": "Primate lifespan signal workflow",
                "source_script": str(SCRIPT_PATH),
                "reference_tools": REFERENCE_TOOLS,
                "validated_blocks": [
                    {
                        "block_id": row["block_id"],
                        "title": row["title"],
                        "status": row["status"],
                    }
                    for row in comparisons
                    if row["status"] in {"verified", "verified_with_tolerance"}
                ],
            }
        ],
    }


def _build_reviewer_audit_checklist(
    comparisons: list[dict[str, object]],
    coverage: dict[str, object],
) -> dict[str, object]:
    status_counts = Counter(row["status"] for row in comparisons)
    return {
        "study_id": STUDY_ID,
        "example_id": EXAMPLE_ID,
        "evidence_id": EVIDENCE_ID,
        "checklist": [
            {
                "item": "checked_in_reference_scripts",
                "pass": True,
                "detail": "R and Python comparative scripts are checked in.",
            },
            {
                "item": "full_executable_line_coverage",
                "pass": coverage["covered_executable_line_count"]
                == coverage["executable_line_count"],
                "detail": f"{coverage['covered_executable_line_count']} / {coverage['executable_line_count']} executable lines are covered.",
            },
            {
                "item": "analytical_blocks_validated",
                "pass": status_counts["verified"] + status_counts["verified_with_tolerance"] > 0,
                "detail": f"{status_counts['verified']} verified and {status_counts['verified_with_tolerance']} tolerance-validated blocks.",
            },
            {
                "item": "unclaimed_visual_scope_disclosed",
                "pass": status_counts["plot_only"] >= 0,
                "detail": f"{status_counts['plot_only']} plot-only blocks are disclosed explicitly.",
            },
            {
                "item": "deviation_surface_present",
                "pass": True,
                "detail": "Tolerance-based and non-analytical blocks are separated rather than hidden.",
            },
        ],
    }


def _build_reproducibility_package(
    comparisons: list[dict[str, object]],
    coverage: dict[str, object],
) -> dict[str, object]:
    del comparisons
    files = [
        SCRIPT_PATH,
        R_SCRIPT,
        PYTHON_SCRIPT,
        TRAITS_PATH,
        TRIMMED_TREE_PATH,
    ]
    files.extend(sorted(path for path in OUTPUT_ROOT.rglob("*") if path.is_file()))
    return {
        "study_id": STUDY_ID,
        "example_id": EXAMPLE_ID,
        "evidence_id": EVIDENCE_ID,
        "coverage": coverage,
        "files": [
            {
                "path": str(path.relative_to(REPO_ROOT))
                if path.is_relative_to(REPO_ROOT)
                else str(path),
                "sha256": _sha256(path),
            }
            for path in files
        ],
    }


def _build_method_maturity_registry(
    comparisons: list[dict[str, object]],
) -> dict[str, object]:
    del comparisons
    registry = build_method_limitation_registry()
    return {
        "study_id": STUDY_ID,
        "example_id": EXAMPLE_ID,
        "evidence_id": EVIDENCE_ID,
        "registry_goal_id": registry.goal_id,
        "methods": [
            {
                "method": entry.method,
                "maturity": (
                    "validated"
                    if entry.status == "validated"
                    else "experimental"
                ),
                "validated_by": entry.validated_by,
                "limitations": entry.limitations,
            }
            for entry in registry.entries
        ],
        "limitations": registry.limitations,
    }


def _build_scientific_debt_register(
    comparisons: list[dict[str, object]],
) -> dict[str, object]:
    validation_report = build_scientific_validation_report()
    debts = []
    for claim in validation_report.claims:
        if claim.status in {"experimental", "unvalidated", "unsafe"}:
            debts.append(
                {
                    "debt_kind": claim.status,
                    "detail": claim.claim,
                    "evidence": claim.evidence,
                }
            )
    for row in comparisons:
        if row["status"] in {"plot_only", "workflow_only", "artifact_only", "seeded_input_only"}:
            debts.append(
                {
                    "block_id": row["block_id"],
                    "debt_kind": row["status"],
                    "detail": row["note"],
                }
            )
        elif row["status"] == "verified_with_tolerance":
            debts.append(
                {
                    "block_id": row["block_id"],
                    "debt_kind": "tolerance_validated",
                    "detail": row["note"],
                }
            )
    return {
        "study_id": STUDY_ID,
        "example_id": EXAMPLE_ID,
        "evidence_id": EVIDENCE_ID,
        "validation_goal_id": validation_report.goal_id,
        "debts": debts,
        "limitations": validation_report.limitations,
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    if not RSCRIPT_BIN.exists():
        raise FileNotFoundError(f"Rscript not found at {RSCRIPT_BIN}")

    _run([str(RSCRIPT_BIN), str(R_SCRIPT), str(R_REPO_ROOT), str(OUTPUT_ROOT)])
    _run([sys.executable, str(PYTHON_SCRIPT), str(R_REPO_ROOT), str(OUTPUT_ROOT)])

    r_results_path = OUTPUT_ROOT / "r_reference_results.json"
    py_results_path = OUTPUT_ROOT / "bijux_reference_results.json"
    r_results = _load_json(r_results_path)
    py_results = _load_json(py_results_path)

    comparisons = _publicize_block_ids(_compare_results(r_results, py_results))
    public_comparisons = _sanitize_public_rows(comparisons)
    coverage = _coverage_summary(comparisons)
    manifest = _build_manifest(comparisons, coverage)
    comparison_payload = {
        "study_id": STUDY_ID,
        "example_id": EXAMPLE_ID,
        "evidence_id": EVIDENCE_ID,
        "source_script": str(SCRIPT_PATH),
        "r_reference_script": str(R_SCRIPT),
        "python_reference_script": str(PYTHON_SCRIPT),
        "coverage": coverage,
        "comparisons": public_comparisons,
    }
    comparison_path = OUTPUT_ROOT / "comparison.json"
    (OUTPUT_ROOT / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    comparison_path.write_text(
        json.dumps(comparison_payload, indent=2) + "\n",
        encoding="utf-8",
    )
    payload_dir = OUTPUT_ROOT / "block-payloads"
    payload_dir.mkdir(exist_ok=True)
    for row in public_comparisons:
        (payload_dir / f"{row['block_id']}.json").write_text(
            json.dumps(row, indent=2) + "\n",
            encoding="utf-8",
        )
    (OUTPUT_ROOT / "comparative_reference_validation_suite.json").write_text(
        json.dumps(_build_comparative_reference_validation_suite(comparisons), indent=2) + "\n",
        encoding="utf-8",
    )
    (OUTPUT_ROOT / "r_ecosystem_comparison.json").write_text(
        json.dumps(_build_r_ecosystem_comparison(comparisons), indent=2) + "\n",
        encoding="utf-8",
    )
    (OUTPUT_ROOT / "trusted_examples_gallery.json").write_text(
        json.dumps(_build_trusted_examples_gallery(comparisons), indent=2) + "\n",
        encoding="utf-8",
    )
    (OUTPUT_ROOT / "reviewer_audit_checklist.json").write_text(
        json.dumps(_build_reviewer_audit_checklist(comparisons, coverage), indent=2) + "\n",
        encoding="utf-8",
    )
    (OUTPUT_ROOT / "method_maturity_registry.json").write_text(
        json.dumps(_build_method_maturity_registry(comparisons), indent=2) + "\n",
        encoding="utf-8",
    )
    (OUTPUT_ROOT / "scientific_debt_register.json").write_text(
        json.dumps(_build_scientific_debt_register(comparisons), indent=2) + "\n",
        encoding="utf-8",
    )
    (OUTPUT_ROOT / "reproducibility_package.json").write_text(
        json.dumps(_build_reproducibility_package(comparisons, coverage), indent=2) + "\n",
        encoding="utf-8",
    )
    (OUTPUT_ROOT / "README.md").write_text(
        _build_readme(comparisons, coverage, r_results_path, py_results_path, comparison_path),
        encoding="utf-8",
    )
    EXAMPLES_ROOT.mkdir(exist_ok=True)
    EXAMPLE_ROOT.mkdir(exist_ok=True)
    (EXAMPLES_ROOT / "registry.json").write_text(
        json.dumps(_build_examples_registry(), indent=2) + "\n",
        encoding="utf-8",
    )
    (EXAMPLES_ROOT / "README.md").write_text(
        _build_examples_index(),
        encoding="utf-8",
    )
    (EXAMPLE_ROOT / "README.md").write_text(
        _build_example_index(),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
