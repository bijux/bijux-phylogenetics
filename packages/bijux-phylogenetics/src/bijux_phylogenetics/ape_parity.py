from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, replace
import json
from importlib import metadata
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from bijux_phylogenetics.clades import extract_tree_clades, extract_tree_set_clades
from bijux_phylogenetics.compare.structural_parity import (
    compare_tree_sets_structurally,
    compare_tree_structurally,
)
from bijux_phylogenetics.distance import compute_pairwise_genetic_distance_matrix
from bijux_phylogenetics.diagnostics.validation import inspect_tree_path
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.core.topology import root_tree_on_outgroup
from bijux_phylogenetics.io.fasta import load_fasta_alignment, translate_coding_alignment
from bijux_phylogenetics.io.newick import (
    dumps_newick,
    load_newick_tree_set,
    write_newick,
    write_newick_tree_set,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.shared_dna_alignment_fixtures import (
    get_shared_dna_alignment_fixture,
)
from bijux_phylogenetics.shared_tree_fixtures import get_shared_tree_fixture
from bijux_phylogenetics.shared_tree_set_fixtures import get_shared_tree_set_fixture


@dataclass(frozen=True, slots=True)
class ApeParityCase:
    """One governed live `ape` parity case."""

    case_id: str
    fixture_kind: str
    fixture_id: str
    function_name: str
    python_function_name: str
    operation: str
    input_fixture: Path
    tolerance: float
    expected_status: str = "ok"
    pairwise_deletion: bool | None = None
    genetic_code_id: int | None = None
    outgroup_taxa: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class ApeParityObservation:
    """One live parity comparison between Bijux and `ape`."""

    case_id: str
    fixture_kind: str
    fixture_id: str
    function_name: str
    python_function_name: str
    input_fixture: Path
    tolerance: float
    r_version: str | None
    ape_version: str | None
    bijux_version: str
    bijux_commit: str | None
    status: str
    passed: bool
    mismatch_reason: str | None
    reproducible_artifact_root: Path | None
    reference_summary: dict[str, object] | None
    bijux_summary: dict[str, object] | None
    reference_error: dict[str, object] | None
    bijux_error: dict[str, object] | None


@dataclass(frozen=True, slots=True)
class ApeParitySummaryRow:
    """One function-level summary across governed `ape` parity cases."""

    function_name: str
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int


@dataclass(slots=True)
class ApeParityReport:
    """Aggregate report for governed live `ape` parity cases."""

    observations: list[ApeParityObservation]
    summary_rows: list[ApeParitySummaryRow]
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int
    all_passed: bool
    limitations: list[str]


def _package_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[4]


def _fixtures_root() -> Path:
    return _package_root() / "tests" / "fixtures"


def _ape_runner_path() -> Path:
    return (
        Path(__file__).resolve().parent
        / "resources"
        / "reference"
        / "ape_parity_runner.R"
    )


def _failure_root() -> Path:
    return _repository_root() / "artifacts" / "ape-parity-failures"


def _reference_environment() -> dict[str, str]:
    environment = dict(os.environ)
    r_library = _repository_root() / "artifacts" / "r-lib"
    if "R_LIBS_USER" not in environment and r_library.is_dir():
        environment["R_LIBS_USER"] = str(r_library)
    return environment


def _bijux_version() -> str:
    try:
        return metadata.version("bijux-phylogenetics")
    except metadata.PackageNotFoundError:
        return "0.1.0"


def _bijux_commit() -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"],
        capture_output=True,
        check=False,
        cwd=_repository_root(),
        text=True,
    )
    if result.returncode != 0:
        return None
    commit = result.stdout.strip()
    return commit or None


def list_ape_parity_cases(fixtures_root: Path | None = None) -> list[ApeParityCase]:
    """Return the governed live `ape` parity cases."""
    root = _fixtures_root() if fixtures_root is None else fixtures_root

    def fixture_path(fixture_kind: str, fixture_id: str) -> Path:
        if fixture_kind == "tree":
            fixture = get_shared_tree_fixture(fixture_id)
        elif fixture_kind == "tree-set":
            fixture = get_shared_tree_set_fixture(fixture_id)
        elif fixture_kind == "dna-alignment":
            fixture = get_shared_dna_alignment_fixture(fixture_id)
        else:
            raise ValueError(f"unsupported ape parity fixture kind '{fixture_kind}'")
        if fixtures_root is None:
            return fixture.path
        return root / fixture.relative_path

    return [
        ApeParityCase(
            case_id="read-tree-balanced-rooted-ultrametric",
            fixture_kind="tree",
            fixture_id="balanced_rooted_ultrametric",
            function_name="ape::read.tree",
            python_function_name="load_tree+extract_tree_clades",
            operation="read-tree-structure",
            input_fixture=fixture_path("tree", "balanced_rooted_ultrametric"),
            tolerance=0.0,
        ),
        ApeParityCase(
            case_id="read-tree-unrooted-branch-length",
            fixture_kind="tree",
            fixture_id="unrooted_branch_length_tree",
            function_name="ape::read.tree",
            python_function_name="load_tree+extract_tree_clades",
            operation="read-tree-structure",
            input_fixture=fixture_path("tree", "unrooted_branch_length_tree"),
            tolerance=0.0,
        ),
        ApeParityCase(
            case_id="read-tree-internal-node-labels",
            fixture_kind="tree",
            fixture_id="internal_node_labels",
            function_name="ape::read.tree",
            python_function_name="load_tree+extract_tree_clades",
            operation="read-tree-structure",
            input_fixture=fixture_path("tree", "internal_node_labels"),
            tolerance=0.0,
        ),
        ApeParityCase(
            case_id="read-tree-support-labels",
            fixture_kind="tree",
            fixture_id="branch_support_labels",
            function_name="ape::read.tree",
            python_function_name="load_tree+extract_tree_clades",
            operation="read-tree-structure",
            input_fixture=fixture_path("tree", "branch_support_labels"),
            tolerance=1e-12,
        ),
        ApeParityCase(
            case_id="read-tree-quoted-taxon-labels",
            fixture_kind="tree",
            fixture_id="quoted_taxon_labels",
            function_name="ape::read.tree",
            python_function_name="load_tree+extract_tree_clades",
            operation="read-tree-structure",
            input_fixture=fixture_path("tree", "quoted_taxon_labels"),
            tolerance=0.0,
        ),
        ApeParityCase(
            case_id="read-tree-multiple-trees",
            fixture_kind="tree-set",
            fixture_id="basic_newick_tree_set",
            function_name="ape::read.tree",
            python_function_name="extract_tree_set_clades",
            operation="read-tree-set-structure",
            input_fixture=fixture_path("tree-set", "basic_newick_tree_set"),
            tolerance=0.0,
        ),
        ApeParityCase(
            case_id="read-tree-malformed-newick",
            fixture_kind="tree",
            fixture_id="malformed_unbalanced_parentheses",
            function_name="ape::read.tree",
            python_function_name="load_tree",
            operation="read-tree-structure",
            input_fixture=fixture_path("tree", "malformed_unbalanced_parentheses"),
            tolerance=0.0,
            expected_status="parse-error",
        ),
        ApeParityCase(
            case_id="write-tree-balanced-rooted-ultrametric",
            fixture_kind="tree",
            fixture_id="balanced_rooted_ultrametric",
            function_name="ape::write.tree",
            python_function_name="write_newick+ape::read.tree",
            operation="write-tree-structure",
            input_fixture=fixture_path("tree", "balanced_rooted_ultrametric"),
            tolerance=0.0,
        ),
        ApeParityCase(
            case_id="write-tree-unrooted-branch-length",
            fixture_kind="tree",
            fixture_id="unrooted_branch_length_tree",
            function_name="ape::write.tree",
            python_function_name="write_newick+ape::read.tree",
            operation="write-tree-structure",
            input_fixture=fixture_path("tree", "unrooted_branch_length_tree"),
            tolerance=0.0,
        ),
        ApeParityCase(
            case_id="write-tree-internal-node-labels",
            fixture_kind="tree",
            fixture_id="internal_node_labels",
            function_name="ape::write.tree",
            python_function_name="write_newick+ape::read.tree",
            operation="write-tree-structure",
            input_fixture=fixture_path("tree", "internal_node_labels"),
            tolerance=0.0,
        ),
        ApeParityCase(
            case_id="write-tree-support-labels",
            fixture_kind="tree",
            fixture_id="branch_support_labels",
            function_name="ape::write.tree",
            python_function_name="write_newick+ape::read.tree",
            operation="write-tree-structure",
            input_fixture=fixture_path("tree", "branch_support_labels"),
            tolerance=1e-12,
        ),
        ApeParityCase(
            case_id="write-tree-quoted-taxon-labels",
            fixture_kind="tree",
            fixture_id="quoted_taxon_labels",
            function_name="ape::write.tree",
            python_function_name="write_newick+ape::read.tree",
            operation="write-tree-structure",
            input_fixture=fixture_path("tree", "quoted_taxon_labels"),
            tolerance=0.0,
        ),
        ApeParityCase(
            case_id="write-tree-multiple-trees",
            fixture_kind="tree-set",
            fixture_id="basic_newick_tree_set",
            function_name="ape::write.tree",
            python_function_name="write_newick_tree_set+ape::read.tree",
            operation="write-tree-set-structure",
            input_fixture=fixture_path("tree-set", "basic_newick_tree_set"),
            tolerance=0.0,
        ),
        ApeParityCase(
            case_id="root-tree-single-outgroup-tip",
            fixture_kind="tree",
            fixture_id="outgroup_rootable_unrooted",
            function_name="ape::root",
            python_function_name="root_tree_on_outgroup",
            operation="root-tree-outgroup",
            input_fixture=fixture_path("tree", "outgroup_rootable_unrooted"),
            tolerance=1e-12,
            outgroup_taxa=("D",),
        ),
        ApeParityCase(
            case_id="root-tree-multiple-outgroup-tips",
            fixture_kind="tree",
            fixture_id="outgroup_rootable_unrooted",
            function_name="ape::root",
            python_function_name="root_tree_on_outgroup",
            operation="root-tree-outgroup",
            input_fixture=fixture_path("tree", "outgroup_rootable_unrooted"),
            tolerance=1e-12,
            outgroup_taxa=("C", "D"),
        ),
        ApeParityCase(
            case_id="root-tree-already-rooted",
            fixture_kind="tree",
            fixture_id="outgroup_rooted_on_d",
            function_name="ape::root",
            python_function_name="root_tree_on_outgroup",
            operation="root-tree-outgroup",
            input_fixture=fixture_path("tree", "outgroup_rooted_on_d"),
            tolerance=1e-12,
            outgroup_taxa=("D",),
        ),
        ApeParityCase(
            case_id="root-tree-missing-outgroup",
            fixture_kind="tree",
            fixture_id="outgroup_rootable_unrooted",
            function_name="ape::root",
            python_function_name="root_tree_on_outgroup",
            operation="root-tree-outgroup",
            input_fixture=fixture_path("tree", "outgroup_rootable_unrooted"),
            tolerance=0.0,
            expected_status="rooting-error",
            outgroup_taxa=("Z",),
        ),
        ApeParityCase(
            case_id="root-tree-non-monophyletic-outgroup",
            fixture_kind="tree",
            fixture_id="outgroup_rootable_unrooted",
            function_name="ape::root",
            python_function_name="root_tree_on_outgroup",
            operation="root-tree-outgroup",
            input_fixture=fixture_path("tree", "outgroup_rootable_unrooted"),
            tolerance=0.0,
            expected_status="rooting-error",
            outgroup_taxa=("B", "D"),
        ),
        ApeParityCase(
            case_id="dna-base-frequency-lowercase",
            fixture_kind="dna-alignment",
            fixture_id="lowercase_aligned_dna",
            function_name="ape::base.freq",
            python_function_name="load_fasta_alignment+ape-style-base-frequency",
            operation="dna-base-frequency",
            input_fixture=fixture_path("dna-alignment", "lowercase_aligned_dna"),
            tolerance=1e-12,
        ),
        ApeParityCase(
            case_id="dna-base-frequency-ambiguity",
            fixture_kind="dna-alignment",
            fixture_id="dna_with_ambiguity",
            function_name="ape::base.freq",
            python_function_name="load_fasta_alignment+ape-style-base-frequency",
            operation="dna-base-frequency",
            input_fixture=fixture_path("dna-alignment", "dna_with_ambiguity"),
            tolerance=1e-12,
        ),
        ApeParityCase(
            case_id="dna-raw-distance-clean",
            fixture_kind="dna-alignment",
            fixture_id="clean_aligned_dna",
            function_name="ape::dist.dna",
            python_function_name="compute_pairwise_genetic_distance_matrix",
            operation="dna-raw-distance",
            input_fixture=fixture_path("dna-alignment", "clean_aligned_dna"),
            tolerance=1e-12,
            pairwise_deletion=False,
        ),
        ApeParityCase(
            case_id="dna-raw-distance-gaps",
            fixture_kind="dna-alignment",
            fixture_id="dna_with_gaps",
            function_name="ape::dist.dna",
            python_function_name="compute_pairwise_genetic_distance_matrix",
            operation="dna-raw-distance",
            input_fixture=fixture_path("dna-alignment", "dna_with_gaps"),
            tolerance=1e-12,
            pairwise_deletion=True,
        ),
        ApeParityCase(
            case_id="dna-raw-distance-identical",
            fixture_kind="dna-alignment",
            fixture_id="identical_sequences",
            function_name="ape::dist.dna",
            python_function_name="compute_pairwise_genetic_distance_matrix",
            operation="dna-raw-distance",
            input_fixture=fixture_path("dna-alignment", "identical_sequences"),
            tolerance=1e-12,
            pairwise_deletion=False,
        ),
        ApeParityCase(
            case_id="dna-raw-distance-high-divergence",
            fixture_kind="dna-alignment",
            fixture_id="high_divergence_sequences",
            function_name="ape::dist.dna",
            python_function_name="compute_pairwise_genetic_distance_matrix",
            operation="dna-raw-distance",
            input_fixture=fixture_path("dna-alignment", "high_divergence_sequences"),
            tolerance=1e-12,
            pairwise_deletion=False,
        ),
        ApeParityCase(
            case_id="dna-raw-distance-missing-data",
            fixture_kind="dna-alignment",
            fixture_id="dna_with_missing_data",
            function_name="ape::dist.dna",
            python_function_name="compute_pairwise_genetic_distance_matrix",
            operation="dna-raw-distance",
            input_fixture=fixture_path("dna-alignment", "dna_with_missing_data"),
            tolerance=1e-12,
            pairwise_deletion=True,
        ),
        ApeParityCase(
            case_id="dna-translation-valid-frame",
            fixture_kind="dna-alignment",
            fixture_id="coding_valid_reading_frame",
            function_name="ape::trans",
            python_function_name="translate_coding_alignment",
            operation="dna-translation",
            input_fixture=fixture_path("dna-alignment", "coding_valid_reading_frame"),
            tolerance=0.0,
            genetic_code_id=1,
        ),
        ApeParityCase(
            case_id="dna-translation-internal-stop",
            fixture_kind="dna-alignment",
            fixture_id="coding_internal_stop",
            function_name="ape::trans",
            python_function_name="translate_coding_alignment",
            operation="dna-translation",
            input_fixture=fixture_path("dna-alignment", "coding_internal_stop"),
            tolerance=0.0,
            genetic_code_id=1,
        ),
        ApeParityCase(
            case_id="dna-translation-terminal-stop",
            fixture_kind="dna-alignment",
            fixture_id="coding_terminal_stop",
            function_name="ape::trans",
            python_function_name="translate_coding_alignment",
            operation="dna-translation",
            input_fixture=fixture_path("dna-alignment", "coding_terminal_stop"),
            tolerance=0.0,
            genetic_code_id=1,
        ),
    ]


def _build_case_lookup(fixtures_root: Path | None = None) -> dict[str, ApeParityCase]:
    return {case.case_id: case for case in list_ape_parity_cases(fixtures_root)}


def _selected_cases(
    *,
    case_ids: list[str] | None,
    fixtures_root: Path | None = None,
) -> list[ApeParityCase]:
    cases = _build_case_lookup(fixtures_root)
    if case_ids is None:
        return list(cases.values())
    selected: list[ApeParityCase] = []
    for case_id in case_ids:
        try:
            selected.append(cases[case_id])
        except KeyError as error:
            supported = ", ".join(sorted(cases))
            raise ValueError(
                f"unsupported ape parity case '{case_id}'; expected one of: {supported}"
            ) from error
    return selected


def _write_case_file(path: Path, case: ApeParityCase) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "case_id": case.case_id,
                "fixture_kind": case.fixture_kind,
                "fixture_id": case.fixture_id,
                "function_name": case.function_name,
                "operation": case.operation,
                "input_fixture": str(case.input_fixture),
                "tolerance": case.tolerance,
                "expected_status": case.expected_status,
                "pairwise_deletion": case.pairwise_deletion,
                "genetic_code_id": case.genetic_code_id,
                "outgroup_taxa": list(case.outgroup_taxa),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def _node_kind_order(node_kind: str) -> int:
    return {"root": 0, "internal": 1, "tip": 2}.get(node_kind, 9)


def _clade_rows_to_parity_rows(rows) -> list[dict[str, object]]:
    parity_rows = [
        {
            "tree_index": "" if row.tree_index is None else row.tree_index,
            "node_kind": row.node_kind,
            "clade_id": row.clade_id,
            "node_label": "" if row.node_label is None else row.node_label,
            "taxon_count": row.taxon_count,
            "taxa": "|".join(row.taxa),
            "support": "" if row.support is None else row.support,
            "branch_length": "" if row.branch_length is None else row.branch_length,
        }
        for row in rows
    ]
    return _sort_parity_rows(parity_rows)


def _build_bijux_tree_structure(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    tree = load_tree(input_fixture)
    inspection = inspect_tree_path(input_fixture)
    clades = extract_tree_clades(input_fixture)
    return _tree_structure_payload(tree, inspection.rooted, clades.rows)


def _tree_structure_payload(
    tree: TreeNode | PhyloTree,
    rooted: bool | None,
    clade_rows,
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    phylo_tree = tree if isinstance(tree, PhyloTree) else PhyloTree(root=tree, rooted=rooted)
    summary = {
        "tree_count": 1,
        "tip_count": phylo_tree.tip_count,
        "internal_node_count": phylo_tree.internal_node_count,
        "edge_count": phylo_tree.tip_count + phylo_tree.internal_node_count - 1,
        "rooted": rooted,
        "tip_labels": phylo_tree.tip_names,
        "branch_length_count": sum(
            1
            for branch_length in phylo_tree.branch_lengths()
            if branch_length is not None
        ),
    }
    return summary, _clade_rows_to_parity_rows(clade_rows), dumps_newick(phylo_tree)


def _build_bijux_tree_set_structure(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]], None]:
    clades = extract_tree_set_clades(input_fixture)
    parity_rows = _clade_rows_to_parity_rows(clades.rows)
    tree_indices = sorted(
        {
            row["tree_index"]
            for row in parity_rows
            if row["tree_index"] != ""
        }
    )
    first_tree_tip_labels = [
        row["node_label"]
        for row in parity_rows
        if row["tree_index"] == 1 and row["node_kind"] == "tip"
    ]
    summary = {
        "tree_count": clades.tree_count,
        "source_format": clades.source_format,
        "tree_indices": tree_indices,
        "shared_tip_labels": sorted(first_tree_tip_labels),
        "unique_tip_label_count": len(first_tree_tip_labels),
    }
    return summary, parity_rows, None


def _build_bijux_root_outgroup_structure(
    input_fixture: Path,
    *,
    outgroup_taxa: tuple[str, ...],
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    rooted_tree, _report = root_tree_on_outgroup(
        input_fixture,
        outgroup_taxa=list(outgroup_taxa),
    )
    with tempfile.TemporaryDirectory(prefix="bijux-ape-root-") as tmpdir:
        rooted_path = Path(tmpdir) / "rooted.nwk"
        write_newick(rooted_path, rooted_tree)
        clades = extract_tree_clades(rooted_path)
    return _tree_structure_payload(rooted_tree, rooted_tree.rooted, clades.rows)


def _materialize_reference_input(case: ApeParityCase, working_root: Path) -> Path:
    reference_input_path = working_root / "bijux-reference-input.nwk"
    if case.operation == "write-tree-structure":
        tree = load_tree(case.input_fixture)
        write_newick(reference_input_path, tree)
        return reference_input_path
    if case.operation == "write-tree-set-structure":
        trees = load_newick_tree_set(case.input_fixture)
        write_newick_tree_set(reference_input_path, trees)
        return reference_input_path
    return case.input_fixture


def _ape_base_frequency_rows(input_fixture: Path) -> list[dict[str, object]]:
    state_order = [
        "a",
        "c",
        "g",
        "t",
        "r",
        "m",
        "w",
        "s",
        "k",
        "y",
        "v",
        "h",
        "d",
        "b",
        "n",
        "-",
        "?",
    ]
    counts = {state: 0 for state in state_order}
    records = load_fasta_alignment(input_fixture)
    total = 0
    for record in records:
        for residue in record.sequence:
            normalized = residue.lower().replace("u", "t")
            if normalized not in counts:
                continue
            counts[normalized] += 1
            total += 1
    if total == 0:
        return [{"state": state, "frequency": 0.0} for state in state_order]
    return [
        {"state": state, "frequency": counts[state] / total}
        for state in state_order
    ]


def _build_bijux_base_frequency_summary(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    records = load_fasta_alignment(input_fixture)
    rows = _ape_base_frequency_rows(input_fixture)
    return {
        "sequence_count": len(records),
        "alignment_length": len(records[0].sequence),
        "state_count": len(rows),
    }, rows


def _build_bijux_distance_rows(
    input_fixture: Path,
    *,
    pairwise_deletion: bool,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_pairwise_genetic_distance_matrix(
        input_fixture,
        model="p-distance",
        gap_handling="pairwise-deletion" if pairwise_deletion else "complete-deletion",
        ambiguity_policy="ignore",
    )
    pair_lookup = {
        (row.left_identifier, row.right_identifier): row.distance for row in report.pairs
    }
    rows: list[dict[str, object]] = []
    for left_identifier in report.identifiers:
        for right_identifier in report.identifiers:
            distance = pair_lookup.get((left_identifier, right_identifier))
            if distance is None:
                distance = pair_lookup.get((right_identifier, left_identifier))
            rows.append(
                {
                    "left_identifier": left_identifier,
                    "right_identifier": right_identifier,
                    "distance": distance,
                }
            )
    return {
        "sequence_count": len(report.identifiers),
        "alignment_length": report.alignment_length,
        "pairwise_deletion": pairwise_deletion,
    }, rows


def _build_bijux_translation_rows(
    input_fixture: Path,
    *,
    genetic_code_id: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    translated, report = translate_coding_alignment(
        input_fixture, genetic_code=genetic_code_id
    )
    return {
        "sequence_count": report.translated_sequence_count,
        "translated_length": report.translated_alignment_length,
        "stop_codon_count": report.stop_codon_count,
    }, [
        {
            "identifier": row.identifier,
            "amino_acid_sequence": row.sequence,
        }
        for row in translated
    ]


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _normalize_newick_label(label: str) -> str:
    if len(label) >= 2 and label.startswith("'") and label.endswith("'"):
        return label[1:-1].replace("''", "'")
    return label


def _normalize_expected_label(
    label: str,
    *,
    expected_tip_labels: set[str] | None,
) -> str:
    normalized = _normalize_newick_label(label)
    if (
        expected_tip_labels
        and normalized not in expected_tip_labels
        and normalized.replace("_", " ") in expected_tip_labels
    ):
        return normalized.replace("_", " ")
    return normalized


def _normalize_joined_labels(
    value: str,
    *,
    expected_tip_labels: set[str] | None,
) -> str:
    if value == "":
        return value
    labels = [
        _normalize_expected_label(label, expected_tip_labels=expected_tip_labels)
        for label in value.split("|")
    ]
    return "|".join(sorted(labels))


def _normalize_reference_summary(summary: dict[str, object]) -> dict[str, object]:
    normalized = dict(summary)
    tip_labels = normalized.get("tip_labels")
    if isinstance(tip_labels, list):
        expected_tip_labels = {
            _normalize_newick_label(str(label)) for label in tip_labels
        }
        normalized["tip_labels"] = [
            _normalize_expected_label(
                str(label),
                expected_tip_labels=expected_tip_labels,
            )
            for label in tip_labels
        ]
    return normalized


def _coerce_table_cell(value: str) -> object:
    if value == "":
        return ""
    if value in {"true", "false"}:
        return value == "true"
    if re.fullmatch(r"-?\d+", value):
        return int(value)
    if re.fullmatch(r"-?(?:\d+\.\d*|\d*\.\d+)(?:[eE][+-]?\d+)?", value):
        return float(value)
    return value


def _sort_parity_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(
        rows,
        key=lambda row: (
            0 if row["tree_index"] == "" else int(row["tree_index"]),
            _node_kind_order(str(row["node_kind"])),
            str(row["clade_id"]),
            str(row["node_label"]),
        ),
    )


def _load_rows_table(
    path: Path,
    *,
    expected_tip_labels: set[str] | None = None,
    sort_rows: bool = False,
) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    normalized_rows = [
        {
            key: _normalize_expected_label(value, expected_tip_labels=expected_tip_labels)
            if key.endswith("label")
            else _normalize_joined_labels(
                value,
                expected_tip_labels=expected_tip_labels,
            )
            if key in {"clade_id", "taxa"}
            else _coerce_table_cell(value)
            for key, value in row.items()
        }
        for row in rows
    ]
    if sort_rows:
        return _sort_parity_rows(normalized_rows)
    return normalized_rows


def _compare_scalar(expected: object, observed: object, *, tolerance: float) -> bool:
    if isinstance(expected, (int, float)) and isinstance(observed, (int, float)):
        return abs(float(expected) - float(observed)) <= tolerance
    return expected == observed


def _compare_json(expected: object, observed: object, *, tolerance: float) -> bool:
    if isinstance(expected, dict) and isinstance(observed, dict):
        if set(expected) != set(observed):
            return False
        return all(
            _compare_json(expected[key], observed[key], tolerance=tolerance)
            for key in expected
        )
    if isinstance(expected, list) and isinstance(observed, list):
        if len(expected) != len(observed):
            return False
        return all(
            _compare_json(left, right, tolerance=tolerance)
            for left, right in zip(expected, observed, strict=True)
        )
    return _compare_scalar(expected, observed, tolerance=tolerance)


def _normalize_tree_labels(
    node: TreeNode,
    *,
    expected_tip_labels: set[str] | None,
) -> None:
    if node.name is not None:
        normalized = _normalize_newick_label(node.name)
        if (
            expected_tip_labels
            and normalized not in expected_tip_labels
            and normalized.replace("_", " ") in expected_tip_labels
        ):
            normalized = normalized.replace("_", " ")
        node.name = normalized
    for child in node.children:
        _normalize_tree_labels(child, expected_tip_labels=expected_tip_labels)


def _canonical_newick(
    path: Path,
    *,
    expected_tip_labels: set[str] | None = None,
) -> str:
    tree = load_tree(path)
    _normalize_tree_labels(tree.root, expected_tip_labels=expected_tip_labels)
    return dumps_newick(tree)


def _copy_if_exists(source: Path, destination: Path) -> None:
    if source.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def _write_rows_table(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def _persist_failure_bundle(
    *,
    failure_root: Path,
    case: ApeParityCase,
    case_file: Path,
    execution_root: Path,
    execution_payload: dict[str, object] | None,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
    reference_error: dict[str, object] | None,
    bijux_error: dict[str, object] | None,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
    bijux_normalized_text: str | None,
    mismatch_reason: str,
) -> Path:
    artifact_root = failure_root / case.case_id
    artifact_root.mkdir(parents=True, exist_ok=True)
    _copy_if_exists(case_file, artifact_root / "case.json")
    _copy_if_exists(
        execution_root.parent / "bijux-reference-input.nwk",
        artifact_root / "bijux-reference-input.nwk",
    )
    _copy_if_exists(execution_root / "reference-execution.json", artifact_root / "reference-execution.json")
    if execution_payload is not None:
        outputs = execution_payload.get("outputs")
        if isinstance(outputs, dict):
            for path_string in outputs.values():
                if isinstance(path_string, str):
                    source = Path(path_string)
                    _copy_if_exists(source, artifact_root / f"reference-{source.name}")
    if execution_payload is not None:
        _write_json(artifact_root / "reference-execution.observed.json", execution_payload)
    if reference_summary is not None:
        _write_json(artifact_root / "reference-summary.observed.json", reference_summary)
    if reference_rows is not None:
        _write_json(artifact_root / "reference-rows.observed.json", reference_rows)
    if bijux_summary is not None:
        _write_json(artifact_root / "bijux-summary.json", bijux_summary)
    if reference_error is not None:
        _write_json(artifact_root / "reference-error.observed.json", reference_error)
    if bijux_error is not None:
        _write_json(artifact_root / "bijux-error.json", bijux_error)
    if bijux_rows is not None:
        _write_json(artifact_root / "bijux-rows.json", bijux_rows)
    if bijux_normalized_text is not None:
        (artifact_root / "bijux-normalized.txt").write_text(
            f"{bijux_normalized_text}\n",
            encoding="utf-8",
        )
    _write_json(
        artifact_root / "comparison.json",
        {
            "case_id": case.case_id,
            "function_name": case.function_name,
            "mismatch_reason": mismatch_reason,
        },
    )
    return artifact_root


def _build_bijux_case_payload(
    case: ApeParityCase,
) -> tuple[dict[str, object], list[dict[str, object]] | None, str | None]:
    if case.operation in {"read-tree-structure", "write-tree-structure"}:
        summary, rows, normalized_text = _build_bijux_tree_structure(case.input_fixture)
        return summary, rows, normalized_text
    if case.operation == "root-tree-outgroup":
        summary, rows, normalized_text = _build_bijux_root_outgroup_structure(
            case.input_fixture,
            outgroup_taxa=case.outgroup_taxa,
        )
        return summary, rows, normalized_text
    if case.operation in {"read-tree-set-structure", "write-tree-set-structure"}:
        summary, rows, normalized_text = _build_bijux_tree_set_structure(case.input_fixture)
        return summary, rows, normalized_text
    if case.operation == "dna-base-frequency":
        summary, rows = _build_bijux_base_frequency_summary(case.input_fixture)
        return summary, rows, None
    if case.operation == "dna-raw-distance":
        if case.pairwise_deletion is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing pairwise deletion policy"
            )
        summary, rows = _build_bijux_distance_rows(
            case.input_fixture,
            pairwise_deletion=case.pairwise_deletion,
        )
        return summary, rows, None
    if case.operation == "dna-translation":
        if case.genetic_code_id is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a genetic code id"
            )
        summary, rows = _build_bijux_translation_rows(
            case.input_fixture,
            genetic_code_id=case.genetic_code_id,
        )
        return summary, rows, None
    raise ValueError(f"unsupported ape parity operation '{case.operation}'")


def _load_reference_case_payload(
    case: ApeParityCase,
    execution_root: Path,
) -> tuple[dict[str, object], list[dict[str, object]] | None, str | None]:
    if case.operation in {"read-tree-structure", "write-tree-structure", "root-tree-outgroup"}:
        summary = _normalize_reference_summary(_load_json(execution_root / "summary.json"))
        expected_tip_labels = {
            str(label) for label in summary.get("tip_labels", [])
        }
        rows = _load_rows_table(
            execution_root / "clades.tsv",
            expected_tip_labels=expected_tip_labels,
            sort_rows=True,
        )
        normalized_text = _canonical_newick(
            execution_root / "normalized-tree.nwk",
            expected_tip_labels=expected_tip_labels,
        )
        return summary, rows, normalized_text
    if case.operation in {"read-tree-set-structure", "write-tree-set-structure"}:
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "clades.tsv", sort_rows=True)
        return summary, rows, None
    if case.operation == "dna-base-frequency":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "base-frequency.tsv")
        return summary, rows, None
    if case.operation == "dna-raw-distance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "distance-matrix.tsv")
        return summary, rows, None
    if case.operation == "dna-translation":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "translation.tsv")
        return summary, rows, None
    raise ValueError(f"unsupported ape parity operation '{case.operation}'")


def _tree_structure_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    expected_tree = load_tree(case.input_fixture)
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _tree_set_structure_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree_set = load_newick_tree_set(execution_root / "normalized-tree-set.nwk")
    expected_tree_set = load_newick_tree_set(case.input_fixture)
    expected_tip_labels = {
        tip_name
        for tree in expected_tree_set
        for tip_name in tree.tip_names
    }
    for tree in reference_tree_set:
        _normalize_tree_labels(tree.root, expected_tip_labels=expected_tip_labels)
    report = compare_tree_sets_structurally(
        expected_tree_set,
        reference_tree_set,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _root_tree_outgroup_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    # Canonical Newick does not preserve rootedness metadata, but ape::root
    # produced this record explicitly as a rooted output for these governed cases.
    reference_tree.rooted = True
    expected_tree, _report = root_tree_on_outgroup(
        case.input_fixture,
        outgroup_taxa=list(case.outgroup_taxa),
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _summary_rows(observations: list[ApeParityObservation]) -> list[ApeParitySummaryRow]:
    rows: list[ApeParitySummaryRow] = []
    for function_name in sorted({item.function_name for item in observations}):
        selected = [item for item in observations if item.function_name == function_name]
        rows.append(
            ApeParitySummaryRow(
                function_name=function_name,
                case_count=len(selected),
                passed_case_count=sum(1 for item in selected if item.status == "passed"),
                failed_case_count=sum(1 for item in selected if item.status == "failed"),
                skipped_case_count=sum(1 for item in selected if item.status == "skipped"),
            )
        )
    return rows


def run_ape_parity_cases(
    *,
    case_ids: list[str] | None = None,
    rscript_executable: str = "Rscript",
    failure_root: Path | None = None,
    fixtures_root: Path | None = None,
) -> ApeParityReport:
    """Run governed live `ape` parity cases through the checked-in R runner."""
    selected = _selected_cases(case_ids=case_ids, fixtures_root=fixtures_root)
    observations: list[ApeParityObservation] = []
    active_failure_root = _failure_root() if failure_root is None else failure_root
    bijux_version = _bijux_version()
    bijux_commit = _bijux_commit()
    for case in selected:
        with tempfile.TemporaryDirectory(prefix=f"bijux-ape-parity-{case.case_id}-") as tmpdir:
            working_root = Path(tmpdir)
            reference_input_path = _materialize_reference_input(case, working_root)
            reference_case = replace(case, input_fixture=reference_input_path)
            case_file = _write_case_file(working_root / "case.json", reference_case)
            execution_root = working_root / "reference"
            execution_root.mkdir(parents=True, exist_ok=True)
            bijux_summary: dict[str, object] | None = None
            bijux_rows: list[dict[str, object]] | None = None
            bijux_normalized_text: str | None = None
            bijux_error: dict[str, object] | None = None
            try:
                (
                    bijux_summary,
                    bijux_rows,
                    bijux_normalized_text,
                ) = _build_bijux_case_payload(case)
            except Exception as error:
                bijux_error = {
                    "error_type": type(error).__name__,
                    "message": str(error),
                }
            execution_payload: dict[str, object] | None = None
            reference_summary: dict[str, object] | None = None
            reference_error: dict[str, object] | None = None
            reference_rows: list[dict[str, object]] | None = None
            reference_normalized_text: str | None = None
            status = "failed"
            mismatch_reason: str | None = None
            artifact_root: Path | None = None
            r_version: str | None = None
            ape_version: str | None = None
            process_stdout = ""
            process_stderr = ""
            try:
                process = subprocess.run(
                    [
                        rscript_executable,
                        str(_ape_runner_path()),
                        str(case_file),
                        str(execution_root),
                    ],
                    capture_output=True,
                    check=False,
                    cwd=_repository_root(),
                    env=_reference_environment(),
                    text=True,
                )
                process_stdout = process.stdout
                process_stderr = process.stderr
            except FileNotFoundError:
                process = None
                status = "skipped"
                mismatch_reason = "rscript_unavailable"
            if process is None:
                pass
            elif process.returncode != 0:
                mismatch_reason = "reference_execution_failed"
            else:
                execution_path = execution_root / "reference-execution.json"
                if not execution_path.exists():
                    mismatch_reason = "reference_execution_failed"
                else:
                    execution_payload = _load_json(execution_path)
                    r_version = execution_payload.get("r_version")  # type: ignore[assignment]
                    ape_version = execution_payload.get("ape_version")  # type: ignore[assignment]
                    execution_status = execution_payload.get("status")
                    if execution_status == "unavailable":
                        status = "skipped"
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason", "ape_package_unavailable"
                            )
                        )
                    elif execution_status != "ok":
                        reference_error = {
                            "error_type": str(
                                execution_payload.get(
                                    "error_type",
                                    execution_payload.get(
                                        "mismatch_reason", "reference_execution_failed"
                                    ),
                                )
                            ),
                            "message": str(execution_payload.get("message", "")),
                        }
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason", "reference_execution_failed"
                            )
                        )
                    else:
                        (
                            reference_summary,
                            reference_rows,
                            reference_normalized_text,
                        ) = _load_reference_case_payload(case, execution_root)
                        if case.operation in {"read-tree-structure", "write-tree-structure"}:
                            mismatch_reason = _tree_structure_mismatch_reason(
                                case,
                                execution_root,
                            )
                        elif case.operation == "root-tree-outgroup":
                            mismatch_reason = _root_tree_outgroup_mismatch_reason(
                                case,
                                execution_root,
                            )
                        elif case.operation in {
                            "read-tree-set-structure",
                            "write-tree-set-structure",
                        }:
                            mismatch_reason = _tree_set_structure_mismatch_reason(
                                case,
                                execution_root,
                            )
                        elif not _compare_json(
                            reference_summary, bijux_summary, tolerance=case.tolerance
                        ):
                            mismatch_reason = "summary_mismatch"
                        elif not _compare_json(
                            reference_rows,
                            bijux_rows,
                            tolerance=case.tolerance,
                        ):
                            mismatch_reason = "rows_mismatch"
                        elif reference_normalized_text != bijux_normalized_text:
                            mismatch_reason = "normalized_text_mismatch"
                        else:
                            status = "passed"
                        if mismatch_reason is None:
                            status = "passed"
            if case.expected_status == "parse-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_parse_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_parse_error_missing"
                elif not bijux_error.get("message") or not reference_error.get("message"):
                    mismatch_reason = "parse_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "rooting-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_rooting_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_rooting_error_missing"
                elif not bijux_error.get("message") or not reference_error.get("message"):
                    mismatch_reason = "rooting_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if status != "passed":
                artifact_root = _persist_failure_bundle(
                    failure_root=active_failure_root,
                    case=case,
                    case_file=case_file,
                    execution_root=execution_root,
                    execution_payload=execution_payload,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                    reference_rows=reference_rows,
                    bijux_rows=bijux_rows,
                    bijux_normalized_text=bijux_normalized_text,
                    mismatch_reason=mismatch_reason or "reference_execution_failed",
                )
                if process_stdout:
                    (artifact_root / "reference-stdout.txt").write_text(
                        process_stdout, encoding="utf-8"
                    )
                if process_stderr:
                    (artifact_root / "reference-stderr.txt").write_text(
                        process_stderr, encoding="utf-8"
                    )
            observations.append(
                ApeParityObservation(
                    case_id=case.case_id,
                    fixture_kind=case.fixture_kind,
                    fixture_id=case.fixture_id,
                    function_name=case.function_name,
                    python_function_name=case.python_function_name,
                    input_fixture=case.input_fixture,
                    tolerance=case.tolerance,
                    r_version=r_version,
                    ape_version=ape_version,
                    bijux_version=bijux_version,
                    bijux_commit=bijux_commit,
                    status=status,
                    passed=status == "passed",
                    mismatch_reason=mismatch_reason,
                    reproducible_artifact_root=artifact_root,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                )
            )
    case_count = len(observations)
    passed_case_count = sum(1 for item in observations if item.status == "passed")
    failed_case_count = sum(1 for item in observations if item.status == "failed")
    skipped_case_count = sum(1 for item in observations if item.status == "skipped")
    return ApeParityReport(
        observations=observations,
        summary_rows=_summary_rows(observations),
        case_count=case_count,
        passed_case_count=passed_case_count,
        failed_case_count=failed_case_count,
        skipped_case_count=skipped_case_count,
        all_passed=case_count > 0
        and passed_case_count == case_count
        and failed_case_count == 0
        and skipped_case_count == 0,
        limitations=[
            "The governed live `ape` parity registry is intentionally narrow until later rounds expand the shared fixture surface.",
            "This harness requires Rscript plus the `ape` and `jsonlite` R packages for live reference execution.",
        ],
    )


def write_ape_parity_summary_table(path: Path, report: ApeParityReport) -> Path:
    """Write one row per governed `ape` function summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "function_name",
                "case_count",
                "passed_case_count",
                "failed_case_count",
                "skipped_case_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.summary_rows:
            writer.writerow(asdict(row))
    return path


def write_ape_parity_observation_table(path: Path, report: ApeParityReport) -> Path:
    """Write one row per governed `ape` parity observation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "fixture_kind",
                "fixture_id",
                "function_name",
                "python_function_name",
                "input_fixture",
                "tolerance",
                "r_version",
                "ape_version",
                "bijux_version",
                "bijux_commit",
                "status",
                "passed",
                "mismatch_reason",
                "reproducible_artifact_root",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for observation in report.observations:
            writer.writerow(
                {
                    "case_id": observation.case_id,
                    "fixture_kind": observation.fixture_kind,
                    "fixture_id": observation.fixture_id,
                    "function_name": observation.function_name,
                    "python_function_name": observation.python_function_name,
                    "input_fixture": str(observation.input_fixture),
                    "tolerance": format(observation.tolerance, ".12g"),
                    "r_version": observation.r_version or "",
                    "ape_version": observation.ape_version or "",
                    "bijux_version": observation.bijux_version,
                    "bijux_commit": observation.bijux_commit or "",
                    "status": observation.status,
                    "passed": str(observation.passed).lower(),
                    "mismatch_reason": observation.mismatch_reason or "",
                    "reproducible_artifact_root": ""
                    if observation.reproducible_artifact_root is None
                    else str(observation.reproducible_artifact_root),
                }
            )
    return path
