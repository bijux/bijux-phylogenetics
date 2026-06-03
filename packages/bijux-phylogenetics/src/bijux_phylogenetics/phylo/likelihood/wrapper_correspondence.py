from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile

from bijux_phylogenetics.compare.topology import (
    compare_branch_score_distance,
    compare_robinson_foulds,
)
from bijux_phylogenetics.io.newick import load_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.likelihood.models import (
    LikelihoodWrapperCorrespondenceObservation,
    LikelihoodWrapperCorrespondenceReport,
    LikelihoodWrapperCorrespondenceSummaryRow,
)
from bijux_phylogenetics.phylo.likelihood.tree_inference import (
    infer_nucleotide_likelihood_tree_from_alignment,
)
from bijux_phylogenetics.phylo.topology import rooted_topology_fingerprint
from bijux_phylogenetics.reports.service.artifacts import write_json_artifact

LIKELIHOOD_WRAPPER_CORRESPONDENCE_STATUSES = (
    "exact-match",
    "tolerance-match",
    "expected-model-assumption-difference",
    "unsupported-case",
    "native-bug",
)
_IQTREE_BRANCH_SCORE_TOLERANCE = 1e-6


def summarize_likelihood_wrapper_correspondence() -> (
    LikelihoodWrapperCorrespondenceReport
):
    """Summarize governed correspondence between native ML inference and wrapper references."""
    native_report = infer_nucleotide_likelihood_tree_from_alignment(
        _alignment_fixture_path(),
        model_name="jc69",
        search_method="nni",
        start_tree_count=3,
        start_tree_seed=17,
        upper_branch_length_bound=1.0,
    )
    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary_root = Path(temporary_directory)
        native_tree_path = temporary_root / "native_ml_tree.nwk"
        write_newick(
            native_tree_path, loads_newick(native_report.best_final_tree_newick)
        )
        observations = [
            _build_iqtree_topology_exact_observation(
                native_report=native_report,
                native_tree_path=native_tree_path,
            ),
            _build_iqtree_branch_length_tolerance_observation(
                native_report=native_report,
                native_tree_path=native_tree_path,
            ),
            _build_fasttree_expected_difference_observation(
                native_report=native_report,
                native_tree_path=native_tree_path,
            ),
            _build_raxml_unsupported_observation(),
        ]
    summary_rows = _build_summary_rows(observations)
    exact_match_case_count = sum(
        1 for observation in observations if observation.status == "exact-match"
    )
    tolerance_match_case_count = sum(
        1 for observation in observations if observation.status == "tolerance-match"
    )
    expected_model_assumption_difference_case_count = sum(
        1
        for observation in observations
        if observation.status == "expected-model-assumption-difference"
    )
    unsupported_case_count = sum(
        1 for observation in observations if observation.status == "unsupported-case"
    )
    native_bug_case_count = sum(
        1 for observation in observations if observation.status == "native-bug"
    )
    blocking_case_count = sum(1 for observation in observations if observation.blocking)
    supported_case_count = sum(
        1 for observation in observations if observation.supported
    )
    return LikelihoodWrapperCorrespondenceReport(
        observations=observations,
        summary_rows=summary_rows,
        case_count=len(observations),
        supported_case_count=supported_case_count,
        exact_match_case_count=exact_match_case_count,
        tolerance_match_case_count=tolerance_match_case_count,
        expected_model_assumption_difference_case_count=(
            expected_model_assumption_difference_case_count
        ),
        unsupported_case_count=unsupported_case_count,
        native_bug_case_count=native_bug_case_count,
        blocking_case_count=blocking_case_count,
        all_supported_cases_clear=blocking_case_count == 0,
    )


def write_likelihood_wrapper_correspondence_summary_table(
    path: Path,
    report: LikelihoodWrapperCorrespondenceReport,
) -> Path:
    """Write one status-level ML wrapper correspondence summary table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "status",
                "case_count",
                "blocking_case_count",
                "case_ids",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.summary_rows:
            writer.writerow(
                {
                    "status": row.status,
                    "case_count": row.case_count,
                    "blocking_case_count": row.blocking_case_count,
                    "case_ids": "|".join(row.case_ids),
                }
            )
    return path


def write_likelihood_wrapper_correspondence_observation_table(
    path: Path,
    report: LikelihoodWrapperCorrespondenceReport,
) -> Path:
    """Write one case-level ML wrapper correspondence ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "wrapper_engine",
                "native_surface",
                "wrapper_surface",
                "comparison_policy",
                "status",
                "supported",
                "blocking",
                "tolerance",
                "rationale",
                "input_fixtures",
                "expected_output",
                "observed_output",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for observation in report.observations:
            writer.writerow(
                {
                    "case_id": observation.case_id,
                    "wrapper_engine": observation.wrapper_engine,
                    "native_surface": observation.native_surface,
                    "wrapper_surface": observation.wrapper_surface,
                    "comparison_policy": observation.comparison_policy,
                    "status": observation.status,
                    "supported": str(observation.supported).lower(),
                    "blocking": str(observation.blocking).lower(),
                    "tolerance": (
                        ""
                        if observation.tolerance is None
                        else format(observation.tolerance, ".15g")
                    ),
                    "rationale": observation.rationale,
                    "input_fixtures": "|".join(observation.input_fixtures),
                    "expected_output": json.dumps(
                        observation.expected_output,
                        sort_keys=True,
                    ),
                    "observed_output": json.dumps(
                        observation.observed_output,
                        sort_keys=True,
                    ),
                }
            )
    return path


def write_likelihood_wrapper_correspondence_run_json(
    path: Path,
    report: LikelihoodWrapperCorrespondenceReport,
) -> Path:
    """Write one machine-readable ML wrapper correspondence payload."""
    return write_json_artifact(
        path,
        {
            "case_count": report.case_count,
            "supported_case_count": report.supported_case_count,
            "exact_match_case_count": report.exact_match_case_count,
            "tolerance_match_case_count": report.tolerance_match_case_count,
            "expected_model_assumption_difference_case_count": (
                report.expected_model_assumption_difference_case_count
            ),
            "unsupported_case_count": report.unsupported_case_count,
            "native_bug_case_count": report.native_bug_case_count,
            "blocking_case_count": report.blocking_case_count,
            "all_supported_cases_clear": report.all_supported_cases_clear,
            "summary_rows": [
                {
                    "status": row.status,
                    "case_count": row.case_count,
                    "blocking_case_count": row.blocking_case_count,
                    "case_ids": row.case_ids,
                }
                for row in report.summary_rows
            ],
            "observations": [
                {
                    "case_id": observation.case_id,
                    "wrapper_engine": observation.wrapper_engine,
                    "native_surface": observation.native_surface,
                    "wrapper_surface": observation.wrapper_surface,
                    "comparison_policy": observation.comparison_policy,
                    "status": observation.status,
                    "supported": observation.supported,
                    "blocking": observation.blocking,
                    "tolerance": observation.tolerance,
                    "rationale": observation.rationale,
                    "input_fixtures": observation.input_fixtures,
                    "expected_output": observation.expected_output,
                    "observed_output": observation.observed_output,
                }
                for observation in report.observations
            ],
        },
    )


def write_likelihood_wrapper_correspondence_artifacts(
    out_dir: Path,
    report: LikelihoodWrapperCorrespondenceReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one ML wrapper correspondence run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = write_likelihood_wrapper_correspondence_summary_table(
        out_dir / "summary.tsv",
        report,
    )
    observations_path = write_likelihood_wrapper_correspondence_observation_table(
        out_dir / "observations.tsv",
        report,
    )
    run_json_path = write_likelihood_wrapper_correspondence_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "summary_path": summary_path,
        "observations_path": observations_path,
        "run_json_path": run_json_path,
    }


def _build_iqtree_topology_exact_observation(
    *,
    native_report,
    native_tree_path: Path,
) -> LikelihoodWrapperCorrespondenceObservation:
    wrapper_tree_path = _wrapper_fixture_root() / "iqtree_small_topology.treefile"
    wrapper_tree = load_newick(wrapper_tree_path)
    wrapper_topology_fingerprint = rooted_topology_fingerprint(wrapper_tree)
    rf_report = compare_robinson_foulds(native_tree_path, wrapper_tree_path)
    status = "exact-match" if rf_report.topology_equal else "native-bug"
    return LikelihoodWrapperCorrespondenceObservation(
        case_id="iqtree-small-topology-reference",
        wrapper_engine="IQ-TREE-style",
        native_surface="infer_nucleotide_likelihood_tree_from_alignment",
        wrapper_surface="checked-in IQ-TREE-style ML topology reference",
        comparison_policy="exact-topology-fingerprint",
        status=status,
        supported=True,
        blocking=status == "native-bug",
        tolerance=None,
        rationale=(
            "The governed IQ-TREE-style reference tree should match the native tiny-alignment ML topology exactly."
        ),
        input_fixtures=[
            str(_alignment_fixture_path()),
            str(wrapper_tree_path),
        ],
        expected_output={
            "wrapper_tree_newick": wrapper_tree.to_newick(),
            "wrapper_topology_fingerprint": wrapper_topology_fingerprint,
        },
        observed_output={
            "native_best_tree_newick": native_report.best_final_tree_newick,
            "native_best_topology_fingerprint": (
                native_report.best_final_topology_fingerprint
            ),
            "robinson_foulds_distance": rf_report.robinson_foulds_distance,
            "normalized_robinson_foulds": rf_report.normalized_robinson_foulds,
            "topology_equal": rf_report.topology_equal,
        },
    )


def _build_iqtree_branch_length_tolerance_observation(
    *,
    native_report,
    native_tree_path: Path,
) -> LikelihoodWrapperCorrespondenceObservation:
    wrapper_tree_path = _wrapper_fixture_root() / "iqtree_small_rounded.treefile"
    wrapper_tree = load_newick(wrapper_tree_path)
    branch_score_report = compare_branch_score_distance(
        native_tree_path, wrapper_tree_path
    )
    wrapper_topology_fingerprint = rooted_topology_fingerprint(wrapper_tree)
    topology_equal = (
        native_report.best_final_topology_fingerprint == wrapper_topology_fingerprint
    )
    branch_score_distance = branch_score_report.branch_score_distance
    status = (
        "tolerance-match"
        if topology_equal
        and branch_score_distance is not None
        and branch_score_distance <= _IQTREE_BRANCH_SCORE_TOLERANCE
        else "native-bug"
    )
    return LikelihoodWrapperCorrespondenceObservation(
        case_id="iqtree-small-rounded-branch-reference",
        wrapper_engine="IQ-TREE-style",
        native_surface="infer_nucleotide_likelihood_tree_from_alignment",
        wrapper_surface="checked-in IQ-TREE-style rounded branch-length reference",
        comparison_policy="branch-score-tolerance",
        status=status,
        supported=True,
        blocking=status == "native-bug",
        tolerance=_IQTREE_BRANCH_SCORE_TOLERANCE,
        rationale=(
            "The governed IQ-TREE-style rounded branch-length reference intentionally rounds numerically negligible leaves and should stay within the declared branch-score tolerance."
        ),
        input_fixtures=[
            str(_alignment_fixture_path()),
            str(wrapper_tree_path),
        ],
        expected_output={
            "wrapper_tree_newick": wrapper_tree.to_newick(),
            "wrapper_topology_fingerprint": wrapper_topology_fingerprint,
        },
        observed_output={
            "native_best_tree_newick": native_report.best_final_tree_newick,
            "native_best_topology_fingerprint": (
                native_report.best_final_topology_fingerprint
            ),
            "branch_score_distance": branch_score_distance,
            "shared_split_count": branch_score_report.shared_split_count,
            "topology_equal": topology_equal,
        },
    )


def _build_fasttree_expected_difference_observation(
    *,
    native_report,
    native_tree_path: Path,
) -> LikelihoodWrapperCorrespondenceObservation:
    wrapper_tree_path = _wrapper_fixture_root() / "fasttree_small_approximate.nwk"
    wrapper_tree = load_newick(wrapper_tree_path)
    wrapper_topology_fingerprint = rooted_topology_fingerprint(wrapper_tree)
    rf_report = compare_robinson_foulds(native_tree_path, wrapper_tree_path)
    branch_score_report = compare_branch_score_distance(
        native_tree_path, wrapper_tree_path
    )
    return LikelihoodWrapperCorrespondenceObservation(
        case_id="fasttree-small-approximate-reference",
        wrapper_engine="FastTree-style",
        native_surface="infer_nucleotide_likelihood_tree_from_alignment",
        wrapper_surface="checked-in FastTree-style approximate ML tree reference",
        comparison_policy="expected-difference",
        status="expected-model-assumption-difference",
        supported=True,
        blocking=False,
        tolerance=None,
        rationale=(
            "FastTree-style approximate inference can share the same tiny-alignment topology while still differing in optimization and SH-like local support semantics, so correspondence must remain an expected assumption difference rather than exact parity."
        ),
        input_fixtures=[
            str(_alignment_fixture_path()),
            str(wrapper_tree_path),
        ],
        expected_output={
            "wrapper_tree_newick": wrapper_tree.to_newick(),
            "wrapper_topology_fingerprint": wrapper_topology_fingerprint,
            "support_labels_present": True,
        },
        observed_output={
            "native_best_tree_newick": native_report.best_final_tree_newick,
            "native_best_topology_fingerprint": (
                native_report.best_final_topology_fingerprint
            ),
            "robinson_foulds_distance": rf_report.robinson_foulds_distance,
            "normalized_robinson_foulds": rf_report.normalized_robinson_foulds,
            "branch_score_distance": branch_score_report.branch_score_distance,
            "topology_equal": rf_report.topology_equal,
        },
    )


def _build_raxml_unsupported_observation() -> (
    LikelihoodWrapperCorrespondenceObservation
):
    return LikelihoodWrapperCorrespondenceObservation(
        case_id="raxml-style-small-reference",
        wrapper_engine="RAxML-style",
        native_surface="none",
        wrapper_surface="governed RAxML-style small ML reference",
        comparison_policy="unsupported",
        status="unsupported-case",
        supported=False,
        blocking=False,
        tolerance=None,
        rationale=(
            "The repository does not yet carry a governed small RAxML-style ML reference corpus or a native correspondence normalization surface for it, so this lane must remain explicitly unsupported."
        ),
        input_fixtures=[],
        expected_output={"governed_fixture_available": False},
        observed_output={"native_wrapper_surface_available": False},
    )


def _build_summary_rows(
    observations: list[LikelihoodWrapperCorrespondenceObservation],
) -> list[LikelihoodWrapperCorrespondenceSummaryRow]:
    rows: list[LikelihoodWrapperCorrespondenceSummaryRow] = []
    for status in LIKELIHOOD_WRAPPER_CORRESPONDENCE_STATUSES:
        selected = [
            observation for observation in observations if observation.status == status
        ]
        rows.append(
            LikelihoodWrapperCorrespondenceSummaryRow(
                status=status,
                case_count=len(selected),
                blocking_case_count=sum(
                    1 for observation in selected if observation.blocking
                ),
                case_ids=[observation.case_id for observation in selected],
            )
        )
    return rows


def _alignment_fixture_path() -> Path:
    return _fixture_root() / "alignments" / "jc69_likelihood_nni_alignment_4_taxa.fasta"


def _wrapper_fixture_root() -> Path:
    return _fixture_root() / "engine_outputs" / "likelihood_wrapper_correspondence"


def _fixture_root() -> Path:
    return Path(__file__).resolve().parents[4] / "tests" / "fixtures"
