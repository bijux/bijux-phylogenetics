from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tempfile

from bijux_phylogenetics.bayesian.beast.logs import summarize_beast_log
from bijux_phylogenetics.bayesian.beast.posterior_trees import (
    summarize_beast_posterior_trees,
)
from bijux_phylogenetics.bayesian.mrbayes.posterior_trees import (
    parse_mrbayes_consensus_tree,
    summarize_mrbayes_posterior_trees,
)
from bijux_phylogenetics.bayesian.posterior_sets.tree_sets import (
    summarize_maximum_clade_credibility_tree,
)
from bijux_phylogenetics.compare.topology import (
    compare_branch_score_distance,
    compare_robinson_foulds,
)
from bijux_phylogenetics.datasets.shared_fixtures.beast_posteriors import (
    get_shared_beast_posterior_fixture,
)
from bijux_phylogenetics.io.newick import write_newick

_STRICT_YULE_REAL_POSTERIOR = "strict_yule_real_posterior"
_MRBAYES_DEFAULT_BURNIN_FRACTION = 0.25
_MRBAYES_PARTITIONED_ANALYSIS = "partitioned-analysis"

BAYESIAN_WRAPPER_CORRESPONDENCE_STATUSES = (
    "exact-match",
    "tolerance-match",
    "expected-model-assumption-difference",
    "unsupported-case",
    "native-bug",
)


@dataclass(frozen=True, slots=True)
class BayesianWrapperCorrespondenceObservation:
    """One governed native-versus-wrapper Bayesian correspondence case."""

    case_id: str
    wrapper_engine: str
    native_surface: str
    wrapper_surface: str
    comparison_policy: str
    status: str
    supported: bool
    blocking: bool
    tolerance: float | None
    rationale: str
    input_fixtures: list[Path]
    expected_output: dict[str, object]
    observed_output: dict[str, object]


@dataclass(frozen=True, slots=True)
class BayesianWrapperCorrespondenceSummaryRow:
    """One status-level summary row across governed correspondence cases."""

    status: str
    case_count: int
    blocking_case_count: int
    case_ids: list[str]


@dataclass(frozen=True, slots=True)
class BayesianWrapperCorrespondenceReport:
    """One governed correspondence report across cached Bayesian wrapper fixtures."""

    observations: list[BayesianWrapperCorrespondenceObservation]
    summary_rows: list[BayesianWrapperCorrespondenceSummaryRow]
    case_count: int
    supported_case_count: int
    exact_match_case_count: int
    tolerance_match_case_count: int
    expected_model_assumption_difference_case_count: int
    unsupported_case_count: int
    native_bug_case_count: int
    blocking_case_count: int
    all_supported_cases_clear: bool


def summarize_bayesian_wrapper_correspondence() -> BayesianWrapperCorrespondenceReport:
    """Summarize governed correspondence between native Bayesian reports and wrapper artifacts."""
    observations = [
        _build_beast_log_parameter_summary_observation(),
        _build_beast_consensus_tree_observation(),
        _build_beast_maximum_clade_credibility_observation(),
        _build_mrbayes_consensus_topology_observation(),
        _build_mrbayes_consensus_branch_length_semantics_observation(),
        _build_revbayes_unsupported_observation(),
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
    return BayesianWrapperCorrespondenceReport(
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


def _build_beast_log_parameter_summary_observation() -> (
    BayesianWrapperCorrespondenceObservation
):
    fixture = get_shared_beast_posterior_fixture(_STRICT_YULE_REAL_POSTERIOR)
    reference = fixture.load_reference()
    observed = summarize_beast_log(
        fixture.posterior_log_path,
        burnin_fraction=fixture.recommended_burnin_fraction,
    )
    expected_output = _flatten_beast_parameter_reference(reference)
    observed_output = _select_expected_parameter_keys(
        _flatten_beast_log_summary(observed),
        expected_output,
    )
    tolerance = 1e-9
    status = (
        "tolerance-match"
        if _numeric_outputs_match(
            expected_output,
            observed_output,
            tolerance=tolerance,
        )
        else "native-bug"
    )
    return BayesianWrapperCorrespondenceObservation(
        case_id="beast-log-parameter-summaries-strict-yule",
        wrapper_engine="BEAST",
        native_surface="summarize_beast_log",
        wrapper_surface="governed BEAST posterior parameter reference",
        comparison_policy="numeric-tolerance",
        status=status,
        supported=True,
        blocking=status == "native-bug",
        tolerance=tolerance,
        rationale=(
            "The governed BEAST posterior log is accepted under numeric tolerance policy; "
            "this checked-in fixture currently matches the cached reference exactly to recorded precision."
        ),
        input_fixtures=[fixture.posterior_log_path, fixture.reference_json_path],
        expected_output=expected_output,
        observed_output=observed_output,
    )


def _build_beast_consensus_tree_observation() -> (
    BayesianWrapperCorrespondenceObservation
):
    fixture = get_shared_beast_posterior_fixture(_STRICT_YULE_REAL_POSTERIOR)
    reference = fixture.load_reference()
    _consensus_tree, observed = summarize_beast_posterior_trees(
        fixture.posterior_trees_path,
        burnin_fraction=fixture.recommended_burnin_fraction,
    )
    expected_output = {
        "consensus_newick": reference.consensus_reference.newick,
        "annotated_node_count": reference.consensus_reference.annotated_node_count,
        "minimum_posterior_probability": (
            reference.consensus_reference.minimum_posterior_probability
        ),
        "maximum_posterior_probability": (
            reference.consensus_reference.maximum_posterior_probability
        ),
    }
    observed_output = {
        "consensus_newick": observed.consensus_newick,
        "annotated_node_count": observed.annotated_node_count,
        "minimum_posterior_probability": observed.minimum_posterior_probability,
        "maximum_posterior_probability": observed.maximum_posterior_probability,
    }
    status = "exact-match" if observed_output == expected_output else "native-bug"
    return BayesianWrapperCorrespondenceObservation(
        case_id="beast-consensus-tree-strict-yule",
        wrapper_engine="BEAST",
        native_surface="summarize_beast_posterior_trees",
        wrapper_surface="checked-in BEAST consensus tree",
        comparison_policy="exact",
        status=status,
        supported=True,
        blocking=status == "native-bug",
        tolerance=None,
        rationale=(
            "The governed BEAST consensus tree is a directly comparable majority-rule posterior summary and should match exactly."
        ),
        input_fixtures=[
            fixture.posterior_trees_path,
            fixture.consensus_tree_path,
            fixture.reference_json_path,
        ],
        expected_output=expected_output,
        observed_output=observed_output,
    )


def _build_beast_maximum_clade_credibility_observation() -> (
    BayesianWrapperCorrespondenceObservation
):
    fixture = get_shared_beast_posterior_fixture(_STRICT_YULE_REAL_POSTERIOR)
    reference = fixture.load_reference()
    _mcc_tree, observed = summarize_maximum_clade_credibility_tree(
        fixture.posterior_trees_path,
        burnin_fraction=fixture.recommended_burnin_fraction,
    )
    expected_output = {
        "selected_tree_index": reference.mcc_reference.selected_tree_index,
        "clade_credibility_score": reference.mcc_reference.clade_credibility_score,
        "mcc_newick": reference.mcc_reference.newick,
    }
    observed_output = {
        "selected_tree_index": observed.selected_tree_index,
        "clade_credibility_score": observed.clade_credibility_score,
        "mcc_newick": observed.mcc_newick,
    }
    status = "exact-match" if observed_output == expected_output else "native-bug"
    return BayesianWrapperCorrespondenceObservation(
        case_id="beast-maximum-clade-credibility-tree-strict-yule",
        wrapper_engine="BEAST",
        native_surface="summarize_maximum_clade_credibility_tree",
        wrapper_surface="checked-in BEAST maximum clade credibility tree",
        comparison_policy="exact",
        status=status,
        supported=True,
        blocking=status == "native-bug",
        tolerance=None,
        rationale=(
            "The governed BEAST MCC tree is a directly comparable posterior tree-set selection and should match exactly."
        ),
        input_fixtures=[
            fixture.posterior_trees_path,
            fixture.mcc_tree_path,
            fixture.reference_json_path,
        ],
        expected_output=expected_output,
        observed_output=observed_output,
    )


def _build_mrbayes_consensus_topology_observation() -> (
    BayesianWrapperCorrespondenceObservation
):
    tree_set_path, consensus_path = _mrbayes_partitioned_fixture_paths()
    native_tree, _native_summary = summarize_mrbayes_posterior_trees(
        tree_set_path,
        burnin_fraction=_MRBAYES_DEFAULT_BURNIN_FRACTION,
    )
    wrapper_tree, wrapper_report = parse_mrbayes_consensus_tree(consensus_path)
    topology_report = _compare_tree_topology(native_tree, wrapper_tree)
    expected_output = {
        "rooted_robinson_foulds_distance": 0,
        "topology_equal": True,
    }
    observed_output = {
        "rooted_robinson_foulds_distance": topology_report.robinson_foulds_distance,
        "normalized_robinson_foulds": topology_report.normalized_robinson_foulds,
        "topology_equal": topology_report.topology_equal,
        "annotated_node_count": wrapper_report.annotated_node_count,
    }
    status = (
        "exact-match"
        if observed_output["rooted_robinson_foulds_distance"] == 0
        and bool(observed_output["topology_equal"])
        else "native-bug"
    )
    return BayesianWrapperCorrespondenceObservation(
        case_id="mrbayes-consensus-topology-partitioned-analysis",
        wrapper_engine="MrBayes",
        native_surface="summarize_mrbayes_posterior_trees",
        wrapper_surface="checked-in MrBayes consensus tree topology",
        comparison_policy="exact",
        status=status,
        supported=True,
        blocking=status == "native-bug",
        tolerance=None,
        rationale=(
            "The governed MrBayes consensus topology is directly comparable to the native majority-rule posterior topology."
        ),
        input_fixtures=[tree_set_path, consensus_path],
        expected_output=expected_output,
        observed_output=observed_output,
    )


def _build_mrbayes_consensus_branch_length_semantics_observation() -> (
    BayesianWrapperCorrespondenceObservation
):
    tree_set_path, consensus_path = _mrbayes_partitioned_fixture_paths()
    native_tree, native_summary = summarize_mrbayes_posterior_trees(
        tree_set_path,
        burnin_fraction=_MRBAYES_DEFAULT_BURNIN_FRACTION,
    )
    wrapper_tree, wrapper_report = parse_mrbayes_consensus_tree(consensus_path)
    topology_report = _compare_tree_topology(native_tree, wrapper_tree)
    branch_report = _compare_tree_branch_score(native_tree, wrapper_tree)
    expected_output = {
        "rooted_robinson_foulds_distance": 0,
        "same_topology_different_branch_lengths": True,
    }
    observed_output = {
        "rooted_robinson_foulds_distance": topology_report.robinson_foulds_distance,
        "branch_score_distance": branch_report.branch_score_distance,
        "same_topology_different_branch_lengths": (
            topology_report.robinson_foulds_distance == 0
            and branch_report.branch_score_distance is not None
            and branch_report.branch_score_distance > 0.0
        ),
        "native_consensus_newick": native_summary.consensus_newick,
        "wrapper_consensus_newick": wrapper_report.consensus_newick,
    }
    status = (
        "expected-model-assumption-difference"
        if observed_output["same_topology_different_branch_lengths"]
        else "native-bug"
    )
    return BayesianWrapperCorrespondenceObservation(
        case_id="mrbayes-consensus-branch-length-semantics-partitioned-analysis",
        wrapper_engine="MrBayes",
        native_surface="summarize_mrbayes_posterior_trees",
        wrapper_surface="checked-in MrBayes consensus branch lengths and node annotations",
        comparison_policy="expected-difference",
        status=status,
        supported=True,
        blocking=status == "native-bug",
        tolerance=None,
        rationale=(
            "MrBayes consensus topology is comparable, but its branch-length and annotation semantics differ from the native plain majority-rule summary and should be recorded as an expected semantic difference."
        ),
        input_fixtures=[tree_set_path, consensus_path],
        expected_output=expected_output,
        observed_output=observed_output,
    )


def _build_revbayes_unsupported_observation() -> (
    BayesianWrapperCorrespondenceObservation
):
    return BayesianWrapperCorrespondenceObservation(
        case_id="revbayes-governed-posterior-corpus",
        wrapper_engine="RevBayes-style",
        native_surface="none",
        wrapper_surface="governed RevBayes-style posterior fixture",
        comparison_policy="unsupported",
        status="unsupported-case",
        supported=False,
        blocking=False,
        tolerance=None,
        rationale=(
            "The repository does not yet carry a governed RevBayes-style posterior corpus or a native RevBayes-wrapper parser surface, so the correspondence lane must stay explicitly unsupported."
        ),
        input_fixtures=[],
        expected_output={"governed_fixture_available": False},
        observed_output={"native_wrapper_surface_available": False},
    )


def _build_summary_rows(
    observations: list[BayesianWrapperCorrespondenceObservation],
) -> list[BayesianWrapperCorrespondenceSummaryRow]:
    rows: list[BayesianWrapperCorrespondenceSummaryRow] = []
    for status in BAYESIAN_WRAPPER_CORRESPONDENCE_STATUSES:
        selected = [
            observation for observation in observations if observation.status == status
        ]
        rows.append(
            BayesianWrapperCorrespondenceSummaryRow(
                status=status,
                case_count=len(selected),
                blocking_case_count=sum(
                    1 for observation in selected if observation.blocking
                ),
                case_ids=[observation.case_id for observation in selected],
            )
        )
    return rows


def _flatten_beast_parameter_reference(reference) -> dict[str, float]:
    flattened: dict[str, float] = {}
    for parameter_name, summary in sorted(reference.parameter_reference.items()):
        flattened[f"{parameter_name}.effective_sample_size"] = (
            summary.effective_sample_size
        )
        flattened[f"{parameter_name}.mean"] = summary.mean
        flattened[f"{parameter_name}.median"] = summary.median
        flattened[f"{parameter_name}.hpd_95_lower"] = summary.hpd_95_lower
        flattened[f"{parameter_name}.hpd_95_upper"] = summary.hpd_95_upper
    return flattened


def _flatten_beast_log_summary(summary) -> dict[str, float]:
    flattened: dict[str, float] = {}
    for row in summary.parameter_summaries:
        flattened[f"{row.parameter}.effective_sample_size"] = row.effective_sample_size
        flattened[f"{row.parameter}.mean"] = row.mean
        flattened[f"{row.parameter}.median"] = row.median
        flattened[f"{row.parameter}.hpd_95_lower"] = row.hpd_95_lower
        flattened[f"{row.parameter}.hpd_95_upper"] = row.hpd_95_upper
    return flattened


def _numeric_outputs_match(
    expected_output: dict[str, float],
    observed_output: dict[str, float],
    *,
    tolerance: float,
) -> bool:
    if set(expected_output) != set(observed_output):
        return False
    return all(
        abs(observed_output[key] - expected_output[key]) <= tolerance
        for key in expected_output
    )


def _select_expected_parameter_keys(
    observed_output: dict[str, float],
    expected_output: dict[str, float],
) -> dict[str, float]:
    return {key: observed_output[key] for key in expected_output}


def _mrbayes_partitioned_fixture_paths() -> tuple[Path, Path]:
    base = _package_root() / "tests" / "fixtures" / "mrbayes"
    return (
        base / f"{_MRBAYES_PARTITIONED_ANALYSIS}.run1.t",
        base / f"{_MRBAYES_PARTITIONED_ANALYSIS}.con.tre",
    )


def _package_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _compare_tree_topology(left_tree, right_tree):
    with tempfile.TemporaryDirectory(
        prefix="bijux-bayesian-wrapper-topology-"
    ) as tmp_dir:
        tmp_root = Path(tmp_dir)
        left_path = write_newick(tmp_root / "left.nwk", left_tree)
        right_path = write_newick(tmp_root / "right.nwk", right_tree)
        return compare_robinson_foulds(left_path, right_path, rf_mode="rooted")


def _compare_tree_branch_score(left_tree, right_tree):
    with tempfile.TemporaryDirectory(
        prefix="bijux-bayesian-wrapper-branch-score-"
    ) as tmp_dir:
        tmp_root = Path(tmp_dir)
        left_path = write_newick(tmp_root / "left.nwk", left_tree)
        right_path = write_newick(tmp_root / "right.nwk", right_tree)
        return compare_branch_score_distance(left_path, right_path)
