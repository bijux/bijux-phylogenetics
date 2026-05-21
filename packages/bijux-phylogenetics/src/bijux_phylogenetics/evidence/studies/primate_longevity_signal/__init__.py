from __future__ import annotations

from datetime import date
import json
from pathlib import Path

from .definitions import CLAIM_DEFINITIONS
from .definitions import EVIDENCE_ID
from .definitions import FAMILY_DEFINITIONS
from .definitions import FRAGMENT_CLASSIFICATIONS
from .definitions import SOURCE_LOCATOR
from .definitions import STUDY_ID
from .registry import build_primate_claim_registry
from .registry import build_primate_family_index
from .registry import build_primate_source_fragment_map

def build_primate_summary_bundle_claims(repo_root: Path) -> dict[str, object]:
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


def refresh_primate_summary_bundle(repo_root: Path) -> dict[str, object]:
    bundle_root = _bundle_root(Path(repo_root))
    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("expected summary bundle manifest to be a JSON object")
    claims_payload = build_primate_summary_bundle_claims(Path(repo_root))
    manifest["claim_ids"] = [
        claim["claim_id"]
        for claim in claims_payload["claims"]
        if isinstance(claim, dict) and isinstance(claim.get("claim_id"), str)
    ]
    freshness = manifest.get("freshness")
    if isinstance(freshness, dict):
        freshness["last_generated_on"] = date.today().isoformat()
    return {"manifest": manifest, "claims": claims_payload}

def build_primate_parity_policy(repo_root: Path) -> dict[str, object]:
    fragment_map = build_primate_source_fragment_map(repo_root)
    policies: list[dict[str, object]] = []
    for family_id in [
        "workflow-contracts",
        "data-preparation",
        "tree-operations",
        "visual-surfaces",
        "simulation-inputs",
        "comparative-signal",
        "ancestral-reconstruction",
        "artifact-provenance",
    ]:
        family_fragments = [
            fragment
            for fragment in fragment_map["fragments"]
            if fragment["concept_family"] == family_id
        ]
        if family_id in {
            "workflow-contracts",
            "visual-surfaces",
            "artifact-provenance",
        }:
            expectation = "not_comparable"
            metric_tolerances = []
            rule = "Tracked for reviewer transparency without a numerical parity claim."
        elif family_id in {"data-preparation", "tree-operations", "simulation-inputs"}:
            expectation = "exact"
            metric_tolerances = [{"metric_kind": "default", "tolerance_abs_diff": 0.0}]
            rule = (
                "Scalar counts, boolean checks, and frozen inputs must match exactly."
            )
        elif family_id == "comparative-signal":
            expectation = "statistical_tolerance"
            metric_tolerances = [
                {"metric_kind": "lambda_value", "tolerance_abs_diff": 0.001},
                {"metric_kind": "log_likelihood", "tolerance_abs_diff": 0.001},
                {"metric_kind": "likelihood_ratio", "tolerance_abs_diff": 0.0001},
                {"metric_kind": "p_value", "tolerance_abs_diff": 1e-12},
                {"metric_kind": "covariance_entry", "tolerance_abs_diff": 1e-8},
            ]
            rule = "Likelihood-based metrics must stay within governed floating-point tolerance while preserving the same scientific conclusion."
        else:
            expectation = "near_exact"
            metric_tolerances = [
                {"metric_kind": "point_estimate", "tolerance_abs_diff": 1e-9},
                {"metric_kind": "interval_bound", "tolerance_abs_diff": 1e-9},
                {"metric_kind": "count", "tolerance_abs_diff": 0.0},
            ]
            rule = "Internal-node estimates may differ only by floating-point noise; derived counts must still match exactly."

        policies.append(
            {
                "family_id": family_id,
                "family_title": FAMILY_DEFINITIONS[family_id]["title"],
                "parity_expectation": expectation,
                "metric_tolerances": metric_tolerances,
                "source_fragments": [
                    fragment["fragment_id"] for fragment in family_fragments
                ],
                "rule": rule,
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": EVIDENCE_ID,
        "policy_count": len(policies),
        "policies": policies,
    }


def _observed_abs_diff(left: object, right: object) -> float:
    if isinstance(left, bool) and isinstance(right, bool):
        return 0.0 if left is right else 1.0
    return abs(float(left) - float(right))


def _row(
    *,
    row_id: str,
    fragment_id: str,
    method_family: str,
    metric_name: str,
    parity_expectation: str,
    r_value: object | None,
    bijux_value: object | None,
    tolerance_abs_diff: float,
    observed_abs_diff: float | None = None,
) -> dict[str, object]:
    diff = (
        observed_abs_diff
        if observed_abs_diff is not None
        else _observed_abs_diff(r_value, bijux_value)
    )
    verdict = "matched" if diff == 0.0 else "matched_with_tolerance"
    return {
        "row_id": row_id,
        "fragment_id": fragment_id,
        "method_family": method_family,
        "metric_name": metric_name,
        "parity_expectation": parity_expectation,
        "r_value": r_value,
        "bijux_value": bijux_value,
        "observed_abs_diff": diff,
        "tolerance_abs_diff": tolerance_abs_diff,
        "verdict": verdict if diff <= tolerance_abs_diff else "mismatch_unexplained",
    }


def build_primate_scalar_parity_table(repo_root: Path) -> dict[str, object]:
    bundle_root = _bundle_root(repo_root)
    payload_root = bundle_root / "block-payloads"

    def load_payload(name: str) -> dict[str, object]:
        return json.loads((payload_root / f"{name}.json").read_text(encoding="utf-8"))[
            "evidence"
        ]

    tree_import = load_payload("tree-import-and-pruning")
    unrooted = load_payload("unrooted-tree-demo")
    extract_clade = load_payload("extract-clade-node-77")
    vector_assembly = load_payload("primate-longevity-vector-assembly")
    tree_join = load_payload("treeio-node-mapping-and-join")
    random_lambda = load_payload("random-signal-lambda-fits")
    primate_lambda = load_payload("primate-lambda-fit")
    lambda_zero = load_payload("lambda-zero-covariance-and-lrt")
    ancestral_points = load_payload("continuous-ancestral-point-estimates")
    ancestral_intervals = load_payload("continuous-ancestral-intervals")
    mrca = load_payload("bonobo-gibbon-mrca-estimate")
    increase_counts = load_payload("lifespan-increase-counts")

    rows: list[dict[str, object]] = [
        _row(
            row_id="tree-import-original-tip-count",
            fragment_id="tree-import-and-pruning",
            method_family="tree-operations",
            metric_name="original_tip_count",
            parity_expectation="exact",
            r_value=tree_import["r"]["original_tip_count"],
            bijux_value=tree_import["bijux"]["original_tree"]["inspect"]["tip_count"],
            tolerance_abs_diff=0.0,
        ),
        _row(
            row_id="tree-import-trimmed-tip-count",
            fragment_id="tree-import-and-pruning",
            method_family="tree-operations",
            metric_name="trimmed_tip_count",
            parity_expectation="exact",
            r_value=tree_import["r"]["trimmed_tip_count"],
            bijux_value=tree_import["bijux"]["trimmed_tree"]["inspect"]["tip_count"],
            tolerance_abs_diff=0.0,
        ),
        _row(
            row_id="unrooted-tip-count",
            fragment_id="unrooted-tree-demo",
            method_family="tree-operations",
            metric_name="tip_count",
            parity_expectation="exact",
            r_value=unrooted["r"]["tip_count"],
            bijux_value=unrooted["bijux"]["tip_count"],
            tolerance_abs_diff=0.0,
        ),
        _row(
            row_id="extract-clade-tip-count",
            fragment_id="extract-clade-node-77",
            method_family="tree-operations",
            metric_name="tip_count",
            parity_expectation="exact",
            r_value=extract_clade["r_tip_count"],
            bijux_value=extract_clade["bijux_tip_count"],
            tolerance_abs_diff=0.0,
        ),
        _row(
            row_id="vector-length",
            fragment_id="primate-longevity-vector-assembly",
            method_family="data-preparation",
            metric_name="vector_length",
            parity_expectation="exact",
            r_value=vector_assembly["r_vector_length"],
            bijux_value=vector_assembly["bijux_vector_length"],
            tolerance_abs_diff=0.0,
        ),
        _row(
            row_id="tree-join-tip-count",
            fragment_id="treeio-node-mapping-and-join",
            method_family="tree-operations",
            metric_name="joined_tip_count",
            parity_expectation="exact",
            r_value=tree_join["r_joined_tip_count"],
            bijux_value=tree_join["bijux_analysis_taxa_count"],
            tolerance_abs_diff=0.0,
        ),
    ]

    for node_name in ["Hylobates_lar", "Pan_paniscus", "Node32"]:
        rows.append(
            _row(
                row_id=f"tree-join-node-id-{node_name.lower()}",
                fragment_id="treeio-node-mapping-and-join",
                method_family="tree-operations",
                metric_name=f"node_id:{node_name}",
                parity_expectation="exact",
                r_value=tree_join["nodeid_examples_r"][node_name],
                bijux_value=tree_join["nodeid_examples_bijux"][node_name],
                tolerance_abs_diff=0.0,
            )
        )

    for dataset_name, dataset in sorted(random_lambda.items()):
        rows.append(
            _row(
                row_id=f"{dataset_name}-lambda",
                fragment_id="random-signal-lambda-fits",
                method_family="comparative-signal",
                metric_name="lambda_value",
                parity_expectation="statistical_tolerance",
                r_value=dataset["r"]["lambda_value"],
                bijux_value=dataset["bijux"]["lambda_value"],
                tolerance_abs_diff=0.001,
            )
        )
        rows.append(
            _row(
                row_id=f"{dataset_name}-log-likelihood",
                fragment_id="random-signal-lambda-fits",
                method_family="comparative-signal",
                metric_name="log_likelihood",
                parity_expectation="statistical_tolerance",
                r_value=dataset["r"]["log_likelihood"],
                bijux_value=dataset["bijux"]["log_likelihood"],
                tolerance_abs_diff=0.001,
            )
        )

    rows.extend(
        [
            _row(
                row_id="primate-lambda-value",
                fragment_id="primate-lambda-fit",
                method_family="comparative-signal",
                metric_name="lambda_value",
                parity_expectation="statistical_tolerance",
                r_value=primate_lambda["r"]["lambda_value"],
                bijux_value=primate_lambda["bijux"]["lambda_value"],
                tolerance_abs_diff=0.001,
            ),
            _row(
                row_id="primate-log-likelihood",
                fragment_id="primate-lambda-fit",
                method_family="comparative-signal",
                metric_name="log_likelihood",
                parity_expectation="statistical_tolerance",
                r_value=primate_lambda["r"]["log_likelihood"],
                bijux_value=primate_lambda["bijux"]["log_likelihood"],
                tolerance_abs_diff=0.001,
            ),
            _row(
                row_id="lambda-zero-likelihood-ratio-diff",
                fragment_id="lambda-zero-covariance-and-lrt",
                method_family="comparative-signal",
                metric_name="likelihood_ratio_abs_diff",
                parity_expectation="statistical_tolerance",
                r_value=None,
                bijux_value=None,
                tolerance_abs_diff=0.0001,
                observed_abs_diff=lambda_zero["likelihood_ratio_abs_diff"],
            ),
            _row(
                row_id="lambda-zero-p-value-diff",
                fragment_id="lambda-zero-covariance-and-lrt",
                method_family="comparative-signal",
                metric_name="p_value_abs_diff",
                parity_expectation="statistical_tolerance",
                r_value=None,
                bijux_value=None,
                tolerance_abs_diff=1e-12,
                observed_abs_diff=lambda_zero["p_value_abs_diff"],
            ),
            _row(
                row_id="lambda-zero-vcv-diff",
                fragment_id="lambda-zero-covariance-and-lrt",
                method_family="comparative-signal",
                metric_name="lambda0_vcv_top3_max_abs_diff",
                parity_expectation="statistical_tolerance",
                r_value=None,
                bijux_value=None,
                tolerance_abs_diff=1e-8,
                observed_abs_diff=lambda_zero["lambda0_vcv_top3_max_abs_diff"],
            ),
            _row(
                row_id="lambda-real-vcv-diff",
                fragment_id="lambda-zero-covariance-and-lrt",
                method_family="comparative-signal",
                metric_name="real_vcv_top3_max_abs_diff",
                parity_expectation="statistical_tolerance",
                r_value=None,
                bijux_value=None,
                tolerance_abs_diff=1e-8,
                observed_abs_diff=lambda_zero["real_vcv_top3_max_abs_diff"],
            ),
            _row(
                row_id="ancestral-point-max-diff",
                fragment_id="continuous-ancestral-point-estimates",
                method_family="ancestral-reconstruction",
                metric_name="point_max_abs_diff",
                parity_expectation="near_exact",
                r_value=None,
                bijux_value=None,
                tolerance_abs_diff=1e-9,
                observed_abs_diff=ancestral_points["point_max_abs_diff"],
            ),
            _row(
                row_id="ancestral-lower-interval-max-diff",
                fragment_id="continuous-ancestral-intervals",
                method_family="ancestral-reconstruction",
                metric_name="lower_95_max_abs_diff",
                parity_expectation="near_exact",
                r_value=None,
                bijux_value=None,
                tolerance_abs_diff=1e-9,
                observed_abs_diff=ancestral_intervals["lower_95_max_abs_diff"],
            ),
            _row(
                row_id="ancestral-upper-interval-max-diff",
                fragment_id="continuous-ancestral-intervals",
                method_family="ancestral-reconstruction",
                metric_name="upper_95_max_abs_diff",
                parity_expectation="near_exact",
                r_value=None,
                bijux_value=None,
                tolerance_abs_diff=1e-9,
                observed_abs_diff=ancestral_intervals["upper_95_max_abs_diff"],
            ),
            _row(
                row_id="mrca-estimate",
                fragment_id="bonobo-gibbon-mrca-estimate",
                method_family="ancestral-reconstruction",
                metric_name="mrca_estimate",
                parity_expectation="near_exact",
                r_value=mrca["r"]["mrca_estimate"][0]["estimate"],
                bijux_value=mrca["bijux"]["mrca_estimate"],
                tolerance_abs_diff=1e-9,
            ),
            _row(
                row_id="lifespan-increase-count",
                fragment_id="lifespan-increase-counts",
                method_family="ancestral-reconstruction",
                metric_name="increase_count",
                parity_expectation="exact",
                r_value=increase_counts["r"]["increase_count"],
                bijux_value=increase_counts["bijux"]["increase_count"],
                tolerance_abs_diff=0.0,
            ),
            _row(
                row_id="lifespan-increase-gt12-count",
                fragment_id="lifespan-increase-counts",
                method_family="ancestral-reconstruction",
                metric_name="increase_gt12_count",
                parity_expectation="exact",
                r_value=increase_counts["r"]["increase_gt12_count"],
                bijux_value=increase_counts["bijux"]["increase_gt12_count"],
                tolerance_abs_diff=0.0,
            ),
        ]
    )

    verdict_counts: dict[str, int] = {}
    for row in rows:
        verdict = row["verdict"]
        verdict_counts[verdict] = verdict_counts.get(verdict, 0) + 1

    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": EVIDENCE_ID,
        "row_count": len(rows),
        "verdict_counts": dict(sorted(verdict_counts.items())),
        "rows": rows,
    }


def _format_markdown_value(value: object | None) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def render_primate_scalar_parity_table_markdown(payload: dict[str, object]) -> str:
    lines = [
        "# Primate Scalar Parity Table",
        "",
        "| Row | Fragment | Metric | R | Bijux | Abs diff | Tolerance | Verdict |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload["rows"]:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["row_id"]),
                    str(row["fragment_id"]),
                    str(row["metric_name"]),
                    _format_markdown_value(row["r_value"]),
                    _format_markdown_value(row["bijux_value"]),
                    _format_markdown_value(row["observed_abs_diff"]),
                    _format_markdown_value(row["tolerance_abs_diff"]),
                    str(row["verdict"]),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            f"Rows: `{payload['row_count']}`",
            "",
        ]
    )
    return "\n".join(lines)
