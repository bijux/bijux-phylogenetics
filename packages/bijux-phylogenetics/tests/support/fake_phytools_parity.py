from __future__ import annotations

from pathlib import Path

from tests.support.fake_external_engines import write_executable


def fake_phytools_rscript(
    path: Path,
    *,
    phytools_available: bool = True,
    summary_overrides: dict[str, dict[str, object]] | None = None,
) -> Path:
    summary_payload = repr(summary_overrides or {})
    return write_executable(
        path,
        """#!/usr/bin/env python3
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

SUMMARIES = {
    "phylosig-lambda-non-ultrametric-strong-signal-twenty-four-taxa": {
        "taxon_count": 24,
        "trait_name": "signal_strong",
        "lambda_value": 1.0,
        "log_likelihood": -27.527955447046294,
        "null_log_likelihood": -38.15175747253939,
        "p_value": 1.0,
    },
    "phylosig-lambda-weak-signal-twenty-four-taxa": {
        "taxon_count": 24,
        "trait_name": "signal_weak",
        "lambda_value": 7.33736799610782e-05,
        "log_likelihood": -35.82077682773808,
        "null_log_likelihood": -35.82025224132399,
        "p_value": 1.0,
    },
    "phylosig-k-strong-signal-twenty-four-taxa": {
        "taxon_count": 24,
        "trait_name": "signal_strong",
        "k": 0.881409594863313,
        "p_value": 0.005025125628140704,
        "permutation_count": 199,
        "permutation_seed": 17,
        "null_distribution_count": 198,
        "simulated_k_minimum": 0.004227447597570447,
        "simulated_k_mean": 0.029614826253456215,
        "simulated_k_maximum": 0.1198277748473274,
    },
    "phylosig-k-weak-signal-twenty-four-taxa": {
        "taxon_count": 24,
        "trait_name": "signal_weak",
        "k": 0.03137456576480431,
        "p_value": 0.18090452261306519,
        "permutation_count": 199,
        "permutation_seed": 17,
        "null_distribution_count": 198,
        "simulated_k_minimum": 0.004692621984576097,
        "simulated_k_mean": 0.021817522911770332,
        "simulated_k_maximum": 0.10646570986048905,
    },
    "phylosig-lambda-primate-longevity": {
        "taxon_count": 90,
        "trait_name": "longevity",
        "lambda_value": 0.99957,
        "log_likelihood": -170.97952,
        "null_log_likelihood": -193.33681,
        "p_value": 0.0,
    },
}

SUMMARIES.update(__SUMMARY_OVERRIDES__)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_summary_table(path: Path, summary: dict[str, object]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\\t")
        writer.writeheader()
        for key, value in summary.items():
            writer.writerow({"metric": key, "value": value})


def write_rows_table(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def find_repo_root() -> Path:
    for candidate in (Path.cwd(), *Path.cwd().parents):
        package_root = candidate / "packages" / "bijux-phylogenetics"
        if package_root.is_dir():
            return package_root
    raise ValueError("could not locate bijux-phylogenetics package root")


def find_repo_python(package_root: Path) -> Path:
    candidates = (
        package_root / ".venv" / "bin" / "python",
        package_root.parent / ".venv" / "bin" / "python",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise ValueError("could not locate repo python for fake phytools parity payload")


def compute_continuous_ancestral_payload(
    case_payload: dict[str, object],
    *,
    estimator: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    package_root = find_repo_root()
    repo_python = find_repo_python(package_root)
    tree_path, traits_path = case_payload["input_fixtures"]
    inline_script = '''
import json
from pathlib import Path
from bijux_phylogenetics.ancestral.continuous import reconstruct_continuous_ancestral_states

payload = json.loads(Path(__PAYLOAD_PATH__).read_text(encoding="utf-8"))
tree_path, traits_path = payload["input_fixtures"]
report = reconstruct_continuous_ancestral_states(
    Path(tree_path),
    Path(traits_path),
    trait=payload["trait_name"],
    taxon_column=payload["taxon_column"],
    model="brownian",
    estimator=__ESTIMATOR__,
)
summary = {
    "taxon_count": report.taxon_count,
    "trait_name": report.trait,
    "internal_node_count": len([row for row in report.estimates if not row.is_tip]),
    "excluded_taxon_count": len(report.dropped_missing_taxa) + len(report.dropped_non_numeric_taxa),
    "excluded_taxa": sorted(report.dropped_missing_taxa + report.dropped_non_numeric_taxa),
    "tree_is_ultrametric": report.brownian_fit_diagnostics.tree_is_ultrametric,
}
if payload["operation"] == "continuous-ancestral-anc-ml":
    summary["sigma_squared"] = report.brownian_fit_diagnostics.residual_sigma_squared
    summary["log_likelihood"] = report.brownian_fit_diagnostics.log_likelihood
rows = [
    {
        "node": row.node,
        "estimate": row.estimate,
        "standard_error": row.standard_error,
        "lower_95_interval": row.lower_95_interval,
        "upper_95_interval": row.upper_95_interval,
    }
    for row in report.estimates
    if not row.is_tip
]
print(json.dumps({"summary": summary, "rows": rows}))
'''
    payload_path = package_root / "artifacts" / "fake-phytools-continuous-payload.json"
    payload_path.parent.mkdir(parents=True, exist_ok=True)
    payload_path.write_text(json.dumps(case_payload), encoding="utf-8")
    command = [
        str(repo_python),
        "-c",
        inline_script
        .replace("__PAYLOAD_PATH__", repr(str(payload_path)))
        .replace("__ESTIMATOR__", repr(estimator)),
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_root / "src")
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(package_root),
        env=env,
    )
    payload = json.loads(result.stdout)
    return payload["summary"], payload["rows"]


def compute_discrete_mk_payload(
    case_payload: dict[str, object],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    package_root = find_repo_root()
    repo_python = find_repo_python(package_root)
    inline_script = '''
import json
from pathlib import Path
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model

payload = json.loads(Path(__PAYLOAD_PATH__).read_text(encoding="utf-8"))
tree_path, traits_path = payload["input_fixtures"]
report = fit_discrete_mk_model(
    Path(tree_path),
    Path(traits_path),
    trait=payload["trait_name"],
    taxon_column=payload["taxon_column"],
    model=payload.get("discrete_model", "equal-rates"),
)
summary = {
    "taxon_count": report.taxon_count,
    "trait_name": report.trait,
    "excluded_taxon_count": len(report.input_audit.pruned_missing_value_taxa),
    "excluded_taxa": list(report.input_audit.pruned_missing_value_taxa),
    "model": report.model,
    "state_count": len(report.input_audit.observed_states),
    "parameter_count": report.parameter_count,
    "log_likelihood": report.log_likelihood,
    "aic": report.aic,
    "aicc": report.aicc,
    "overparameterized": report.overparameterized,
    "baseline_model": None if report.baseline_comparison is None else report.baseline_comparison.baseline_model,
    "preferred_model_by_aic": None if report.baseline_comparison is None else report.baseline_comparison.preferred_model_by_aic,
}
rows = [
    {
        "source_state": row.source_state,
        "target_state": row.target_state,
        "transition_allowed": row.transition_allowed,
        "step_distance": row.step_distance,
        "rate": row.rate,
    }
    for row in report.transition_rate_rows
]
print(json.dumps({"summary": summary, "rows": rows}))
'''
    payload_path = package_root / "artifacts" / "fake-phytools-discrete-payload.json"
    payload_path.parent.mkdir(parents=True, exist_ok=True)
    payload_path.write_text(json.dumps(case_payload), encoding="utf-8")
    command = [
        str(repo_python),
        "-c",
        inline_script.replace("__PAYLOAD_PATH__", repr(str(payload_path))),
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_root / "src")
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(package_root),
        env=env,
    )
    payload = json.loads(result.stdout)
    return payload["summary"], payload["rows"]


def compute_discrete_ancestral_payload(
    case_payload: dict[str, object],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    package_root = find_repo_root()
    repo_python = find_repo_python(package_root)
    inline_script = '''
import json
from pathlib import Path
from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states

payload = json.loads(Path(__PAYLOAD_PATH__).read_text(encoding="utf-8"))
tree_path, traits_path = payload["input_fixtures"]
report = reconstruct_discrete_ancestral_states(
    Path(tree_path),
    Path(traits_path),
    trait=payload["trait_name"],
    taxon_column=payload["taxon_column"],
    model=payload.get("discrete_model", "equal-rates"),
    root_prior_mode=payload.get("root_prior_mode", "equal"),
)
summary = {
    "taxon_count": report.taxon_count,
    "trait_name": report.trait,
    "excluded_taxon_count": len(report.dropped_missing_taxa),
    "excluded_taxa": list(report.dropped_missing_taxa),
    "model": report.model,
    "state_count": len(report.observed_states),
    "internal_node_count": len([row for row in report.estimates if not row.is_tip]),
    "root_prior_mode": report.root_prior_mode,
    "phytools_rerooting_method_comparable": report.rerooting_method_compatibility.comparable,
}
rows = [
    {
        "node": row.node,
        "state": state,
        "probability": probability,
    }
    for row in report.estimates
    if not row.is_tip
    for state, probability in row.state_probabilities.items()
]
print(json.dumps({"summary": summary, "rows": rows}))
'''
    payload_path = package_root / "artifacts" / "fake-phytools-rerooting-payload.json"
    payload_path.parent.mkdir(parents=True, exist_ok=True)
    payload_path.write_text(json.dumps(case_payload), encoding="utf-8")
    command = [
        str(repo_python),
        "-c",
        inline_script.replace("__PAYLOAD_PATH__", repr(str(payload_path))),
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_root / "src")
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(package_root),
        env=env,
    )
    payload = json.loads(result.stdout)
    return payload["summary"], payload["rows"]


def compute_stochastic_map_payload(
    case_payload: dict[str, object],
    *,
    include_branch_occupancy: bool,
    count_only: bool,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    package_root = find_repo_root()
    sys.path.insert(0, str(package_root / "src"))
    from bijux_phylogenetics.ancestral.common import load_discrete_dataset
    from bijux_phylogenetics.comparative.discrete_evolution import (
        count_discrete_stochastic_map_transitions,
        simulate_discrete_stochastic_maps,
        summarize_discrete_stochastic_map_density,
        summarize_discrete_stochastic_maps,
    )

    tree_path, traits_path = case_payload["input_fixtures"]
    dataset = load_discrete_dataset(
        Path(tree_path),
        Path(traits_path),
        trait=case_payload["trait_name"],
        taxon_column=case_payload["taxon_column"],
    )
    collection = simulate_discrete_stochastic_maps(
        Path(tree_path),
        Path(traits_path),
        trait=case_payload["trait_name"],
        taxon_column=case_payload["taxon_column"],
        model=case_payload.get("discrete_model", "equal-rates"),
        replicates=case_payload.get("stochastic_map_replicate_count", 128),
        seed=case_payload.get("stochastic_map_seed", 1),
    )
    summary_report = (
        collection.summary
        if case_payload["operation"] == "discrete-stochastic-map"
        else summarize_discrete_stochastic_maps(collection)
    )
    count_report = count_discrete_stochastic_map_transitions(collection)
    density_report = summarize_discrete_stochastic_map_density(
        collection,
        resolution=100
        if case_payload.get("density_resolution") is None
        else case_payload["density_resolution"],
        focal_state=case_payload.get("focal_state"),
    )
    summary = {
        "taxon_count": len(dataset.taxa),
        "trait_name": collection.trait,
        "excluded_taxon_count": len(dataset.dropped_missing_taxa),
        "excluded_taxa": list(dataset.dropped_missing_taxa),
        "model": collection.model,
        "state_count": len(collection.fit_audit.state_order),
        "parameter_count": collection.fit_audit.parameter_count,
        "log_likelihood": collection.fit_audit.log_likelihood,
        "aic": collection.fit_audit.aic,
        "aicc": collection.fit_audit.aicc,
        "overparameterized": collection.fit_audit.overparameterized,
        "baseline_model": collection.fit_audit.baseline_model,
        "preferred_model_by_aic": collection.fit_audit.preferred_model_by_aic,
        "requested_replicate_count": collection.replicates,
        "successful_replicate_count": (
            count_report.replicate_count
            if count_only
            else summary_report.replicate_count
        ),
        "simulation_failure_count": (
            len(collection.failures)
            if count_only
            else summary_report.simulation_failure_count
        ),
        "seed": collection.seed,
        "mean_total_transition_count": (
            count_report.mean_total_transition_count
            if count_only
            else summary_report.mean_total_transition_count
        ),
        "lower_95_total_transition_count": (
            count_report.lower_95_total_transition_count
            if count_only
            else summary_report.lower_95_total_transition_count
        ),
        "upper_95_total_transition_count": (
            count_report.upper_95_total_transition_count
            if count_only
            else summary_report.upper_95_total_transition_count
        ),
    }
    if case_payload["operation"] == "discrete-stochastic-map":
        summary["conditioned_on_node_estimates"] = (
            collection.conditioned_on_node_estimates
        )
    elif case_payload["operation"] == "discrete-stochastic-map-description":
        summary["branch_count"] = len(collection.maps[0].branch_histories)
    elif case_payload["operation"] == "discrete-stochastic-map-density":
        summary["branch_count"] = len(density_report.branch_rows)
        summary["focal_state"] = density_report.focal_state
        summary["baseline_state"] = density_report.baseline_state
        summary["resolution"] = density_report.resolution
        summary["total_tree_depth"] = density_report.total_tree_depth
    rows = (
        sorted(
            [
                {
                    "label": f"{row.parent_node}->{row.child_node}",
                    "mean_posterior_probability": row.mean_posterior_probability,
                    "minimum_posterior_probability": (
                        row.minimum_posterior_probability
                    ),
                    "maximum_posterior_probability": (
                        row.maximum_posterior_probability
                    ),
                    "uncertainty": row.uncertainty,
                    "slice_count": row.slice_count,
                }
                for row in density_report.branch_rows
            ],
            key=lambda row: row["label"],
        )
        if case_payload["operation"] == "discrete-stochastic-map-density"
        else sorted(
            [
                {
                    "row_kind": "transition_count",
                    "label": row.transition,
                    "mean_value": row.mean_count,
                    "lower_95_interval": row.lower_95_interval,
                    "upper_95_interval": row.upper_95_interval,
                    "presence_fraction": row.presence_fraction,
                }
                for row in (
                    count_report.aggregate_rows if count_only else summary_report.rows
                )
            ]
            + (
                []
                if count_only
                else [
                    {
                        "row_kind": "state_time",
                        "label": row.state,
                        "mean_value": row.mean_time,
                        "lower_95_interval": row.lower_95_interval,
                        "upper_95_interval": row.upper_95_interval,
                        "presence_fraction": 1.0,
                    }
                    for row in summary_report.state_time_rows
                ]
            ),
            key=lambda row: (row["row_kind"], row["label"]),
        )
    )
    if (
        include_branch_occupancy
        and not count_only
        and case_payload["operation"] != "discrete-stochastic-map-density"
    ):
        rows = sorted(
            rows
            + [
                {
                    "row_kind": "branch_state_occupancy",
                    "label": f"{row.parent_node}->{row.child_node}:{row.state}",
                    "mean_value": row.mean_time,
                    "lower_95_interval": row.lower_95_interval,
                    "upper_95_interval": row.upper_95_interval,
                    "presence_fraction": row.presence_fraction,
                }
                for row in summary_report.branch_occupancy_rows
            ],
            key=lambda row: (row["row_kind"], row["label"]),
        )
    return summary, rows


def compute_discrete_history_payload(
    case_payload: dict[str, object],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    package_root = find_repo_root()
    repo_python = find_repo_python(package_root)
    inline_script = '''
import json
from pathlib import Path
from bijux_phylogenetics.simulation import (
    DiscreteHistoryRateRow,
    simulate_discrete_histories,
)

payload = json.loads(Path(__PAYLOAD_PATH__).read_text(encoding="utf-8"))
tree_path = payload["input_fixtures"][0]
report = simulate_discrete_histories(
    Path(tree_path),
    states=list(payload["simulation_states"]),
    rate_rows=[
        DiscreteHistoryRateRow(
            source_state=row["source_state"],
            target_state=row["target_state"],
            rate=row["rate"],
        )
        for row in payload["simulation_rate_rows"]
    ],
    root_state=payload["simulation_root_state"],
    root_state_probabilities=payload["simulation_root_state_probabilities"],
    replicates=payload["simulation_replicate_count"],
    seed=payload["simulation_seed"],
)
summary = {
    "taxon_count": report.tip_count,
    "trait_name": payload["trait_name"],
    "branch_count": report.branch_count,
    "state_count": len(report.states),
    "requested_replicate_count": report.replicate_count,
    "successful_replicate_count": report.replicate_count,
    "fixed_root_state": report.fixed_root_state,
    "root_prior_probabilities": (
        None
        if report.fixed_root_state is not None
        else report.root_state_probabilities
    ),
    "seed": report.seed,
    "mean_total_transition_count": report.mean_total_transition_count,
    "lower_95_total_transition_count": report.lower_95_total_transition_count,
    "upper_95_total_transition_count": report.upper_95_total_transition_count,
}
rows = [
    {
        "row_kind": row.row_kind,
        "label": row.label,
        "mean_value": row.mean_value,
        "lower_95_interval": row.lower_95_interval,
        "upper_95_interval": row.upper_95_interval,
        "presence_fraction": row.presence_fraction,
    }
    for row in report.rows
]
print(json.dumps({"summary": summary, "rows": rows}))
'''
    payload_path = package_root / "artifacts" / "fake-phytools-sim-history-payload.json"
    payload_path.parent.mkdir(parents=True, exist_ok=True)
    payload_path.write_text(json.dumps(case_payload), encoding="utf-8")
    command = [
        str(repo_python),
        "-c",
        inline_script.replace("__PAYLOAD_PATH__", repr(str(payload_path))),
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_root / "src")
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(package_root),
        env=env,
    )
    payload = json.loads(result.stdout)
    return payload["summary"], payload["rows"]


def compute_continuous_brownian_payload(
    case_payload: dict[str, object],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    package_root = find_repo_root()
    repo_python = find_repo_python(package_root)
    inline_script = '''
import json
from pathlib import Path
from bijux_phylogenetics.simulation import simulate_brownian_trait_collection

payload = json.loads(Path(__PAYLOAD_PATH__).read_text(encoding="utf-8"))
tree_path = payload["input_fixtures"][0]
report = simulate_brownian_trait_collection(
    Path(tree_path),
    root_state=payload["continuous_root_state"],
    sigma_squared=payload["continuous_sigma_squared"],
    replicates=payload["continuous_replicate_count"],
    seed=payload["continuous_seed"],
)
summary = {
    "taxon_count": report.tip_count,
    "branch_count": report.branch_count,
    "requested_replicate_count": report.replicate_count,
    "successful_replicate_count": report.replicate_count,
    "seed": report.seed,
    "root_state": report.root_state,
    "sigma_squared": report.sigma_squared,
}
rows = [
    {
        "row_kind": row.row_kind,
        "label": row.label,
        "mean_value": "" if row.mean_value is None else row.mean_value,
        "standard_deviation": "" if row.standard_deviation is None else row.standard_deviation,
        "minimum": "" if row.minimum is None else row.minimum,
        "median": "" if row.median is None else row.median,
        "maximum": "" if row.maximum is None else row.maximum,
        "covariance": "" if row.covariance is None else row.covariance,
        "correlation": "" if row.correlation is None else row.correlation,
    }
    for row in report.rows
]
print(json.dumps({"summary": summary, "rows": rows}))
'''
    payload_path = package_root / "artifacts" / "fake-phytools-fastbm-payload.json"
    payload_path.parent.mkdir(parents=True, exist_ok=True)
    payload_path.write_text(json.dumps(case_payload), encoding="utf-8")
    command = [
        str(repo_python),
        "-c",
        inline_script.replace("__PAYLOAD_PATH__", repr(str(payload_path))),
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_root / "src")
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(package_root),
        env=env,
    )
    payload = json.loads(result.stdout)
    return payload["summary"], payload["rows"]


def compute_correlated_continuous_brownian_payload(
    case_payload: dict[str, object],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    package_root = find_repo_root()
    repo_python = find_repo_python(package_root)
    inline_script = '''
import json
from pathlib import Path
from bijux_phylogenetics.simulation import simulate_correlated_brownian_trait_collection

payload = json.loads(Path(__PAYLOAD_PATH__).read_text(encoding="utf-8"))
tree_path = payload["input_fixtures"][0]
report = simulate_correlated_brownian_trait_collection(
    Path(tree_path),
    trait_names=payload["continuous_trait_names"],
    evolutionary_covariance_matrix=payload["continuous_covariance_matrix"],
    root_states=payload["continuous_root_states"],
    replicates=payload["continuous_replicate_count"],
    seed=payload["continuous_seed"],
)
summary = {
    "taxon_count": report.tip_count,
    "branch_count": report.branch_count,
    "trait_count": len(report.trait_names),
    "requested_replicate_count": report.replicate_count,
    "successful_replicate_count": report.replicate_count,
    "seed": report.seed,
}
rows = [
    {
        "row_kind": row.row_kind,
        "label": row.label,
        "mean_value": "" if row.mean_value is None else row.mean_value,
        "standard_deviation": "" if row.standard_deviation is None else row.standard_deviation,
        "minimum": "" if row.minimum is None else row.minimum,
        "median": "" if row.median is None else row.median,
        "maximum": "" if row.maximum is None else row.maximum,
        "covariance": "" if row.covariance is None else row.covariance,
        "correlation": "" if row.correlation is None else row.correlation,
    }
    for row in report.rows
]
print(json.dumps({"summary": summary, "rows": rows}))
'''
    payload_path = package_root / "artifacts" / "fake-phytools-simcorrs-payload.json"
    payload_path.parent.mkdir(parents=True, exist_ok=True)
    payload_path.write_text(json.dumps(case_payload), encoding="utf-8")
    command = [
        str(repo_python),
        "-c",
        inline_script.replace("__PAYLOAD_PATH__", repr(str(payload_path))),
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_root / "src")
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(package_root),
        env=env,
    )
    payload = json.loads(result.stdout)
    return payload["summary"], payload["rows"]


def compute_pgls_payload(
    case_payload: dict[str, object],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    package_root = find_repo_root()
    repo_python = find_repo_python(package_root)
    inline_script = '''
import json
from pathlib import Path
from bijux_phylogenetics.comparative.pgls import (
    build_pgls_model_matrix,
    inspect_pgls_inputs,
    run_pgls,
)

payload = json.loads(Path(__PAYLOAD_PATH__).read_text(encoding="utf-8"))
tree_path, traits_path = payload["input_fixtures"]
input_report = inspect_pgls_inputs(
    Path(tree_path),
    Path(traits_path),
    formula=payload["comparative_formula"],
    taxon_column=payload["taxon_column"],
)
model_matrix = build_pgls_model_matrix(
    Path(tree_path),
    Path(traits_path),
    formula=payload["comparative_formula"],
    taxon_column=payload["taxon_column"],
)
report = run_pgls(
    Path(tree_path),
    Path(traits_path),
    formula=payload["comparative_formula"],
    taxon_column=payload["taxon_column"],
    lambda_value=payload.get("comparative_lambda_value", 1.0),
)
rows = []
for coefficient in report.coefficients:
    rows.extend([
        {"row_kind": "coefficient_estimate", "label": coefficient.name, "value": coefficient.estimate},
        {"row_kind": "coefficient_standard_error", "label": coefficient.name, "value": coefficient.standard_error},
        {"row_kind": "coefficient_p_value", "label": coefficient.name, "value": coefficient.p_value},
    ])
for matrix_row in model_matrix.rows:
    for column_name, value in matrix_row.encoded_values.items():
        rows.append(
            {
                "row_kind": "model_matrix",
                "label": f"{matrix_row.taxon}:{column_name}",
                "value": value,
            }
        )
rows = sorted(rows, key=lambda row: (row["row_kind"], row["label"]))
summary = {
    "taxon_count": report.taxon_count,
    "trait_name": report.response,
    "formula": report.formula.formula,
    "analysis_taxon_count": len(report.taxa),
    "coefficient_count": len(report.coefficients),
    "model_matrix_row_count": len(model_matrix.rows),
    "model_matrix_column_count": len(model_matrix.encoded_columns),
    "categorical_predictor_count": len(input_report.categorical_predictors),
    "interaction_term_count": len(report.interaction_terms),
    "lambda_value": report.lambda_value,
    "lambda_estimation_mode": report.lambda_fit.mode,
    "log_likelihood": report.log_likelihood,
    "aic": report.aic,
}
print(json.dumps({"summary": summary, "rows": rows}))
'''
    payload_path = package_root / "artifacts" / "fake-phytools-pgls-payload.json"
    payload_path.parent.mkdir(parents=True, exist_ok=True)
    payload_path.write_text(json.dumps(case_payload), encoding="utf-8")
    command = [
        str(repo_python),
        "-c",
        inline_script.replace("__PAYLOAD_PATH__", repr(str(payload_path))),
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_root / "src")
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(package_root),
        env=env,
    )
    payload = json.loads(result.stdout)
    return payload["summary"], payload["rows"]


def compute_phylogenetic_residual_payload(
    case_payload: dict[str, object],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    package_root = find_repo_root()
    repo_python = find_repo_python(package_root)
    inline_script = '''
import json
from pathlib import Path
from bijux_phylogenetics.comparative.regression import summarize_phylogenetic_residuals

payload = json.loads(Path(__PAYLOAD_PATH__).read_text(encoding="utf-8"))
tree_path, traits_path = payload["input_fixtures"]
method = (
    "brownian"
    if payload.get("comparative_lambda_value") == 1.0
    else "lambda"
)
report = summarize_phylogenetic_residuals(
    Path(tree_path),
    Path(traits_path),
    response=payload["trait_name"],
    predictor=payload["comparative_predictors"][0],
    taxon_column=payload["taxon_column"],
    method=method,
)
summary = {
    "taxon_count": report.analyzed_taxon_count,
    "trait_name": report.response,
    "predictor_name": report.predictor,
    "method": report.method,
    "excluded_taxon_count": len(report.excluded_taxa),
    "excluded_taxa": [row.taxon for row in report.excluded_taxa],
}
if method == "lambda":
    summary["lambda_value"] = report.lambda_value
    summary["log_likelihood"] = report.log_likelihood
rows = [
    {
        "row_kind": "coefficient_estimate",
        "label": row.name,
        "value": row.estimate,
    }
    for row in report.coefficient_rows
] + [
    {
        "row_kind": "taxon_value",
        "label": row.taxon,
        "observed_value": row.observed_value,
        "fitted_value": row.fitted_value,
        "residual": row.residual,
    }
    for row in report.taxon_rows
]
rows = sorted(rows, key=lambda row: (row["row_kind"], row["label"]))
print(json.dumps({"summary": summary, "rows": rows}))
'''
    payload_path = package_root / "artifacts" / "fake-phytools-phyl-resid-payload.json"
    payload_path.parent.mkdir(parents=True, exist_ok=True)
    payload_path.write_text(json.dumps(case_payload), encoding="utf-8")
    command = [
        str(repo_python),
        "-c",
        inline_script.replace("__PAYLOAD_PATH__", repr(str(payload_path))),
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_root / "src")
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(package_root),
        env=env,
    )
    payload = json.loads(result.stdout)
    return payload["summary"], payload["rows"]


def compute_phylogenetic_anova_payload(
    case_payload: dict[str, object],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    package_root = find_repo_root()
    repo_python = find_repo_python(package_root)
    inline_script = '''
import json
from pathlib import Path
from bijux_phylogenetics.comparative.regression import summarize_phylogenetic_anova

payload = json.loads(Path(__PAYLOAD_PATH__).read_text(encoding="utf-8"))
tree_path, traits_path = payload["input_fixtures"]
report = summarize_phylogenetic_anova(
    Path(tree_path),
    Path(traits_path),
    response=payload["trait_name"],
    group=payload["comparative_predictors"][0],
    taxon_column=payload["taxon_column"],
    simulations=payload.get("permutation_count", 199),
    seed=payload.get("permutation_seed", 1),
)
summary = {
    "taxon_count": report.analyzed_taxon_count,
    "trait_name": report.response,
    "group_column": report.group,
    "excluded_taxon_count": len(report.excluded_taxa),
    "excluded_taxa": [row.taxon for row in report.excluded_taxa],
    "group_count": report.group_count,
    "simulation_count": report.simulation_count,
    "seed": report.seed,
    "pairwise_adjustment_method": report.pairwise_adjustment_method,
    "brownian_sigma_squared": report.brownian_sigma_squared,
    "sum_of_squares_between": report.sum_of_squares_between,
    "sum_of_squares_within": report.sum_of_squares_within,
    "mean_square_between": report.mean_square_between,
    "mean_square_within": report.mean_square_within,
    "f_statistic": report.f_statistic,
    "p_value": report.p_value,
    "low_sample_group_count": report.low_sample_group_count,
}
rows = [
    {
        "row_kind": "group_summary",
        "label": row.group,
        "taxon_count": row.taxon_count,
        "taxa": ",".join(row.taxa),
        "mean": row.mean,
        "variance": row.variance,
        "minimum": row.minimum,
        "maximum": row.maximum,
    }
    for row in report.group_rows
] + [
    {
        "row_kind": "pairwise_comparison",
        "label": f"{row.left_group}|{row.right_group}",
        "left_taxon_count": row.left_taxon_count,
        "right_taxon_count": row.right_taxon_count,
        "observed_t_statistic": row.observed_t_statistic,
        "uncorrected_p_value": row.uncorrected_p_value,
        "adjusted_p_value": row.adjusted_p_value,
    }
    for row in report.pairwise_rows
]
rows = sorted(rows, key=lambda row: (row["row_kind"], row["label"]))
print(json.dumps({"summary": summary, "rows": rows}))
'''
    payload_path = package_root / "artifacts" / "fake-phytools-phyl-anova-payload.json"
    payload_path.parent.mkdir(parents=True, exist_ok=True)
    payload_path.write_text(json.dumps(case_payload), encoding="utf-8")
    command = [
        str(repo_python),
        "-c",
        inline_script.replace("__PAYLOAD_PATH__", repr(str(payload_path))),
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_root / "src")
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(package_root),
        env=env,
    )
    payload = json.loads(result.stdout)
    return payload["summary"], payload["rows"]


case_payload = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
output_root = Path(sys.argv[3])
output_root.mkdir(parents=True, exist_ok=True)
execution_path = output_root / "reference-execution.json"
summary_path = output_root / "reference-summary.json"
summary_table_path = output_root / "reference-summary.tsv"
fitmk_rows_path = output_root / "fitmk-rate-matrix.tsv"
stochastic_map_rows_path = output_root / "stochastic-map-summary-rows.tsv"
rerooting_rows_path = output_root / "rerooting-method-node-probabilities.tsv"
fast_anc_rows_path = output_root / "fast-anc-node-estimates.tsv"
anc_ml_rows_path = output_root / "anc-ml-node-estimates.tsv"
fastbm_rows_path = output_root / "fastbm-summary-rows.tsv"
simcorrs_rows_path = output_root / "simcorrs-summary-rows.tsv"
pgls_rows_path = output_root / "pgls-summary-rows.tsv"
phyl_resid_rows_path = output_root / "phyl-resid-summary-rows.tsv"
phyl_anova_rows_path = output_root / "phyl-anova-summary-rows.tsv"

if not __PHYTOOLS_AVAILABLE__:
    write_json(
        execution_path,
        {
            "status": "unavailable",
            "mismatch_reason": "phytools_package_unavailable",
            "r_version": "4.6.0",
            "phytools_version": None,
        },
    )
    raise SystemExit(0)

case_id = case_payload["case_id"]
if case_id not in SUMMARIES and case_payload["operation"] not in {
    "discrete-fit-mk",
    "discrete-stochastic-map",
    "discrete-stochastic-map-count",
    "discrete-stochastic-map-description",
    "discrete-stochastic-map-density",
    "simulate-discrete-history",
    "simulate-continuous-brownian",
    "simulate-continuous-correlated-brownian",
    "comparative-pgls-brownian",
    "phylogenetic-residuals",
    "phylogenetic-anova",
    "discrete-ancestral-rerooting",
    "continuous-ancestral-fast-anc",
    "continuous-ancestral-anc-ml",
}:
    write_json(
        execution_path,
        {
            "status": "failed",
            "mismatch_reason": "unsupported_operation",
            "error_type": "ValueError",
            "message": "unsupported phytools parity case: " + case_id,
            "r_version": "4.6.0",
            "phytools_version": "2.5.2",
        },
    )
    raise SystemExit(0)

summary = SUMMARIES.get(case_id)
rows = None
if case_payload["operation"] == "discrete-fit-mk":
    summary, rows = compute_discrete_mk_payload(case_payload)
elif case_payload["operation"] == "discrete-stochastic-map":
    summary, rows = compute_stochastic_map_payload(
        case_payload,
        include_branch_occupancy=False,
        count_only=False,
    )
elif case_payload["operation"] == "discrete-stochastic-map-count":
    summary, rows = compute_stochastic_map_payload(
        case_payload,
        include_branch_occupancy=False,
        count_only=True,
    )
elif case_payload["operation"] == "discrete-stochastic-map-description":
    summary, rows = compute_stochastic_map_payload(
        case_payload,
        include_branch_occupancy=True,
        count_only=False,
    )
elif case_payload["operation"] == "discrete-stochastic-map-density":
    summary, rows = compute_stochastic_map_payload(
        case_payload,
        include_branch_occupancy=False,
        count_only=False,
    )
elif case_payload["operation"] == "simulate-discrete-history":
    summary, rows = compute_discrete_history_payload(case_payload)
elif case_payload["operation"] == "simulate-continuous-brownian":
    summary, rows = compute_continuous_brownian_payload(case_payload)
elif case_payload["operation"] == "simulate-continuous-correlated-brownian":
    summary, rows = compute_correlated_continuous_brownian_payload(case_payload)
elif case_payload["operation"] == "comparative-pgls-brownian":
    summary, rows = compute_pgls_payload(case_payload)
elif case_payload["operation"] == "phylogenetic-residuals":
    summary, rows = compute_phylogenetic_residual_payload(case_payload)
elif case_payload["operation"] == "phylogenetic-anova":
    summary, rows = compute_phylogenetic_anova_payload(case_payload)
elif case_payload["operation"] == "discrete-ancestral-rerooting":
    summary, rows = compute_discrete_ancestral_payload(case_payload)
elif case_payload["operation"] == "continuous-ancestral-fast-anc":
    summary, rows = compute_continuous_ancestral_payload(
        case_payload,
        estimator="fast-anc",
    )
elif case_payload["operation"] == "continuous-ancestral-anc-ml":
    summary, rows = compute_continuous_ancestral_payload(
        case_payload,
        estimator="anc-ml",
    )
write_json(
    execution_path,
    {
        "status": "ok",
        "r_version": "4.6.0",
        "phytools_version": "2.5.2",
    },
)
write_json(summary_path, summary)
write_summary_table(summary_table_path, summary)
if rows is not None:
    if case_payload["operation"] == "discrete-fit-mk":
        write_rows_table(fitmk_rows_path, rows)
    elif case_payload["operation"] in {
        "discrete-stochastic-map",
        "discrete-stochastic-map-count",
        "discrete-stochastic-map-description",
        "discrete-stochastic-map-density",
        "simulate-discrete-history",
    }:
        write_rows_table(stochastic_map_rows_path, rows)
    elif case_payload["operation"] == "simulate-continuous-brownian":
        write_rows_table(fastbm_rows_path, rows)
    elif case_payload["operation"] == "simulate-continuous-correlated-brownian":
        write_rows_table(simcorrs_rows_path, rows)
    elif case_payload["operation"] == "comparative-pgls-brownian":
        write_rows_table(pgls_rows_path, rows)
    elif case_payload["operation"] == "phylogenetic-residuals":
        write_rows_table(phyl_resid_rows_path, rows)
    elif case_payload["operation"] == "phylogenetic-anova":
        write_rows_table(phyl_anova_rows_path, rows)
    elif case_payload["operation"] == "discrete-ancestral-rerooting":
        write_rows_table(rerooting_rows_path, rows)
    elif case_payload["operation"] == "continuous-ancestral-fast-anc":
        write_rows_table(fast_anc_rows_path, rows)
    elif case_payload["operation"] == "continuous-ancestral-anc-ml":
        write_rows_table(anc_ml_rows_path, rows)
""".replace("__SUMMARY_OVERRIDES__", summary_payload).replace(
            "__PHYTOOLS_AVAILABLE__", repr(phytools_available)
        ),
    )
