from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.reports.service.artifacts import write_json_artifact

from .models import FitchScoreReport, WagnerScoreReport


def write_fitch_steps_table(path: Path, report: FitchScoreReport) -> Path:
    """Write one deterministic per-character Fitch tree-length table."""
    return write_taxon_rows(
        path,
        columns=["character_id", "step_count", "observed_states"],
        rows=[
            {
                "character_id": row.character_id,
                "step_count": row.step_count,
                "observed_states": "|".join(row.observed_states),
            }
            for row in report.step_rows
        ],
    )


def write_fitch_node_state_set_table(path: Path, report: FitchScoreReport) -> Path:
    """Write one deterministic per-node Fitch candidate-state table."""
    return write_taxon_rows(
        path,
        columns=[
            "character_id",
            "node",
            "node_name",
            "descendant_taxa",
            "is_tip",
            "observed_state",
            "state_set",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "node": row.node,
                "node_name": row.node_name,
                "descendant_taxa": "|".join(row.descendant_taxa),
                "is_tip": row.is_tip,
                "observed_state": row.observed_state,
                "state_set": "|".join(row.state_set),
            }
            for row in report.node_state_rows
        ],
    )


def write_fitch_run_json(path: Path, report: FitchScoreReport) -> Path:
    """Write one machine-readable unordered Fitch run payload."""
    return write_json_artifact(
        path,
        {
            "algorithm": report.algorithm,
            "tree_path": None if report.tree_path is None else str(report.tree_path),
            "matrix_path": None
            if report.matrix_path is None
            else str(report.matrix_path),
            "taxon_column": report.taxon_column,
            "taxon_count": report.taxon_count,
            "character_count": report.character_count,
            "total_steps": report.total_steps,
            "step_rows": [
                {
                    "character_id": row.character_id,
                    "step_count": row.step_count,
                    "observed_states": row.observed_states,
                }
                for row in report.step_rows
            ],
            "node_state_rows": [
                {
                    "character_id": row.character_id,
                    "node": row.node,
                    "node_name": row.node_name,
                    "descendant_taxa": row.descendant_taxa,
                    "state_set": row.state_set,
                    "is_tip": row.is_tip,
                    "observed_state": row.observed_state,
                }
                for row in report.node_state_rows
            ],
        },
    )


def write_fitch_artifacts(
    out_dir: Path,
    report: FitchScoreReport,
) -> dict[str, Path]:
    """Write the governed unordered Fitch artifact family."""
    out_dir.mkdir(parents=True, exist_ok=True)
    steps_path = write_fitch_steps_table(out_dir / "steps.tsv", report)
    node_state_sets_path = write_fitch_node_state_set_table(
        out_dir / "node_state_sets.tsv",
        report,
    )
    run_json_path = write_fitch_run_json(out_dir / "run.json", report)
    return {
        "steps_path": steps_path,
        "node_state_sets_path": node_state_sets_path,
        "run_json_path": run_json_path,
    }


def write_wagner_steps_table(path: Path, report: WagnerScoreReport) -> Path:
    """Write one deterministic per-character Wagner weighted-step table."""
    return write_taxon_rows(
        path,
        columns=[
            "character_id",
            "weighted_step_count",
            "observed_states",
            "state_order",
            "optimal_root_states",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "weighted_step_count": row.weighted_step_count,
                "observed_states": "|".join(row.observed_states),
                "state_order": "|".join(row.state_order),
                "optimal_root_states": "|".join(row.optimal_root_states),
            }
            for row in report.step_rows
        ],
    )


def write_wagner_node_cost_table(path: Path, report: WagnerScoreReport) -> Path:
    """Write one deterministic per-node ordered Wagner cost table."""
    return write_taxon_rows(
        path,
        columns=[
            "character_id",
            "node",
            "node_name",
            "descendant_taxa",
            "state",
            "cost",
            "is_optimal_state",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "node": row.node,
                "node_name": row.node_name,
                "descendant_taxa": "|".join(row.descendant_taxa),
                "state": row.state,
                "cost": row.cost,
                "is_optimal_state": row.is_optimal_state,
            }
            for row in report.node_cost_rows
        ],
    )


def write_wagner_run_json(path: Path, report: WagnerScoreReport) -> Path:
    """Write one machine-readable ordered Wagner run payload."""
    return write_json_artifact(
        path,
        {
            "algorithm": report.algorithm,
            "tree_path": None if report.tree_path is None else str(report.tree_path),
            "matrix_path": None
            if report.matrix_path is None
            else str(report.matrix_path),
            "taxon_column": report.taxon_column,
            "taxon_count": report.taxon_count,
            "character_count": report.character_count,
            "total_cost": report.total_cost,
            "step_rows": [
                {
                    "character_id": row.character_id,
                    "weighted_step_count": row.weighted_step_count,
                    "observed_states": row.observed_states,
                    "state_order": row.state_order,
                    "optimal_root_states": row.optimal_root_states,
                }
                for row in report.step_rows
            ],
            "node_cost_rows": [
                {
                    "character_id": row.character_id,
                    "node": row.node,
                    "node_name": row.node_name,
                    "descendant_taxa": row.descendant_taxa,
                    "state": row.state,
                    "cost": row.cost,
                    "is_optimal_state": row.is_optimal_state,
                }
                for row in report.node_cost_rows
            ],
        },
    )


def write_wagner_artifacts(
    out_dir: Path,
    report: WagnerScoreReport,
) -> dict[str, Path]:
    """Write the governed ordered Wagner artifact family."""
    out_dir.mkdir(parents=True, exist_ok=True)
    steps_path = write_wagner_steps_table(out_dir / "steps.tsv", report)
    node_costs_path = write_wagner_node_cost_table(out_dir / "node_costs.tsv", report)
    run_json_path = write_wagner_run_json(out_dir / "run.json", report)
    return {
        "steps_path": steps_path,
        "node_costs_path": node_costs_path,
        "run_json_path": run_json_path,
    }
