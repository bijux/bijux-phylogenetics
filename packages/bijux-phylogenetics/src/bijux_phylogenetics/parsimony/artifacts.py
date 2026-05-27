from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.reports.service.artifacts import write_json_artifact

from .models import (
    CaminSokalScoreReport,
    DolloScoreReport,
    FitchScoreReport,
    ParsimonyConsistencyIndexReport,
    ParsimonyReconstructionReport,
    ParsimonyTreeLengthReport,
    SankoffScoreReport,
    WagnerScoreReport,
)


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


def write_sankoff_steps_table(path: Path, report: SankoffScoreReport) -> Path:
    """Write one deterministic per-character Sankoff minimum-cost table."""
    return write_taxon_rows(
        path,
        columns=["character_id", "minimum_cost", "observed_states", "matrix_states"],
        rows=[
            {
                "character_id": row.character_id,
                "minimum_cost": row.minimum_cost,
                "observed_states": "|".join(row.observed_states),
                "matrix_states": "|".join(row.matrix_states),
            }
            for row in report.step_rows
        ],
    )


def write_sankoff_node_cost_table(path: Path, report: SankoffScoreReport) -> Path:
    """Write one deterministic per-node Sankoff cost-vector table."""
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


def write_sankoff_node_selection_table(path: Path, report: SankoffScoreReport) -> Path:
    """Write one deterministic per-node Sankoff optimal-state table."""
    return write_taxon_rows(
        path,
        columns=[
            "character_id",
            "node",
            "node_name",
            "descendant_taxa",
            "optimal_states",
            "tie_states",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "node": row.node,
                "node_name": row.node_name,
                "descendant_taxa": "|".join(row.descendant_taxa),
                "optimal_states": "|".join(row.optimal_states),
                "tie_states": "|".join(row.tie_states),
            }
            for row in report.node_selection_rows
        ],
    )


def write_sankoff_run_json(path: Path, report: SankoffScoreReport) -> Path:
    """Write one machine-readable Sankoff run payload."""
    return write_json_artifact(
        path,
        {
            "algorithm": report.algorithm,
            "tree_path": None if report.tree_path is None else str(report.tree_path),
            "matrix_path": None
            if report.matrix_path is None
            else str(report.matrix_path),
            "cost_matrix_path": None
            if report.cost_matrix_path is None
            else str(report.cost_matrix_path),
            "taxon_column": report.taxon_column,
            "taxon_count": report.taxon_count,
            "character_count": report.character_count,
            "total_cost": report.total_cost,
            "step_rows": [
                {
                    "character_id": row.character_id,
                    "minimum_cost": row.minimum_cost,
                    "observed_states": row.observed_states,
                    "matrix_states": row.matrix_states,
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
            "node_selection_rows": [
                {
                    "character_id": row.character_id,
                    "node": row.node,
                    "node_name": row.node_name,
                    "descendant_taxa": row.descendant_taxa,
                    "optimal_states": row.optimal_states,
                    "tie_states": row.tie_states,
                }
                for row in report.node_selection_rows
            ],
        },
    )


def write_sankoff_artifacts(
    out_dir: Path,
    report: SankoffScoreReport,
) -> dict[str, Path]:
    """Write the governed Sankoff artifact family."""
    out_dir.mkdir(parents=True, exist_ok=True)
    steps_path = write_sankoff_steps_table(out_dir / "steps.tsv", report)
    node_costs_path = write_sankoff_node_cost_table(out_dir / "node_costs.tsv", report)
    node_selection_path = write_sankoff_node_selection_table(
        out_dir / "selected_states.tsv",
        report,
    )
    run_json_path = write_sankoff_run_json(out_dir / "run.json", report)
    return {
        "steps_path": steps_path,
        "node_costs_path": node_costs_path,
        "node_selection_path": node_selection_path,
        "run_json_path": run_json_path,
    }


def write_dollo_steps_table(path: Path, report: DolloScoreReport) -> Path:
    """Write one deterministic per-character Dollo summary table."""
    return write_taxon_rows(
        path,
        columns=[
            "character_id",
            "derived_taxon_count",
            "gain_node",
            "gain_node_name",
            "gain_descendant_taxa",
            "total_losses",
            "impossible_state_warning",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "derived_taxon_count": row.derived_taxon_count,
                "gain_node": row.gain_node,
                "gain_node_name": row.gain_node_name,
                "gain_descendant_taxa": "|".join(row.gain_descendant_taxa),
                "total_losses": row.total_losses,
                "impossible_state_warning": row.impossible_state_warning,
            }
            for row in report.step_rows
        ],
    )


def write_dollo_branch_change_table(path: Path, report: DolloScoreReport) -> Path:
    """Write one deterministic Dollo gain/loss branch table."""
    return write_taxon_rows(
        path,
        columns=["character_id", "change_kind", "node", "node_name", "descendant_taxa"],
        rows=[
            {
                "character_id": row.character_id,
                "change_kind": row.change_kind,
                "node": row.node,
                "node_name": row.node_name,
                "descendant_taxa": "|".join(row.descendant_taxa),
            }
            for row in report.branch_change_rows
        ],
    )


def write_dollo_run_json(path: Path, report: DolloScoreReport) -> Path:
    """Write one machine-readable Dollo run payload."""
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
            "total_gains": report.total_gains,
            "total_losses": report.total_losses,
            "step_rows": [
                {
                    "character_id": row.character_id,
                    "derived_taxon_count": row.derived_taxon_count,
                    "gain_node": row.gain_node,
                    "gain_node_name": row.gain_node_name,
                    "gain_descendant_taxa": row.gain_descendant_taxa,
                    "total_losses": row.total_losses,
                    "impossible_state_warning": row.impossible_state_warning,
                }
                for row in report.step_rows
            ],
            "branch_change_rows": [
                {
                    "character_id": row.character_id,
                    "change_kind": row.change_kind,
                    "node": row.node,
                    "node_name": row.node_name,
                    "descendant_taxa": row.descendant_taxa,
                }
                for row in report.branch_change_rows
            ],
        },
    )


def write_dollo_artifacts(
    out_dir: Path,
    report: DolloScoreReport,
) -> dict[str, Path]:
    """Write the governed Dollo artifact family."""
    out_dir.mkdir(parents=True, exist_ok=True)
    steps_path = write_dollo_steps_table(out_dir / "steps.tsv", report)
    branch_changes_path = write_dollo_branch_change_table(
        out_dir / "branch_changes.tsv",
        report,
    )
    run_json_path = write_dollo_run_json(out_dir / "run.json", report)
    return {
        "steps_path": steps_path,
        "branch_changes_path": branch_changes_path,
        "run_json_path": run_json_path,
    }


def write_camin_sokal_steps_table(path: Path, report: CaminSokalScoreReport) -> Path:
    """Write one deterministic per-character Camin-Sokal summary table."""
    return write_taxon_rows(
        path,
        columns=["character_id", "derived_taxon_count", "gain_count", "root_state"],
        rows=[
            {
                "character_id": row.character_id,
                "derived_taxon_count": row.derived_taxon_count,
                "gain_count": row.gain_count,
                "root_state": row.root_state,
            }
            for row in report.step_rows
        ],
    )


def write_camin_sokal_branch_change_table(
    path: Path,
    report: CaminSokalScoreReport,
) -> Path:
    """Write one deterministic irreversible Camin-Sokal branch-change table."""
    return write_taxon_rows(
        path,
        columns=["character_id", "change_kind", "node", "node_name", "descendant_taxa"],
        rows=[
            {
                "character_id": row.character_id,
                "change_kind": row.change_kind,
                "node": row.node,
                "node_name": row.node_name,
                "descendant_taxa": "|".join(row.descendant_taxa),
            }
            for row in report.branch_change_rows
        ],
    )


def write_camin_sokal_run_json(path: Path, report: CaminSokalScoreReport) -> Path:
    """Write one machine-readable Camin-Sokal run payload."""
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
            "root_state": report.root_state,
            "total_gains": report.total_gains,
            "step_rows": [
                {
                    "character_id": row.character_id,
                    "derived_taxon_count": row.derived_taxon_count,
                    "gain_count": row.gain_count,
                    "root_state": row.root_state,
                }
                for row in report.step_rows
            ],
            "branch_change_rows": [
                {
                    "character_id": row.character_id,
                    "change_kind": row.change_kind,
                    "node": row.node,
                    "node_name": row.node_name,
                    "descendant_taxa": row.descendant_taxa,
                }
                for row in report.branch_change_rows
            ],
        },
    )


def write_camin_sokal_artifacts(
    out_dir: Path,
    report: CaminSokalScoreReport,
) -> dict[str, Path]:
    """Write the governed Camin-Sokal artifact family."""
    out_dir.mkdir(parents=True, exist_ok=True)
    steps_path = write_camin_sokal_steps_table(out_dir / "steps.tsv", report)
    branch_changes_path = write_camin_sokal_branch_change_table(
        out_dir / "branch_changes.tsv",
        report,
    )
    run_json_path = write_camin_sokal_run_json(out_dir / "run.json", report)
    return {
        "steps_path": steps_path,
        "branch_changes_path": branch_changes_path,
        "run_json_path": run_json_path,
    }


def write_parsimony_reconstruction_steps_table(
    path: Path,
    report: ParsimonyReconstructionReport,
) -> Path:
    """Write one deterministic per-character parsimony reconstruction summary table."""
    return write_taxon_rows(
        path,
        columns=["character_id", "step_count", "observed_states", "root_state"],
        rows=[
            {
                "character_id": row.character_id,
                "step_count": row.step_count,
                "observed_states": "|".join(row.observed_states),
                "root_state": row.root_state,
            }
            for row in report.step_rows
        ],
    )


def write_parsimony_reconstruction_node_state_table(
    path: Path,
    report: ParsimonyReconstructionReport,
) -> Path:
    """Write one deterministic resolved-node-state table for a parsimony reconstruction."""
    return write_taxon_rows(
        path,
        columns=[
            "character_id",
            "node",
            "node_name",
            "descendant_taxa",
            "resolved_state",
            "is_tip",
            "observed_state",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "node": row.node,
                "node_name": row.node_name,
                "descendant_taxa": "|".join(row.descendant_taxa),
                "resolved_state": row.resolved_state,
                "is_tip": row.is_tip,
                "observed_state": row.observed_state,
            }
            for row in report.node_state_rows
        ],
    )


def write_parsimony_reconstruction_branch_change_table(
    path: Path,
    report: ParsimonyReconstructionReport,
) -> Path:
    """Write one deterministic branch-change table for a parsimony reconstruction."""
    return write_taxon_rows(
        path,
        columns=[
            "character_id",
            "parent_node",
            "parent_state",
            "node",
            "node_name",
            "descendant_taxa",
            "change_from",
            "change_to",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "parent_node": row.parent_node,
                "parent_state": row.parent_state,
                "node": row.node,
                "node_name": row.node_name,
                "descendant_taxa": "|".join(row.descendant_taxa),
                "change_from": row.change_from,
                "change_to": row.change_to,
            }
            for row in report.branch_change_rows
        ],
    )


def write_parsimony_reconstruction_run_json(
    path: Path,
    report: ParsimonyReconstructionReport,
) -> Path:
    """Write one machine-readable parsimony reconstruction payload."""
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
                    "root_state": row.root_state,
                }
                for row in report.step_rows
            ],
            "node_state_rows": [
                {
                    "character_id": row.character_id,
                    "node": row.node,
                    "node_name": row.node_name,
                    "descendant_taxa": row.descendant_taxa,
                    "resolved_state": row.resolved_state,
                    "is_tip": row.is_tip,
                    "observed_state": row.observed_state,
                }
                for row in report.node_state_rows
            ],
            "branch_change_rows": [
                {
                    "character_id": row.character_id,
                    "parent_node": row.parent_node,
                    "parent_state": row.parent_state,
                    "node": row.node,
                    "node_name": row.node_name,
                    "descendant_taxa": row.descendant_taxa,
                    "change_from": row.change_from,
                    "change_to": row.change_to,
                }
                for row in report.branch_change_rows
            ],
        },
    )


def write_parsimony_reconstruction_artifacts(
    out_dir: Path,
    report: ParsimonyReconstructionReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one parsimony reconstruction."""
    out_dir.mkdir(parents=True, exist_ok=True)
    steps_path = write_parsimony_reconstruction_steps_table(out_dir / "steps.tsv", report)
    node_states_path = write_parsimony_reconstruction_node_state_table(
        out_dir / "resolved_states.tsv",
        report,
    )
    branch_changes_path = write_parsimony_reconstruction_branch_change_table(
        out_dir / "branch_changes.tsv",
        report,
    )
    run_json_path = write_parsimony_reconstruction_run_json(out_dir / "run.json", report)
    return {
        "steps_path": steps_path,
        "node_states_path": node_states_path,
        "branch_changes_path": branch_changes_path,
        "run_json_path": run_json_path,
    }


def write_parsimony_tree_length_scores_table(
    path: Path,
    report: ParsimonyTreeLengthReport,
) -> Path:
    """Write one deterministic per-character tree-length table."""
    return write_taxon_rows(
        path,
        columns=["character_id", "raw_score", "character_weight", "weighted_score"],
        rows=[
            {
                "character_id": row.character_id,
                "raw_score": row.raw_score,
                "character_weight": row.character_weight,
                "weighted_score": row.weighted_score,
            }
            for row in report.step_rows
        ],
    )


def write_parsimony_tree_length_run_json(
    path: Path,
    report: ParsimonyTreeLengthReport,
) -> Path:
    """Write one machine-readable tree-length payload."""
    return write_json_artifact(
        path,
        {
            "algorithm": report.algorithm,
            "method": report.method,
            "tree_path": None if report.tree_path is None else str(report.tree_path),
            "matrix_path": None
            if report.matrix_path is None
            else str(report.matrix_path),
            "cost_matrix_path": None
            if report.cost_matrix_path is None
            else str(report.cost_matrix_path),
            "weights_path": None
            if report.weights_path is None
            else str(report.weights_path),
            "taxon_column": report.taxon_column,
            "taxon_count": report.taxon_count,
            "character_count": report.character_count,
            "raw_total_score": report.raw_total_score,
            "total_score": report.total_score,
            "step_rows": [
                {
                    "character_id": row.character_id,
                    "raw_score": row.raw_score,
                    "character_weight": row.character_weight,
                    "weighted_score": row.weighted_score,
                }
                for row in report.step_rows
            ],
        },
    )


def write_parsimony_tree_length_artifacts(
    out_dir: Path,
    report: ParsimonyTreeLengthReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one tree-length run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    scores_path = write_parsimony_tree_length_scores_table(
        out_dir / "character_scores.tsv",
        report,
    )
    run_json_path = write_parsimony_tree_length_run_json(out_dir / "run.json", report)
    return {
        "scores_path": scores_path,
        "run_json_path": run_json_path,
    }


def write_parsimony_consistency_index_table(
    path: Path,
    report: ParsimonyConsistencyIndexReport,
) -> Path:
    """Write one deterministic per-character consistency-index table."""
    return write_taxon_rows(
        path,
        columns=[
            "character_id",
            "character_kind",
            "observed_states",
            "minimum_possible_steps",
            "observed_steps",
            "consistency_index",
            "undefined_reason",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "character_kind": row.character_kind,
                "observed_states": "|".join(row.observed_states),
                "minimum_possible_steps": row.minimum_possible_steps,
                "observed_steps": row.observed_steps,
                "consistency_index": row.consistency_index,
                "undefined_reason": row.undefined_reason,
            }
            for row in report.character_rows
        ],
    )


def write_parsimony_consistency_run_json(
    path: Path,
    report: ParsimonyConsistencyIndexReport,
) -> Path:
    """Write one machine-readable consistency-index payload."""
    return write_json_artifact(
        path,
        {
            "algorithm": report.algorithm,
            "method": report.method,
            "tree_path": None if report.tree_path is None else str(report.tree_path),
            "matrix_path": None
            if report.matrix_path is None
            else str(report.matrix_path),
            "taxon_column": report.taxon_column,
            "taxon_count": report.taxon_count,
            "character_count": report.character_count,
            "included_character_count": report.included_character_count,
            "excluded_character_count": report.excluded_character_count,
            "minimum_possible_steps_total": report.minimum_possible_steps_total,
            "observed_steps_total": report.observed_steps_total,
            "consistency_index": report.consistency_index,
            "undefined_reason": report.undefined_reason,
            "character_rows": [
                {
                    "character_id": row.character_id,
                    "character_kind": row.character_kind,
                    "observed_states": row.observed_states,
                    "minimum_possible_steps": row.minimum_possible_steps,
                    "observed_steps": row.observed_steps,
                    "consistency_index": row.consistency_index,
                    "undefined_reason": row.undefined_reason,
                }
                for row in report.character_rows
            ],
        },
    )


def write_parsimony_consistency_artifacts(
    out_dir: Path,
    report: ParsimonyConsistencyIndexReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one consistency-index run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    indices_path = write_parsimony_consistency_index_table(
        out_dir / "character_indices.tsv",
        report,
    )
    run_json_path = write_parsimony_consistency_run_json(out_dir / "run.json", report)
    return {
        "indices_path": indices_path,
        "run_json_path": run_json_path,
    }
