from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.newick import (
    loads_newick,
    write_newick,
    write_newick_tree_set,
)
from bijux_phylogenetics.reports.service.artifacts import write_json_artifact
from bijux_phylogenetics.trees import (
    compute_clade_frequency_table,
    compute_consensus_tree,
    write_clade_frequency_table,
)

from .models import (
    CaminSokalScoreReport,
    DolloScoreReport,
    FitchScoreReport,
    ParsimonyBootstrapReport,
    ParsimonyBremerSupportReport,
    ParsimonyConsistencyIndexReport,
    ParsimonyEqualBestConsensusReport,
    ParsimonyJackknifeReport,
    ParsimonyNniSearchReport,
    ParsimonyPlacementReport,
    ParsimonyRatchetReport,
    ParsimonyReconstructionReport,
    ParsimonyRescaledConsistencyIndexReport,
    ParsimonyRetentionIndexReport,
    ParsimonySprSearchReport,
    ParsimonyTreeLengthReport,
    SankoffScoreReport,
    WagnerScoreReport,
)


def write_fitch_steps_table(path: Path, report: FitchScoreReport) -> Path:
    """Write one deterministic per-character Fitch tree-length table."""
    return write_taxon_rows(
        path,
        columns=[
            "character_id",
            "step_count",
            "observed_states",
            "character_weight",
            "weighted_score",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "step_count": row.step_count,
                "observed_states": "|".join(row.observed_states),
                "character_weight": row.character_weight,
                "weighted_score": row.weighted_score,
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
            "weights_path": None
            if report.weights_path is None
            else str(report.weights_path),
            "total_weighted_score": report.total_weighted_score,
            "step_rows": [
                {
                    "character_id": row.character_id,
                    "step_count": row.step_count,
                    "observed_states": row.observed_states,
                    "character_weight": row.character_weight,
                    "weighted_score": row.weighted_score,
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
            "character_weight",
            "weighted_score",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "weighted_step_count": row.weighted_step_count,
                "observed_states": "|".join(row.observed_states),
                "state_order": "|".join(row.state_order),
                "optimal_root_states": "|".join(row.optimal_root_states),
                "character_weight": row.character_weight,
                "weighted_score": row.weighted_score,
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
    validation_warnings = getattr(report, "validation_warnings", [])
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
            "weights_path": None
            if report.weights_path is None
            else str(report.weights_path),
            "total_weighted_score": report.total_weighted_score,
            "validation_warnings": [
                {
                    "code": warning.code,
                    "message": warning.message,
                    "details": warning.details,
                }
                for warning in validation_warnings
            ],
            "step_rows": [
                {
                    "character_id": row.character_id,
                    "weighted_step_count": row.weighted_step_count,
                    "observed_states": row.observed_states,
                    "state_order": row.state_order,
                    "optimal_root_states": row.optimal_root_states,
                    "character_weight": row.character_weight,
                    "weighted_score": row.weighted_score,
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
        columns=[
            "character_id",
            "minimum_cost",
            "observed_states",
            "matrix_states",
            "character_weight",
            "weighted_score",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "minimum_cost": row.minimum_cost,
                "observed_states": "|".join(row.observed_states),
                "matrix_states": "|".join(row.matrix_states),
                "character_weight": row.character_weight,
                "weighted_score": row.weighted_score,
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
            "weights_path": None
            if report.weights_path is None
            else str(report.weights_path),
            "total_weighted_score": report.total_weighted_score,
            "validation_warnings": [
                {
                    "code": warning.code,
                    "message": warning.message,
                    "details": warning.details,
                }
                for warning in report.validation_warnings
            ],
            "step_rows": [
                {
                    "character_id": row.character_id,
                    "minimum_cost": row.minimum_cost,
                    "observed_states": row.observed_states,
                    "matrix_states": row.matrix_states,
                    "character_weight": row.character_weight,
                    "weighted_score": row.weighted_score,
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
            "step_count",
            "character_weight",
            "weighted_score",
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
                "step_count": row.step_count,
                "character_weight": row.character_weight,
                "weighted_score": row.weighted_score,
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
            "weights_path": None
            if report.weights_path is None
            else str(report.weights_path),
            "total_weighted_score": report.total_weighted_score,
            "step_rows": [
                {
                    "character_id": row.character_id,
                    "derived_taxon_count": row.derived_taxon_count,
                    "gain_node": row.gain_node,
                    "gain_node_name": row.gain_node_name,
                    "gain_descendant_taxa": row.gain_descendant_taxa,
                    "total_losses": row.total_losses,
                    "impossible_state_warning": row.impossible_state_warning,
                    "step_count": row.step_count,
                    "character_weight": row.character_weight,
                    "weighted_score": row.weighted_score,
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
        columns=[
            "character_id",
            "derived_taxon_count",
            "gain_count",
            "root_state",
            "character_weight",
            "weighted_score",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "derived_taxon_count": row.derived_taxon_count,
                "gain_count": row.gain_count,
                "root_state": row.root_state,
                "character_weight": row.character_weight,
                "weighted_score": row.weighted_score,
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
            "weights_path": None
            if report.weights_path is None
            else str(report.weights_path),
            "total_weighted_score": report.total_weighted_score,
            "step_rows": [
                {
                    "character_id": row.character_id,
                    "derived_taxon_count": row.derived_taxon_count,
                    "gain_count": row.gain_count,
                    "root_state": row.root_state,
                    "character_weight": row.character_weight,
                    "weighted_score": row.weighted_score,
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
        columns=[
            "character_id",
            "step_count",
            "observed_states",
            "root_state",
            "character_weight",
            "weighted_score",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "step_count": row.step_count,
                "observed_states": "|".join(row.observed_states),
                "root_state": row.root_state,
                "character_weight": row.character_weight,
                "weighted_score": row.weighted_score,
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


def write_parsimony_ancestral_state_table(
    path: Path,
    report: ParsimonyReconstructionReport,
) -> Path:
    """Write one deterministic ancestral-state export table for a parsimony reconstruction."""
    return write_taxon_rows(
        path,
        columns=[
            "node_id",
            "clade_id",
            "character_id",
            "possible_states",
            "chosen_state",
            "method",
            "ambiguous",
        ],
        rows=[
            {
                "node_id": row.node_id,
                "clade_id": row.clade_id,
                "character_id": row.character_id,
                "possible_states": "|".join(row.possible_states),
                "chosen_state": row.chosen_state,
                "method": row.method,
                "ambiguous": row.ambiguous,
            }
            for row in report.ancestral_state_rows
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
            "branch_id",
            "character_id",
            "parent_node",
            "child_node",
            "child_node_name",
            "child_descendant_taxa",
            "change_from",
            "change_to",
            "ambiguous",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "character_id": row.character_id,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_node_name": row.child_node_name,
                "child_descendant_taxa": "|".join(row.child_descendant_taxa),
                "change_from": row.change_from,
                "change_to": row.change_to,
                "ambiguous": row.ambiguous,
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
            "weights_path": None
            if report.weights_path is None
            else str(report.weights_path),
            "total_weighted_score": report.total_weighted_score,
            "step_rows": [
                {
                    "character_id": row.character_id,
                    "step_count": row.step_count,
                    "observed_states": row.observed_states,
                    "root_state": row.root_state,
                    "character_weight": row.character_weight,
                    "weighted_score": row.weighted_score,
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
            "ancestral_state_rows": [
                {
                    "node_id": row.node_id,
                    "clade_id": row.clade_id,
                    "character_id": row.character_id,
                    "possible_states": row.possible_states,
                    "chosen_state": row.chosen_state,
                    "method": row.method,
                    "ambiguous": row.ambiguous,
                }
                for row in report.ancestral_state_rows
            ],
            "branch_change_rows": [
                {
                    "character_id": row.character_id,
                    "branch_id": row.branch_id,
                    "parent_node": row.parent_node,
                    "parent_state": row.parent_state,
                    "child_node": row.child_node,
                    "child_node_name": row.child_node_name,
                    "child_descendant_taxa": row.child_descendant_taxa,
                    "change_from": row.change_from,
                    "change_to": row.change_to,
                    "ambiguous": row.ambiguous,
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
    ancestral_states_path = write_parsimony_ancestral_state_table(
        out_dir / "ancestral_states.tsv",
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
        "ancestral_states_path": ancestral_states_path,
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


def write_parsimony_retention_index_table(
    path: Path,
    report: ParsimonyRetentionIndexReport,
) -> Path:
    """Write one deterministic per-character retention-index table."""
    return write_taxon_rows(
        path,
        columns=[
            "character_id",
            "character_kind",
            "observed_states",
            "minimum_possible_steps",
            "maximum_possible_steps",
            "observed_steps",
            "retention_index",
            "undefined_reason",
        ],
        rows=[
            {
                "character_id": row.character_id,
                "character_kind": row.character_kind,
                "observed_states": "|".join(row.observed_states),
                "minimum_possible_steps": row.minimum_possible_steps,
                "maximum_possible_steps": row.maximum_possible_steps,
                "observed_steps": row.observed_steps,
                "retention_index": row.retention_index,
                "undefined_reason": row.undefined_reason,
            }
            for row in report.character_rows
        ],
    )


def write_parsimony_retention_run_json(
    path: Path,
    report: ParsimonyRetentionIndexReport,
) -> Path:
    """Write one machine-readable retention-index payload."""
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
            "maximum_possible_steps_total": report.maximum_possible_steps_total,
            "observed_steps_total": report.observed_steps_total,
            "retention_index": report.retention_index,
            "undefined_reason": report.undefined_reason,
            "character_rows": [
                {
                    "character_id": row.character_id,
                    "character_kind": row.character_kind,
                    "observed_states": row.observed_states,
                    "minimum_possible_steps": row.minimum_possible_steps,
                    "maximum_possible_steps": row.maximum_possible_steps,
                    "observed_steps": row.observed_steps,
                    "retention_index": row.retention_index,
                    "undefined_reason": row.undefined_reason,
                }
                for row in report.character_rows
            ],
        },
    )


def write_parsimony_retention_artifacts(
    out_dir: Path,
    report: ParsimonyRetentionIndexReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one retention-index run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    indices_path = write_parsimony_retention_index_table(
        out_dir / "character_indices.tsv",
        report,
    )
    run_json_path = write_parsimony_retention_run_json(out_dir / "run.json", report)
    return {
        "indices_path": indices_path,
        "run_json_path": run_json_path,
    }


def write_parsimony_rescaled_consistency_index_table(
    path: Path,
    report: ParsimonyRescaledConsistencyIndexReport,
) -> Path:
    """Write one deterministic per-character rescaled-consistency table."""
    return write_taxon_rows(
        path,
        columns=["character_id", "ci", "ri", "rc", "undefined_reason"],
        rows=[
            {
                "character_id": row.character_id,
                "ci": row.ci,
                "ri": row.ri,
                "rc": row.rc,
                "undefined_reason": row.undefined_reason,
            }
            for row in report.character_rows
        ],
    )


def write_parsimony_rescaled_consistency_run_json(
    path: Path,
    report: ParsimonyRescaledConsistencyIndexReport,
) -> Path:
    """Write one machine-readable rescaled-consistency payload."""
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
            "ci": report.ci,
            "ri": report.ri,
            "rc": report.rc,
            "undefined_reason": report.undefined_reason,
            "character_rows": [
                {
                    "character_id": row.character_id,
                    "ci": row.ci,
                    "ri": row.ri,
                    "rc": row.rc,
                    "undefined_reason": row.undefined_reason,
                }
                for row in report.character_rows
            ],
        },
    )


def write_parsimony_rescaled_consistency_artifacts(
    out_dir: Path,
    report: ParsimonyRescaledConsistencyIndexReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one rescaled-consistency run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    indices_path = write_parsimony_rescaled_consistency_index_table(
        out_dir / "character_indices.tsv",
        report,
    )
    run_json_path = write_parsimony_rescaled_consistency_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "indices_path": indices_path,
        "run_json_path": run_json_path,
    }


def write_parsimony_bootstrap_replicate_scores_table(
    path: Path,
    report: ParsimonyBootstrapReport,
) -> Path:
    """Write one deterministic exact-search bootstrap replicate score table."""
    return write_taxon_rows(
        path,
        columns=[
            "replicate_index",
            "best_score",
            "optimal_tree_count",
            "tree_newick",
        ],
        rows=[
            {
                "replicate_index": row.replicate_index,
                "best_score": row.best_score,
                "optimal_tree_count": row.optimal_tree_count,
                "tree_newick": row.tree_newick,
            }
            for row in report.replicate_rows
        ],
    )


def write_parsimony_bootstrap_replicate_draws_table(
    path: Path,
    report: ParsimonyBootstrapReport,
) -> Path:
    """Write one deterministic bootstrap character-draw ledger."""
    return write_taxon_rows(
        path,
        columns=["replicate_index", "draw_index", "source_character_id"],
        rows=[
            {
                "replicate_index": row.replicate_index,
                "draw_index": draw_index,
                "source_character_id": source_character_id,
            }
            for row in report.replicate_rows
            for draw_index, source_character_id in enumerate(
                row.sampled_character_ids,
                start=1,
            )
        ],
    )


def write_parsimony_bootstrap_clade_support_table(
    path: Path,
    report: ParsimonyBootstrapReport,
) -> Path:
    """Write one deterministic reference-tree bootstrap support table."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "node_name",
            "descendant_taxa",
            "supporting_tree_count",
            "clade_frequency",
            "support_percent",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "node_name": row.node_name,
                "descendant_taxa": "|".join(row.descendant_taxa),
                "supporting_tree_count": row.supporting_tree_count,
                "clade_frequency": row.clade_frequency,
                "support_percent": row.support_percent,
            }
            for row in report.clade_support_rows
        ],
    )


def write_parsimony_bootstrap_run_json(
    path: Path,
    report: ParsimonyBootstrapReport,
) -> Path:
    """Write one machine-readable exact-search parsimony bootstrap payload."""
    return write_json_artifact(
        path,
        {
            "algorithm": report.algorithm,
            "method": report.method,
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
            "replicate_count": report.replicate_count,
            "random_seed": report.random_seed,
            "candidate_tree_count": report.candidate_tree_count,
            "max_exact_taxa": report.max_exact_taxa,
            "reference_score": report.reference_score,
            "reference_optimal_tree_count": report.reference_optimal_tree_count,
            "reference_tree_newick": report.reference_tree_newick,
            "replicate_rows": [
                {
                    "replicate_index": row.replicate_index,
                    "sampled_character_ids": row.sampled_character_ids,
                    "best_score": row.best_score,
                    "optimal_tree_count": row.optimal_tree_count,
                    "tree_newick": row.tree_newick,
                }
                for row in report.replicate_rows
            ],
            "clade_support_rows": [
                {
                    "branch_id": row.branch_id,
                    "node_name": row.node_name,
                    "descendant_taxa": row.descendant_taxa,
                    "supporting_tree_count": row.supporting_tree_count,
                    "clade_frequency": row.clade_frequency,
                    "support_percent": row.support_percent,
                }
                for row in report.clade_support_rows
            ],
        },
    )


def write_parsimony_bootstrap_artifacts(
    out_dir: Path,
    report: ParsimonyBootstrapReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one parsimony bootstrap run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    reference_tree_path = write_newick(
        out_dir / "reference_tree.nwk",
        loads_newick(report.reference_tree_newick),
    )
    replicate_trees_path = write_newick_tree_set(
        out_dir / "replicate_trees.nwk",
        [loads_newick(row.tree_newick) for row in report.replicate_rows],
    )
    replicate_scores_path = write_parsimony_bootstrap_replicate_scores_table(
        out_dir / "replicate_scores.tsv",
        report,
    )
    replicate_draws_path = write_parsimony_bootstrap_replicate_draws_table(
        out_dir / "replicate_draws.tsv",
        report,
    )
    clade_support_path = write_parsimony_bootstrap_clade_support_table(
        out_dir / "clade_support.tsv",
        report,
    )
    consensus_tree, clade_frequency_report = (
        compute_consensus_tree(replicate_trees_path)[0],
        compute_clade_frequency_table(replicate_trees_path),
    )
    consensus_tree_path = write_newick(out_dir / "consensus_tree.nwk", consensus_tree)
    clade_frequencies_path = write_clade_frequency_table(
        out_dir / "clade_frequencies.tsv",
        clade_frequency_report,
    )
    run_json_path = write_parsimony_bootstrap_run_json(out_dir / "run.json", report)
    return {
        "reference_tree_path": reference_tree_path,
        "replicate_trees_path": replicate_trees_path,
        "replicate_scores_path": replicate_scores_path,
        "replicate_draws_path": replicate_draws_path,
        "clade_support_path": clade_support_path,
        "consensus_tree_path": consensus_tree_path,
        "clade_frequencies_path": clade_frequencies_path,
        "run_json_path": run_json_path,
    }


def write_parsimony_jackknife_replicate_scores_table(
    path: Path,
    report: ParsimonyJackknifeReport,
) -> Path:
    """Write one deterministic exact-search jackknife replicate score table."""
    return write_taxon_rows(
        path,
        columns=[
            "replicate_index",
            "retained_character_count",
            "best_score",
            "optimal_tree_count",
            "tree_newick",
        ],
        rows=[
            {
                "replicate_index": row.replicate_index,
                "retained_character_count": row.retained_character_count,
                "best_score": row.best_score,
                "optimal_tree_count": row.optimal_tree_count,
                "tree_newick": row.tree_newick,
            }
            for row in report.replicate_rows
        ],
    )


def write_parsimony_jackknife_retained_characters_table(
    path: Path,
    report: ParsimonyJackknifeReport,
) -> Path:
    """Write one deterministic jackknife retained-character ledger."""
    return write_taxon_rows(
        path,
        columns=["replicate_index", "retained_index", "source_character_id"],
        rows=[
            {
                "replicate_index": row.replicate_index,
                "retained_index": retained_index,
                "source_character_id": source_character_id,
            }
            for row in report.replicate_rows
            for retained_index, source_character_id in enumerate(
                row.retained_character_ids,
                start=1,
            )
        ],
    )


def write_parsimony_jackknife_clade_support_table(
    path: Path,
    report: ParsimonyJackknifeReport,
) -> Path:
    """Write one deterministic reference-tree jackknife support table."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "node_name",
            "descendant_taxa",
            "supporting_tree_count",
            "clade_frequency",
            "support_percent",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "node_name": row.node_name,
                "descendant_taxa": "|".join(row.descendant_taxa),
                "supporting_tree_count": row.supporting_tree_count,
                "clade_frequency": row.clade_frequency,
                "support_percent": row.support_percent,
            }
            for row in report.clade_support_rows
        ],
    )


def write_parsimony_jackknife_run_json(
    path: Path,
    report: ParsimonyJackknifeReport,
) -> Path:
    """Write one machine-readable exact-search parsimony jackknife payload."""
    return write_json_artifact(
        path,
        {
            "algorithm": report.algorithm,
            "method": report.method,
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
            "replicate_count": report.replicate_count,
            "random_seed": report.random_seed,
            "retain_probability": report.retain_probability,
            "candidate_tree_count": report.candidate_tree_count,
            "max_exact_taxa": report.max_exact_taxa,
            "reference_score": report.reference_score,
            "reference_optimal_tree_count": report.reference_optimal_tree_count,
            "reference_tree_newick": report.reference_tree_newick,
            "replicate_rows": [
                {
                    "replicate_index": row.replicate_index,
                    "retained_character_count": row.retained_character_count,
                    "retained_character_ids": row.retained_character_ids,
                    "best_score": row.best_score,
                    "optimal_tree_count": row.optimal_tree_count,
                    "tree_newick": row.tree_newick,
                }
                for row in report.replicate_rows
            ],
            "clade_support_rows": [
                {
                    "branch_id": row.branch_id,
                    "node_name": row.node_name,
                    "descendant_taxa": row.descendant_taxa,
                    "supporting_tree_count": row.supporting_tree_count,
                    "clade_frequency": row.clade_frequency,
                    "support_percent": row.support_percent,
                }
                for row in report.clade_support_rows
            ],
        },
    )


def write_parsimony_jackknife_artifacts(
    out_dir: Path,
    report: ParsimonyJackknifeReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one parsimony jackknife run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    reference_tree_path = write_newick(
        out_dir / "reference_tree.nwk",
        loads_newick(report.reference_tree_newick),
    )
    replicate_trees_path = write_newick_tree_set(
        out_dir / "replicate_trees.nwk",
        [loads_newick(row.tree_newick) for row in report.replicate_rows],
    )
    replicate_scores_path = write_parsimony_jackknife_replicate_scores_table(
        out_dir / "replicate_scores.tsv",
        report,
    )
    retained_characters_path = write_parsimony_jackknife_retained_characters_table(
        out_dir / "retained_characters.tsv",
        report,
    )
    clade_support_path = write_parsimony_jackknife_clade_support_table(
        out_dir / "clade_support.tsv",
        report,
    )
    consensus_tree, _consensus_report = compute_consensus_tree(replicate_trees_path)
    consensus_tree_path = write_newick(out_dir / "consensus_tree.nwk", consensus_tree)
    clade_frequency_report = compute_clade_frequency_table(replicate_trees_path)
    clade_frequencies_path = write_clade_frequency_table(
        out_dir / "clade_frequencies.tsv",
        clade_frequency_report,
    )
    run_json_path = write_parsimony_jackknife_run_json(out_dir / "run.json", report)
    return {
        "reference_tree_path": reference_tree_path,
        "replicate_trees_path": replicate_trees_path,
        "replicate_scores_path": replicate_scores_path,
        "retained_characters_path": retained_characters_path,
        "clade_support_path": clade_support_path,
        "consensus_tree_path": consensus_tree_path,
        "clade_frequencies_path": clade_frequencies_path,
        "run_json_path": run_json_path,
    }


def write_parsimony_equal_best_scores_table(
    path: Path,
    report: ParsimonyEqualBestConsensusReport,
) -> Path:
    """Write one deterministic ledger of retained equally optimal trees."""
    return write_taxon_rows(
        path,
        columns=["tree_index", "total_score", "tree_newick"],
        rows=[
            {
                "tree_index": row.tree_index,
                "total_score": row.total_score,
                "tree_newick": row.tree_newick,
            }
            for row in report.equal_best_tree_rows
        ],
    )


def write_parsimony_equal_best_consensus_run_json(
    path: Path,
    report: ParsimonyEqualBestConsensusReport,
) -> Path:
    """Write one machine-readable exact equal-best consensus payload."""
    return write_json_artifact(
        path,
        {
            "algorithm": report.algorithm,
            "method": report.method,
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
            "candidate_tree_count": report.candidate_tree_count,
            "max_exact_taxa": report.max_exact_taxa,
            "max_retained_equal_best_trees": report.max_retained_equal_best_trees,
            "best_score": report.best_score,
            "equal_best_tree_count": report.equal_best_tree_count,
            "retained_equal_best_tree_count": report.retained_equal_best_tree_count,
            "retained_all_equal_best_trees": report.retained_all_equal_best_trees,
            "strict_consensus": None
            if report.strict_consensus is None
            else {
                "consensus_method": report.strict_consensus.consensus_method,
                "consensus_threshold": report.strict_consensus.consensus_threshold,
                "tree_count": report.strict_consensus.tree_count,
                "included_clade_count": report.strict_consensus.included_clade_count,
                "consensus_newick": report.strict_consensus.consensus_newick,
            },
            "majority_consensus": None
            if report.majority_consensus is None
            else {
                "consensus_method": report.majority_consensus.consensus_method,
                "consensus_threshold": report.majority_consensus.consensus_threshold,
                "tree_count": report.majority_consensus.tree_count,
                "included_clade_count": report.majority_consensus.included_clade_count,
                "consensus_newick": report.majority_consensus.consensus_newick,
            },
            "equal_best_tree_rows": [
                {
                    "tree_index": row.tree_index,
                    "total_score": row.total_score,
                    "tree_newick": row.tree_newick,
                }
                for row in report.equal_best_tree_rows
            ],
        },
    )


def write_parsimony_equal_best_consensus_artifacts(
    out_dir: Path,
    report: ParsimonyEqualBestConsensusReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one exact equal-best consensus run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    equal_best_trees = [loads_newick(row.tree_newick) for row in report.equal_best_tree_rows]
    equal_best_trees_path = write_newick_tree_set(
        out_dir / "equal_best_trees.nwk",
        equal_best_trees,
    )
    equal_best_scores_path = write_parsimony_equal_best_scores_table(
        out_dir / "equal_best_scores.tsv",
        report,
    )
    run_json_path = write_parsimony_equal_best_consensus_run_json(
        out_dir / "run.json",
        report,
    )
    outputs = {
        "equal_best_trees_path": equal_best_trees_path,
        "equal_best_scores_path": equal_best_scores_path,
        "run_json_path": run_json_path,
    }
    if report.strict_consensus is not None:
        outputs["strict_consensus_tree_path"] = write_newick(
            out_dir / "strict_consensus_tree.nwk",
            loads_newick(report.strict_consensus.consensus_newick),
        )
    if report.majority_consensus is not None:
        outputs["majority_consensus_tree_path"] = write_newick(
            out_dir / "majority_consensus_tree.nwk",
            loads_newick(report.majority_consensus.consensus_newick),
        )
        clade_frequency_report = compute_clade_frequency_table(equal_best_trees_path)
        outputs["clade_frequencies_path"] = write_clade_frequency_table(
            out_dir / "clade_frequencies.tsv",
            clade_frequency_report,
        )
    return outputs


def write_parsimony_placement_summary_table(
    path: Path,
    report: ParsimonyPlacementReport,
) -> Path:
    """Write one deterministic summary row per placed parsimony query."""
    return write_taxon_rows(
        path,
        columns=[
            "query_id",
            "character_count",
            "best_edge_id",
            "best_child_name",
            "best_descendant_taxa",
            "best_total_steps",
            "best_additional_steps",
            "best_total_weighted_score",
            "best_additional_weighted_score",
            "candidate_placement_count",
            "equally_best_placement_count",
        ],
        rows=[
            {
                "query_id": row.query_id,
                "character_count": row.character_count,
                "best_edge_id": row.best_edge_id,
                "best_child_name": row.best_child_name,
                "best_descendant_taxa": "|".join(row.best_descendant_taxa),
                "best_total_steps": row.best_total_steps,
                "best_additional_steps": row.best_additional_steps,
                "best_total_weighted_score": row.best_total_weighted_score,
                "best_additional_weighted_score": row.best_additional_weighted_score,
                "candidate_placement_count": row.candidate_placement_count,
                "equally_best_placement_count": row.equally_best_placement_count,
            }
            for row in report.query_summaries
        ],
    )


def write_parsimony_placement_alternative_table(
    path: Path,
    report: ParsimonyPlacementReport,
) -> Path:
    """Write one deterministic ledger of scored parsimony edge placements."""
    return write_taxon_rows(
        path,
        columns=[
            "query_id",
            "placement_rank",
            "edge_id",
            "child_name",
            "descendant_taxa",
            "total_steps",
            "additional_steps",
            "total_weighted_score",
            "additional_weighted_score",
            "is_equally_best",
        ],
        rows=[
            {
                "query_id": row.query_id,
                "placement_rank": row.placement_rank,
                "edge_id": row.edge_id,
                "child_name": row.child_name,
                "descendant_taxa": "|".join(row.descendant_taxa),
                "total_steps": row.total_steps,
                "additional_steps": row.additional_steps,
                "total_weighted_score": row.total_weighted_score,
                "additional_weighted_score": row.additional_weighted_score,
                "is_equally_best": row.is_equally_best,
            }
            for row in report.alternative_rows
        ],
    )


def write_parsimony_placement_tree_set(
    path: Path,
    report: ParsimonyPlacementReport,
) -> Path:
    """Write the complete equally optimal placement tree set across all queries."""
    return write_newick_tree_set(
        path,
        [
            loads_newick(row.placed_tree_newick)
            for row in report.alternative_rows
            if row.is_equally_best
        ],
    )


def write_parsimony_placement_run_json(
    path: Path,
    report: ParsimonyPlacementReport,
) -> Path:
    """Write one machine-readable parsimony placement payload."""
    return write_json_artifact(
        path,
        {
            "algorithm": report.algorithm,
            "method": report.method,
            "tree_path": None if report.tree_path is None else str(report.tree_path),
            "matrix_path": None
            if report.matrix_path is None
            else str(report.matrix_path),
            "query_matrix_path": None
            if report.query_matrix_path is None
            else str(report.query_matrix_path),
            "taxon_column": report.taxon_column,
            "reference_taxon_count": report.reference_taxon_count,
            "character_count": report.character_count,
            "edge_count": report.edge_count,
            "query_count": report.query_count,
            "reference_total_steps": report.reference_total_steps,
            "weights_path": None
            if report.weights_path is None
            else str(report.weights_path),
            "reference_total_weighted_score": report.reference_total_weighted_score,
            "query_summaries": [
                {
                    "query_id": row.query_id,
                    "character_count": row.character_count,
                    "best_edge_id": row.best_edge_id,
                    "best_child_name": row.best_child_name,
                    "best_descendant_taxa": row.best_descendant_taxa,
                    "best_total_steps": row.best_total_steps,
                    "best_additional_steps": row.best_additional_steps,
                    "best_total_weighted_score": row.best_total_weighted_score,
                    "best_additional_weighted_score": row.best_additional_weighted_score,
                    "candidate_placement_count": row.candidate_placement_count,
                    "equally_best_placement_count": row.equally_best_placement_count,
                    "selected_best_tree_newick": row.selected_best_tree_newick,
                }
                for row in report.query_summaries
            ],
            "alternative_rows": [
                {
                    "query_id": row.query_id,
                    "placement_rank": row.placement_rank,
                    "edge_id": row.edge_id,
                    "child_name": row.child_name,
                    "descendant_taxa": row.descendant_taxa,
                    "total_steps": row.total_steps,
                    "additional_steps": row.additional_steps,
                    "total_weighted_score": row.total_weighted_score,
                    "additional_weighted_score": row.additional_weighted_score,
                    "is_equally_best": row.is_equally_best,
                    "placed_tree_newick": row.placed_tree_newick,
                }
                for row in report.alternative_rows
            ],
        },
    )


def write_parsimony_placement_artifacts(
    out_dir: Path,
    report: ParsimonyPlacementReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one parsimony placement run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    equally_best_tree_path = write_parsimony_placement_tree_set(
        out_dir / "equally_best_placements.nwk",
        report,
    )
    summary_path = write_parsimony_placement_summary_table(
        out_dir / "summary.tsv",
        report,
    )
    alternative_path = write_parsimony_placement_alternative_table(
        out_dir / "alternative_placements.tsv",
        report,
    )
    run_json_path = write_parsimony_placement_run_json(out_dir / "run.json", report)
    return {
        "equally_best_tree_path": equally_best_tree_path,
        "summary_path": summary_path,
        "alternative_path": alternative_path,
        "run_json_path": run_json_path,
    }


def write_parsimony_bremer_support_table(
    path: Path,
    report: ParsimonyBremerSupportReport,
) -> Path:
    """Write one deterministic Bremer support table for a rooted reference tree."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "node_name",
            "descendant_taxa",
            "shortest_lacking_score",
            "decay_index",
            "shortest_lacking_tree_count",
            "shortest_lacking_tree_newick",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "node_name": row.node_name,
                "descendant_taxa": "|".join(row.descendant_taxa),
                "shortest_lacking_score": row.shortest_lacking_score,
                "decay_index": row.decay_index,
                "shortest_lacking_tree_count": row.shortest_lacking_tree_count,
                "shortest_lacking_tree_newick": row.shortest_lacking_tree_newick,
            }
            for row in report.bremer_rows
        ],
    )


def write_parsimony_bremer_support_run_json(
    path: Path,
    report: ParsimonyBremerSupportReport,
) -> Path:
    """Write one machine-readable exact Bremer support payload."""
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
            "candidate_tree_count": report.candidate_tree_count,
            "max_exact_taxa": report.max_exact_taxa,
            "reference_tree_newick": report.reference_tree_newick,
            "reference_tree_score": report.reference_tree_score,
            "optimal_score": report.optimal_score,
            "optimal_tree_count": report.optimal_tree_count,
            "optimal_tree_newick": report.optimal_tree_newick,
            "reference_tree_score_delta_from_optimal": (
                report.reference_tree_score_delta_from_optimal
            ),
            "reference_tree_is_optimal": report.reference_tree_is_optimal,
            "bremer_rows": [
                {
                    "branch_id": row.branch_id,
                    "node_name": row.node_name,
                    "descendant_taxa": row.descendant_taxa,
                    "shortest_lacking_score": row.shortest_lacking_score,
                    "decay_index": row.decay_index,
                    "shortest_lacking_tree_count": row.shortest_lacking_tree_count,
                    "shortest_lacking_tree_newick": row.shortest_lacking_tree_newick,
                }
                for row in report.bremer_rows
            ],
        },
    )


def write_parsimony_bremer_support_artifacts(
    out_dir: Path,
    report: ParsimonyBremerSupportReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one exact Bremer support run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    reference_tree_path = write_newick(
        out_dir / "reference_tree.nwk",
        loads_newick(report.reference_tree_newick),
    )
    optimal_tree_path = write_newick(
        out_dir / "optimal_tree.nwk",
        loads_newick(report.optimal_tree_newick),
    )
    bremer_support_path = write_parsimony_bremer_support_table(
        out_dir / "bremer_support.tsv",
        report,
    )
    run_json_path = write_parsimony_bremer_support_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "reference_tree_path": reference_tree_path,
        "optimal_tree_path": optimal_tree_path,
        "bremer_support_path": bremer_support_path,
        "run_json_path": run_json_path,
    }


def write_parsimony_nni_trace_table(
    path: Path,
    report: ParsimonyNniSearchReport,
) -> Path:
    """Write one deterministic rooted NNI search trace table."""
    return write_taxon_rows(
        path,
        columns=[
            "event_index",
            "event_kind",
            "iteration",
            "score_before",
            "score_after",
            "score_delta",
            "tree_before_newick",
            "tree_after_newick",
            "pivot_branch_id",
            "sibling_clade_id",
            "exchanged_clade_id",
            "stopping_reason",
        ],
        rows=[
            {
                "event_index": row.event_index,
                "event_kind": row.event_kind,
                "iteration": row.iteration,
                "score_before": row.score_before,
                "score_after": row.score_after,
                "score_delta": row.score_delta,
                "tree_before_newick": row.tree_before_newick,
                "tree_after_newick": row.tree_after_newick,
                "pivot_branch_id": row.pivot_branch_id,
                "sibling_clade_id": row.sibling_clade_id,
                "exchanged_clade_id": row.exchanged_clade_id,
                "stopping_reason": row.stopping_reason,
            }
            for row in report.trace_rows
        ],
    )


def write_parsimony_nni_run_json(
    path: Path,
    report: ParsimonyNniSearchReport,
) -> Path:
    """Write one machine-readable rooted NNI search payload."""
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
            "start_tree_newick": report.start_tree_newick,
            "start_score": report.start_score,
            "final_tree_newick": report.final_tree_newick,
            "final_score": report.final_score,
            "accepted_move_count": report.accepted_move_count,
            "evaluated_neighbor_count": report.evaluated_neighbor_count,
            "stopping_reason": report.stopping_reason,
            "trace_rows": [
                {
                    "event_index": row.event_index,
                    "event_kind": row.event_kind,
                    "iteration": row.iteration,
                    "score_before": row.score_before,
                    "score_after": row.score_after,
                    "score_delta": row.score_delta,
                    "tree_before_newick": row.tree_before_newick,
                    "tree_after_newick": row.tree_after_newick,
                    "pivot_branch_id": row.pivot_branch_id,
                    "sibling_clade_id": row.sibling_clade_id,
                    "exchanged_clade_id": row.exchanged_clade_id,
                    "stopping_reason": row.stopping_reason,
                }
                for row in report.trace_rows
            ],
        },
    )


def write_parsimony_nni_artifacts(
    out_dir: Path,
    report: ParsimonyNniSearchReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one rooted NNI search run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    start_tree_path = write_newick(
        out_dir / "start_tree.nwk",
        loads_newick(report.start_tree_newick),
    )
    final_tree_path = write_newick(
        out_dir / "final_tree.nwk",
        loads_newick(report.final_tree_newick),
    )
    trace_path = write_parsimony_nni_trace_table(out_dir / "search_trace.tsv", report)
    run_json_path = write_parsimony_nni_run_json(out_dir / "run.json", report)
    return {
        "start_tree_path": start_tree_path,
        "final_tree_path": final_tree_path,
        "trace_path": trace_path,
        "run_json_path": run_json_path,
    }


def write_parsimony_spr_trace_table(
    path: Path,
    report: ParsimonySprSearchReport,
) -> Path:
    """Write one deterministic rooted SPR search trace table."""
    return write_taxon_rows(
        path,
        columns=[
            "event_index",
            "event_kind",
            "iteration",
            "score_before",
            "score_after",
            "score_delta",
            "tree_before_newick",
            "tree_after_newick",
            "pruned_clade_id",
            "regraft_target_branch_id",
            "stopping_reason",
        ],
        rows=[
            {
                "event_index": row.event_index,
                "event_kind": row.event_kind,
                "iteration": row.iteration,
                "score_before": row.score_before,
                "score_after": row.score_after,
                "score_delta": row.score_delta,
                "tree_before_newick": row.tree_before_newick,
                "tree_after_newick": row.tree_after_newick,
                "pruned_clade_id": row.pruned_clade_id,
                "regraft_target_branch_id": row.regraft_target_branch_id,
                "stopping_reason": row.stopping_reason,
            }
            for row in report.trace_rows
        ],
    )


def write_parsimony_spr_run_json(
    path: Path,
    report: ParsimonySprSearchReport,
) -> Path:
    """Write one machine-readable rooted SPR search payload."""
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
            "start_tree_newick": report.start_tree_newick,
            "start_score": report.start_score,
            "final_tree_newick": report.final_tree_newick,
            "final_score": report.final_score,
            "accepted_move_count": report.accepted_move_count,
            "evaluated_neighbor_count": report.evaluated_neighbor_count,
            "stopping_reason": report.stopping_reason,
            "trace_rows": [
                {
                    "event_index": row.event_index,
                    "event_kind": row.event_kind,
                    "iteration": row.iteration,
                    "score_before": row.score_before,
                    "score_after": row.score_after,
                    "score_delta": row.score_delta,
                    "tree_before_newick": row.tree_before_newick,
                    "tree_after_newick": row.tree_after_newick,
                    "pruned_clade_id": row.pruned_clade_id,
                    "regraft_target_branch_id": row.regraft_target_branch_id,
                    "stopping_reason": row.stopping_reason,
                }
                for row in report.trace_rows
            ],
        },
    )


def write_parsimony_spr_artifacts(
    out_dir: Path,
    report: ParsimonySprSearchReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one rooted SPR search run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    start_tree_path = write_newick(
        out_dir / "start_tree.nwk",
        loads_newick(report.start_tree_newick),
    )
    final_tree_path = write_newick(
        out_dir / "final_tree.nwk",
        loads_newick(report.final_tree_newick),
    )
    trace_path = write_parsimony_spr_trace_table(out_dir / "search_trace.tsv", report)
    run_json_path = write_parsimony_spr_run_json(out_dir / "run.json", report)
    return {
        "start_tree_path": start_tree_path,
        "final_tree_path": final_tree_path,
        "trace_path": trace_path,
        "run_json_path": run_json_path,
    }


def write_parsimony_ratchet_cycle_table(
    path: Path,
    report: ParsimonyRatchetReport,
) -> Path:
    """Write one deterministic parsimony ratchet cycle ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "cycle_index",
            "start_score",
            "start_tree_newick",
            "perturbed_character_ids",
            "perturbation_factor",
            "perturbed_score",
            "perturbed_tree_newick",
            "perturbed_accepted_move_count",
            "restored_score",
            "restored_tree_newick",
            "restored_accepted_move_count",
            "best_score_after_cycle",
            "best_tree_after_cycle",
            "best_tree_improved",
        ],
        rows=[
            {
                "cycle_index": row.cycle_index,
                "start_score": row.start_score,
                "start_tree_newick": row.start_tree_newick,
                "perturbed_character_ids": "|".join(row.perturbed_character_ids),
                "perturbation_factor": row.perturbation_factor,
                "perturbed_score": row.perturbed_score,
                "perturbed_tree_newick": row.perturbed_tree_newick,
                "perturbed_accepted_move_count": row.perturbed_accepted_move_count,
                "restored_score": row.restored_score,
                "restored_tree_newick": row.restored_tree_newick,
                "restored_accepted_move_count": row.restored_accepted_move_count,
                "best_score_after_cycle": row.best_score_after_cycle,
                "best_tree_after_cycle": row.best_tree_after_cycle,
                "best_tree_improved": row.best_tree_improved,
            }
            for row in report.cycle_rows
        ],
    )


def write_parsimony_ratchet_best_tree_history_table(
    path: Path,
    report: ParsimonyRatchetReport,
) -> Path:
    """Write one deterministic ratchet best-tree history ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "history_index",
            "cycle_index",
            "best_score",
            "best_tree_newick",
        ],
        rows=[
            {
                "history_index": row.history_index,
                "cycle_index": row.cycle_index,
                "best_score": row.best_score,
                "best_tree_newick": row.best_tree_newick,
            }
            for row in report.best_tree_history_rows
        ],
    )


def write_parsimony_ratchet_run_json(
    path: Path,
    report: ParsimonyRatchetReport,
) -> Path:
    """Write one machine-readable parsimony ratchet payload."""
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
            "cycle_count": report.cycle_count,
            "random_seed": report.random_seed,
            "perturbed_character_count": report.perturbed_character_count,
            "perturbation_factor": report.perturbation_factor,
            "start_tree_newick": report.start_tree_newick,
            "start_score": report.start_score,
            "final_tree_newick": report.final_tree_newick,
            "final_score": report.final_score,
            "best_tree_newick": report.best_tree_newick,
            "best_score": report.best_score,
            "cycle_rows": [
                {
                    "cycle_index": row.cycle_index,
                    "start_score": row.start_score,
                    "start_tree_newick": row.start_tree_newick,
                    "perturbed_character_ids": row.perturbed_character_ids,
                    "perturbation_factor": row.perturbation_factor,
                    "perturbed_score": row.perturbed_score,
                    "perturbed_tree_newick": row.perturbed_tree_newick,
                    "perturbed_accepted_move_count": row.perturbed_accepted_move_count,
                    "restored_score": row.restored_score,
                    "restored_tree_newick": row.restored_tree_newick,
                    "restored_accepted_move_count": row.restored_accepted_move_count,
                    "best_score_after_cycle": row.best_score_after_cycle,
                    "best_tree_after_cycle": row.best_tree_after_cycle,
                    "best_tree_improved": row.best_tree_improved,
                }
                for row in report.cycle_rows
            ],
            "best_tree_history_rows": [
                {
                    "history_index": row.history_index,
                    "cycle_index": row.cycle_index,
                    "best_score": row.best_score,
                    "best_tree_newick": row.best_tree_newick,
                }
                for row in report.best_tree_history_rows
            ],
        },
    )


def write_parsimony_ratchet_artifacts(
    out_dir: Path,
    report: ParsimonyRatchetReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one parsimony ratchet run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    start_tree_path = write_newick(
        out_dir / "start_tree.nwk",
        loads_newick(report.start_tree_newick),
    )
    final_tree_path = write_newick(
        out_dir / "final_tree.nwk",
        loads_newick(report.final_tree_newick),
    )
    best_tree_path = write_newick(
        out_dir / "best_tree.nwk",
        loads_newick(report.best_tree_newick),
    )
    cycle_history_path = write_parsimony_ratchet_cycle_table(
        out_dir / "cycle_history.tsv",
        report,
    )
    best_tree_history_path = write_parsimony_ratchet_best_tree_history_table(
        out_dir / "best_tree_history.tsv",
        report,
    )
    run_json_path = write_parsimony_ratchet_run_json(out_dir / "run.json", report)
    return {
        "start_tree_path": start_tree_path,
        "final_tree_path": final_tree_path,
        "best_tree_path": best_tree_path,
        "cycle_history_path": cycle_history_path,
        "best_tree_history_path": best_tree_history_path,
        "run_json_path": run_json_path,
    }
