from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict
import json
import math
from pathlib import Path

from bijux_phylogenetics.bayesian.beast import validate_fossil_calibration_table
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .models import (
    DatingCalibrationConstraintIssue,
    DatingCalibrationConstraintReport,
    DatingCalibrationConstraintRow,
    DatingCalibrationNodeWindowRow,
)

_DATE_TOLERANCE = 1e-9


def solve_dating_calibration_constraints(
    tree_path: Path,
    calibration_path: Path,
) -> DatingCalibrationConstraintReport:
    """Resolve dating calibrations onto one rooted tree and detect contradictory bounds."""
    validate_tree_path(tree_path, require_rooted=True)
    tree = load_tree(tree_path)
    tree.rooted = True
    validation_report = validate_fossil_calibration_table(tree_path, calibration_path)
    node_by_descendant_taxa = {
        frozenset(node.descendant_taxa): node
        for node in tree.iter_nodes(order="preorder")
    }

    issue_rows: list[DatingCalibrationConstraintIssue] = []
    constraint_rows: list[DatingCalibrationConstraintRow] = []
    window_state_by_node_id: dict[str, _NodeWindowState] = {}
    calibration_ids_by_node_id: dict[str, list[str]] = defaultdict(list)
    issue_codes_by_node_id: dict[str, set[str]] = defaultdict(set)
    issue_codes_by_calibration_id: dict[str, set[str]] = defaultdict(set)

    for issue in validation_report.issues:
        issue_rows.append(
            DatingCalibrationConstraintIssue(
                scope_kind="calibration",
                scope_id=issue.calibration_id,
                code=issue.code,
                message=issue.message,
                related_node_ids=[],
                related_calibration_ids=[issue.calibration_id],
            )
        )
        issue_codes_by_calibration_id[issue.calibration_id].add(issue.code)

    for calibration in validation_report.calibrations:
        if not calibration.valid:
            continue
        node = node_by_descendant_taxa.get(frozenset(calibration.taxa))
        if node is None or node.node_id is None:
            issue = DatingCalibrationConstraintIssue(
                scope_kind="calibration",
                scope_id=calibration.calibration_id,
                code="unmapped-valid-calibration",
                message="valid calibration could not be mapped to one stable tree node",
                related_node_ids=[],
                related_calibration_ids=[calibration.calibration_id],
            )
            issue_rows.append(issue)
            issue_codes_by_calibration_id[calibration.calibration_id].add(issue.code)
            continue
        node_kind = (
            "root" if node is tree.root else ("tip" if node.is_leaf() else "internal")
        )
        fixed_date = _fixed_date(
            minimum_bound=calibration.minimum_age,
            maximum_bound=calibration.maximum_age,
        )
        constraint_rows.append(
            DatingCalibrationConstraintRow(
                calibration_id=calibration.calibration_id,
                target_kind=calibration.target_kind,
                target_label=calibration.target_label,
                descendant_taxa=sorted(calibration.taxa),
                node_id=node.node_id,
                node_kind=node_kind,
                minimum_bound=calibration.minimum_age,
                maximum_bound=calibration.maximum_age,
                fixed_date=fixed_date,
                contradictory=False,
                issue_codes=[],
            )
        )
        calibration_ids_by_node_id[node.node_id].append(calibration.calibration_id)
        state = window_state_by_node_id.get(node.node_id)
        if state is None:
            state = _NodeWindowState(
                node_id=node.node_id,
                node_kind=node_kind,
                node_label=node.name,
                descendant_taxa=sorted(node.descendant_taxa),
                minimum_bound=calibration.minimum_age,
                maximum_bound=calibration.maximum_age,
            )
            window_state_by_node_id[node.node_id] = state
        else:
            if calibration.minimum_age is not None:
                state.minimum_bound = (
                    calibration.minimum_age
                    if state.minimum_bound is None
                    else max(state.minimum_bound, calibration.minimum_age)
                )
            if calibration.maximum_age is not None:
                state.maximum_bound = (
                    calibration.maximum_age
                    if state.maximum_bound is None
                    else min(state.maximum_bound, calibration.maximum_age)
                )

    for node_id, state in window_state_by_node_id.items():
        if (
            state.minimum_bound is not None
            and state.maximum_bound is not None
            and state.minimum_bound > (state.maximum_bound + _DATE_TOLERANCE)
        ):
            calibration_ids = sorted(calibration_ids_by_node_id[node_id])
            issue = DatingCalibrationConstraintIssue(
                scope_kind="node",
                scope_id=node_id,
                code="node-bounds-conflict",
                message="calibrations on one node impose one empty closed date window",
                related_node_ids=[node_id],
                related_calibration_ids=calibration_ids,
            )
            issue_rows.append(issue)
            issue_codes_by_node_id[node_id].add(issue.code)
            for calibration_id in calibration_ids:
                issue_codes_by_calibration_id[calibration_id].add(issue.code)

    effective_lower_by_node_id = _propagate_effective_lower_bounds(
        tree,
        window_state_by_node_id=window_state_by_node_id,
    )
    effective_upper_by_node_id = _propagate_effective_upper_bounds(
        tree,
        window_state_by_node_id=window_state_by_node_id,
    )

    for parent, child in tree.iter_edges():
        parent_node_id = parent.node_id or ""
        child_node_id = child.node_id or ""
        if (
            parent_node_id not in window_state_by_node_id
            and child_node_id not in window_state_by_node_id
        ):
            continue
        parent_lower = effective_lower_by_node_id[parent_node_id]
        child_upper = effective_upper_by_node_id[child_node_id]
        if parent_lower is None or child_upper is None:
            continue
        if parent_lower >= (child_upper - _DATE_TOLERANCE):
            related_node_ids = [parent_node_id, child_node_id]
            related_calibration_ids = sorted(
                set(calibration_ids_by_node_id.get(parent_node_id, []))
                | set(calibration_ids_by_node_id.get(child_node_id, []))
            )
            issue = DatingCalibrationConstraintIssue(
                scope_kind="edge",
                scope_id=f"{parent_node_id}->{child_node_id}",
                code="chronology-conflict",
                message="ancestor and descendant calibration windows admit no chronological ordering",
                related_node_ids=related_node_ids,
                related_calibration_ids=related_calibration_ids,
            )
            issue_rows.append(issue)
            for node_id in related_node_ids:
                issue_codes_by_node_id[node_id].add(issue.code)
            for calibration_id in related_calibration_ids:
                issue_codes_by_calibration_id[calibration_id].add(issue.code)

    contradictory_node_ids = sorted(
        node_id for node_id, codes in issue_codes_by_node_id.items() if codes
    )
    contradictory_calibration_ids = sorted(
        calibration_id
        for calibration_id, codes in issue_codes_by_calibration_id.items()
        if codes
    )

    constraint_rows = [
        DatingCalibrationConstraintRow(
            calibration_id=row.calibration_id,
            target_kind=row.target_kind,
            target_label=row.target_label,
            descendant_taxa=row.descendant_taxa,
            node_id=row.node_id,
            node_kind=row.node_kind,
            minimum_bound=row.minimum_bound,
            maximum_bound=row.maximum_bound,
            fixed_date=row.fixed_date,
            contradictory=row.calibration_id in contradictory_calibration_ids,
            issue_codes=sorted(issue_codes_by_calibration_id[row.calibration_id]),
        )
        for row in constraint_rows
    ]
    node_window_rows = [
        DatingCalibrationNodeWindowRow(
            node_id=node_id,
            node_kind=state.node_kind,
            node_label=state.node_label,
            descendant_taxa=state.descendant_taxa,
            calibration_ids=sorted(calibration_ids_by_node_id[node_id]),
            minimum_bound=state.minimum_bound,
            maximum_bound=state.maximum_bound,
            effective_lower_bound=effective_lower_by_node_id[node_id],
            effective_upper_bound=effective_upper_by_node_id[node_id],
            contradictory=node_id in contradictory_node_ids,
            issue_codes=sorted(issue_codes_by_node_id[node_id]),
        )
        for node_id, state in sorted(
            window_state_by_node_id.items(),
            key=lambda item: item[0],
        )
    ]

    return DatingCalibrationConstraintReport(
        tree_newick=dumps_newick(tree),
        taxa=sorted(tree.tip_names),
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        tree_path=str(tree_path),
        calibration_path=str(calibration_path),
        calibration_count=validation_report.calibration_count,
        valid_calibration_count=validation_report.valid_calibration_count,
        invalid_calibration_count=validation_report.invalid_calibration_count,
        resolved_calibration_count=len(constraint_rows),
        contradictory_calibration_count=len(contradictory_calibration_ids),
        contradictory_node_count=len(contradictory_node_ids),
        feasible=not issue_rows,
        constraint_rows=sorted(constraint_rows, key=lambda row: row.calibration_id),
        node_window_rows=node_window_rows,
        issue_rows=sorted(
            issue_rows,
            key=lambda row: (row.scope_kind, row.scope_id, row.code),
        ),
    )


def require_feasible_dating_calibration_constraints(
    tree_path: Path,
    calibration_path: Path,
) -> DatingCalibrationConstraintReport:
    """Resolve dating calibrations and raise when the set is invalid or contradictory."""
    report = solve_dating_calibration_constraints(tree_path, calibration_path)
    if report.feasible:
        return report
    first_issue = report.issue_rows[0]
    raise PhylogeneticsError(
        f"dating calibrations are infeasible: {first_issue.message}",
        code="dating_calibration_error",
    )


def write_dating_calibration_constraint_summary_tsv(
    path: Path,
    report: DatingCalibrationConstraintReport,
) -> Path:
    """Write one summary row for one dating calibration constraint solve."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "tree_path",
        "calibration_path",
        "tip_count",
        "internal_node_count",
        "calibration_count",
        "valid_calibration_count",
        "invalid_calibration_count",
        "resolved_calibration_count",
        "contradictory_calibration_count",
        "contradictory_node_count",
        "feasible",
    ]
    values = [
        report.tree_path,
        report.calibration_path,
        str(report.tip_count),
        str(report.internal_node_count),
        str(report.calibration_count),
        str(report.valid_calibration_count),
        str(report.invalid_calibration_count),
        str(report.resolved_calibration_count),
        str(report.contradictory_calibration_count),
        str(report.contradictory_node_count),
        str(report.feasible).lower(),
    ]
    path.write_text(
        "\n".join(["\t".join(columns), "\t".join(values)]) + "\n",
        encoding="utf-8",
    )
    return path


def write_dating_calibration_constraints_tsv(
    path: Path,
    report: DatingCalibrationConstraintReport,
) -> Path:
    """Write one resolved calibration row per valid input calibration."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "calibration_id",
        "target_kind",
        "target_label",
        "descendant_taxa",
        "node_id",
        "node_kind",
        "minimum_bound",
        "maximum_bound",
        "fixed_date",
        "contradictory",
        "issue_codes",
    ]
    lines = ["\t".join(columns)]
    lines.extend(
        "\t".join(
            [
                row.calibration_id,
                row.target_kind,
                row.target_label,
                "|".join(row.descendant_taxa),
                row.node_id,
                row.node_kind,
                _format_optional_float(row.minimum_bound),
                _format_optional_float(row.maximum_bound),
                _format_optional_float(row.fixed_date),
                str(row.contradictory).lower(),
                "|".join(row.issue_codes),
            ]
        )
        for row in report.constraint_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_dating_calibration_node_windows_tsv(
    path: Path,
    report: DatingCalibrationConstraintReport,
) -> Path:
    """Write one aggregated date window row per calibrated node."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "node_id",
        "node_kind",
        "node_label",
        "descendant_taxa",
        "calibration_ids",
        "minimum_bound",
        "maximum_bound",
        "effective_lower_bound",
        "effective_upper_bound",
        "contradictory",
        "issue_codes",
    ]
    lines = ["\t".join(columns)]
    lines.extend(
        "\t".join(
            [
                row.node_id,
                row.node_kind,
                row.node_label or "",
                "|".join(row.descendant_taxa),
                "|".join(row.calibration_ids),
                _format_optional_float(row.minimum_bound),
                _format_optional_float(row.maximum_bound),
                _format_optional_float(row.effective_lower_bound),
                _format_optional_float(row.effective_upper_bound),
                str(row.contradictory).lower(),
                "|".join(row.issue_codes),
            ]
        )
        for row in report.node_window_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_dating_calibration_issues_tsv(
    path: Path,
    report: DatingCalibrationConstraintReport,
) -> Path:
    """Write one contradiction or validation issue row."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "scope_kind",
        "scope_id",
        "code",
        "message",
        "related_node_ids",
        "related_calibration_ids",
    ]
    lines = ["\t".join(columns)]
    lines.extend(
        "\t".join(
            [
                row.scope_kind,
                row.scope_id,
                row.code,
                row.message,
                "|".join(row.related_node_ids),
                "|".join(row.related_calibration_ids),
            ]
        )
        for row in report.issue_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_dating_calibration_constraint_run_json(
    path: Path,
    report: DatingCalibrationConstraintReport,
) -> Path:
    """Write the full dating calibration constraint report as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_dating_calibration_constraint_artifacts(
    out_dir: Path,
    report: DatingCalibrationConstraintReport,
) -> dict[str, Path]:
    """Write governed artifact outputs for one dating calibration constraint solve."""
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = write_dating_calibration_constraint_summary_tsv(
        out_dir / "summary.tsv",
        report,
    )
    constraints_path = write_dating_calibration_constraints_tsv(
        out_dir / "constraints.tsv",
        report,
    )
    node_windows_path = write_dating_calibration_node_windows_tsv(
        out_dir / "node_windows.tsv",
        report,
    )
    issues_path = write_dating_calibration_issues_tsv(
        out_dir / "issues.tsv",
        report,
    )
    run_json_path = write_dating_calibration_constraint_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "summary_path": summary_path,
        "constraints_path": constraints_path,
        "node_windows_path": node_windows_path,
        "issues_path": issues_path,
        "run_json_path": run_json_path,
    }


class _NodeWindowState:
    def __init__(
        self,
        *,
        node_id: str,
        node_kind: str,
        node_label: str | None,
        descendant_taxa: list[str],
        minimum_bound: float | None,
        maximum_bound: float | None,
    ) -> None:
        self.node_id = node_id
        self.node_kind = node_kind
        self.node_label = node_label
        self.descendant_taxa = descendant_taxa
        self.minimum_bound = minimum_bound
        self.maximum_bound = maximum_bound


def _propagate_effective_lower_bounds(
    tree,
    *,
    window_state_by_node_id: dict[str, _NodeWindowState],
) -> dict[str, float | None]:
    propagated: dict[str, float | None] = {}

    def visit(node, inherited_lower: float | None) -> None:
        node_id = node.node_id or ""
        own_lower = None
        state = window_state_by_node_id.get(node_id)
        if state is not None:
            own_lower = state.minimum_bound
        effective_lower = own_lower
        if inherited_lower is not None:
            effective_lower = (
                inherited_lower
                if effective_lower is None
                else max(effective_lower, inherited_lower)
            )
        propagated[node_id] = effective_lower
        for child in node.children:
            visit(child, effective_lower)

    visit(tree.root, None)
    return propagated


def _propagate_effective_upper_bounds(
    tree,
    *,
    window_state_by_node_id: dict[str, _NodeWindowState],
) -> dict[str, float | None]:
    propagated: dict[str, float | None] = {}

    def visit(node) -> float | None:
        node_id = node.node_id or ""
        own_upper = None
        state = window_state_by_node_id.get(node_id)
        if state is not None:
            own_upper = state.maximum_bound
        descendant_uppers = [visit(child) for child in node.children]
        child_upper = (
            min(upper for upper in descendant_uppers if upper is not None)
            if any(upper is not None for upper in descendant_uppers)
            else None
        )
        effective_upper = own_upper
        if child_upper is not None:
            effective_upper = (
                child_upper
                if effective_upper is None
                else min(effective_upper, child_upper)
            )
        propagated[node_id] = effective_upper
        return effective_upper

    visit(tree.root)
    return propagated


def _fixed_date(
    *,
    minimum_bound: float | None,
    maximum_bound: float | None,
) -> float | None:
    if minimum_bound is None or maximum_bound is None:
        return None
    if math.isclose(
        minimum_bound,
        maximum_bound,
        rel_tol=0.0,
        abs_tol=_DATE_TOLERANCE,
    ):
        return minimum_bound
    return None


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".15g")
