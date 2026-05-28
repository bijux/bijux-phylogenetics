from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.bayesian.beast import validate_fossil_calibration_table
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError, UnrootedTreeError

from .models import DatingCalibrationAnchor

_DATE_TOLERANCE = 1e-9


def load_fixed_dating_calibrations(
    tree_path: Path,
    calibration_path: Path,
) -> list[DatingCalibrationAnchor]:
    """Load fixed internal-node calibrations resolved onto one rooted tree."""
    validate_tree_path(tree_path, require_rooted=True)
    tree = load_tree(tree_path)
    tree.rooted = True
    if tree.rooted is not True:
        raise UnrootedTreeError(
            "dating calibrations require one rooted tree",
            code="dating_calibration_error",
        )
    validation_report = validate_fossil_calibration_table(tree_path, calibration_path)
    if validation_report.invalid_calibration_count:
        raise PhylogeneticsError(
            "dating calibration table contains invalid targets or age bounds",
            code="dating_calibration_error",
        )

    node_by_descendant_taxa = {
        frozenset(node.descendant_taxa): node for node in tree.iter_nodes(order="preorder")
    }
    anchors: list[DatingCalibrationAnchor] = []
    date_by_node_id: dict[str, float] = {}
    calibration_id_by_node_id: dict[str, str] = {}
    for calibration in validation_report.calibrations:
        if calibration.minimum_age is None or calibration.maximum_age is None:
            raise PhylogeneticsError(
                "smoothing cross-validation requires fixed calibrations with both minimum_age and maximum_age",
                code="dating_calibration_error",
            )
        if not math.isclose(
            calibration.minimum_age,
            calibration.maximum_age,
            rel_tol=0.0,
            abs_tol=_DATE_TOLERANCE,
        ):
            raise PhylogeneticsError(
                "smoothing cross-validation requires fixed calibrations with identical minimum_age and maximum_age",
                code="dating_calibration_error",
            )
        descendant_taxa = frozenset(calibration.taxa)
        node = node_by_descendant_taxa.get(descendant_taxa)
        if node is None or node.node_id is None:
            raise PhylogeneticsError(
                f"calibration '{calibration.calibration_id}' could not be mapped to one stable tree node",
                code="dating_calibration_error",
            )
        if node.is_leaf():
            raise PhylogeneticsError(
                "smoothing cross-validation requires internal-node calibrations rather than tip calibrations",
                code="dating_calibration_error",
            )
        fixed_date = float(calibration.minimum_age)
        existing_date = date_by_node_id.get(node.node_id)
        if existing_date is not None and not math.isclose(
            existing_date,
            fixed_date,
            rel_tol=0.0,
            abs_tol=_DATE_TOLERANCE,
        ):
            existing_calibration_id = calibration_id_by_node_id[node.node_id]
            raise PhylogeneticsError(
                "multiple fixed calibrations target the same node with different dates: "
                f"{existing_calibration_id}, {calibration.calibration_id}",
                code="dating_calibration_error",
            )
        if existing_date is not None:
            raise PhylogeneticsError(
                "smoothing cross-validation requires each fixed calibration to target a unique internal node",
                code="dating_calibration_error",
            )
        date_by_node_id[node.node_id] = fixed_date
        calibration_id_by_node_id[node.node_id] = calibration.calibration_id
        anchors.append(
            DatingCalibrationAnchor(
                calibration_id=calibration.calibration_id,
                target_kind=calibration.target_kind,
                target_label=calibration.target_label,
                descendant_taxa=sorted(descendant_taxa),
                node_id=node.node_id,
                node_kind="root" if node is tree.root else "internal",
                fixed_date=fixed_date,
            )
        )
    if len(anchors) < 2:
        raise PhylogeneticsError(
            "smoothing cross-validation requires at least two fixed internal-node calibrations",
            code="dating_calibration_error",
        )
    return sorted(anchors, key=lambda row: row.calibration_id)
