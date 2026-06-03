from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .constraints import require_feasible_dating_calibration_constraints
from .models import DatingCalibrationAnchor

_DATE_TOLERANCE = 1e-9


def load_fixed_dating_calibrations(
    tree_path: Path,
    calibration_path: Path,
) -> list[DatingCalibrationAnchor]:
    """Load fixed internal-node calibrations resolved onto one rooted tree."""
    report = require_feasible_dating_calibration_constraints(
        tree_path, calibration_path
    )
    anchors: list[DatingCalibrationAnchor] = []
    date_by_node_id: dict[str, float] = {}
    calibration_id_by_node_id: dict[str, str] = {}
    for calibration in report.constraint_rows:
        if calibration.minimum_bound is None or calibration.maximum_bound is None:
            raise PhylogeneticsError(
                "smoothing cross-validation requires fixed calibrations with both minimum_age and maximum_age",
                code="dating_calibration_error",
            )
        if calibration.fixed_date is None:
            raise PhylogeneticsError(
                "smoothing cross-validation requires fixed calibrations with identical minimum_age and maximum_age",
                code="dating_calibration_error",
            )
        if calibration.node_kind == "tip":
            raise PhylogeneticsError(
                "smoothing cross-validation requires internal-node calibrations rather than tip calibrations",
                code="dating_calibration_error",
            )
        fixed_date = float(calibration.fixed_date)
        existing_date = date_by_node_id.get(calibration.node_id)
        if existing_date is not None and not math.isclose(
            existing_date,
            fixed_date,
            rel_tol=0.0,
            abs_tol=_DATE_TOLERANCE,
        ):
            existing_calibration_id = calibration_id_by_node_id[calibration.node_id]
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
        date_by_node_id[calibration.node_id] = fixed_date
        calibration_id_by_node_id[calibration.node_id] = calibration.calibration_id
        anchors.append(
            DatingCalibrationAnchor(
                calibration_id=calibration.calibration_id,
                target_kind=calibration.target_kind,
                target_label=calibration.target_label,
                descendant_taxa=calibration.descendant_taxa,
                node_id=calibration.node_id,
                node_kind=calibration.node_kind,
                fixed_date=fixed_date,
            )
        )
    if len(anchors) < 2:
        raise PhylogeneticsError(
            "smoothing cross-validation requires at least two fixed internal-node calibrations",
            code="dating_calibration_error",
        )
    return sorted(anchors, key=lambda row: row.calibration_id)
