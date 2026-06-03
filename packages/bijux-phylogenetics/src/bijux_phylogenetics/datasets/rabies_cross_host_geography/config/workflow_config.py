from __future__ import annotations

import json
from pathlib import Path

from ..models import (
    _ALIGNMENT_MODE,
    _BOOTSTRAP_CONSENSUS_THRESHOLD,
    _BOOTSTRAP_REPLICATES,
    _BOOTSTRAP_ROBUST_SUPPORT_THRESHOLD,
    _CLADE_METADATA_COLUMNS,
    _COMPARATIVE_BRANCH_LENGTH_FLOOR,
    _COMPARATIVE_FORMULA,
    _COMPARATIVE_RESPONSE,
    _DATASET_ID,
    _DATASET_LABEL,
    _GEOGRAPHY_MODEL,
    _GEOGRAPHY_TRAIT,
    _HOST_MODEL,
    _HOST_TRAIT,
    _IQTREE_SEED,
    _IQTREE_THREADS,
    _MAX_BOOTSTRAP_TREE_COUNT,
    _MAX_REPORT_TABLE_ROWS,
    _MEMORY_WARNING_THRESHOLD_BYTES,
    _OUTGROUP_TAXA,
    _SEQUENCE_TYPE,
    _TRIM_GAP_THRESHOLD,
    _TRIMMING_MODE,
    _WORKFLOW_PREFIX,
    _WORKFLOW_TIMEOUT_SECONDS,
    RabiesCrossHostGeographyPanelWorkflowConfig,
)
from .resources import _default_workflow_config_path


def _load_workflow_config(
    config_path: Path | None,
) -> RabiesCrossHostGeographyPanelWorkflowConfig:
    resolved_path = (
        _default_workflow_config_path()
        if config_path is None
        else Path(config_path).expanduser().resolve()
    )
    payload = json.loads(resolved_path.read_text(encoding="utf-8"))
    dataset_root = resolved_path.parent
    dataset_id = payload.get("dataset_id", _DATASET_ID)
    if dataset_id != _DATASET_ID:
        raise ValueError(
            f"workflow config dataset_id must be '{_DATASET_ID}', got '{dataset_id}'"
        )
    return RabiesCrossHostGeographyPanelWorkflowConfig(
        config_path=resolved_path,
        dataset_id=dataset_id,
        label=payload.get("label", _DATASET_LABEL),
        sequences_path=dataset_root / payload.get("sequences_path", "sequences.fasta"),
        metadata_path=dataset_root / payload.get("metadata_path", "metadata.csv"),
        centroids_path=dataset_root
        / payload.get("centroids_path", "region-centroids.csv"),
        sequence_type=payload.get("sequence_type", _SEQUENCE_TYPE),
        workflow_prefix=payload.get("workflow_prefix", _WORKFLOW_PREFIX),
        host_trait=payload.get("host_trait", _HOST_TRAIT),
        geography_trait=payload.get("geography_trait", _GEOGRAPHY_TRAIT),
        host_model=payload.get("host_model", _HOST_MODEL),
        geography_model=payload.get("geography_model", _GEOGRAPHY_MODEL),
        outgroup_taxa=tuple(payload.get("outgroup_taxa", list(_OUTGROUP_TAXA))),
        iqtree_seed=int(payload.get("iqtree_seed", _IQTREE_SEED)),
        iqtree_threads=int(payload.get("iqtree_threads", _IQTREE_THREADS)),
        bootstrap_replicates=int(
            payload.get("bootstrap_replicates", _BOOTSTRAP_REPLICATES)
        ),
        timeout_seconds=(
            None
            if payload.get("timeout_seconds", _WORKFLOW_TIMEOUT_SECONDS) is None
            else float(payload.get("timeout_seconds", _WORKFLOW_TIMEOUT_SECONDS))
        ),
        max_bootstrap_tree_count=(
            None
            if payload.get("max_bootstrap_tree_count", _MAX_BOOTSTRAP_TREE_COUNT)
            is None
            else int(payload.get("max_bootstrap_tree_count", _MAX_BOOTSTRAP_TREE_COUNT))
        ),
        max_report_table_rows=(
            None
            if payload.get("max_report_table_rows", _MAX_REPORT_TABLE_ROWS) is None
            else int(payload.get("max_report_table_rows", _MAX_REPORT_TABLE_ROWS))
        ),
        memory_warning_threshold_bytes=(
            None
            if payload.get(
                "memory_warning_threshold_bytes",
                _MEMORY_WARNING_THRESHOLD_BYTES,
            )
            is None
            else int(
                payload.get(
                    "memory_warning_threshold_bytes",
                    _MEMORY_WARNING_THRESHOLD_BYTES,
                )
            )
        ),
        alignment_mode=payload.get("alignment_mode", _ALIGNMENT_MODE),
        trimming_mode=payload.get("trimming_mode", _TRIMMING_MODE),
        trim_gap_threshold=float(
            payload.get("trim_gap_threshold", _TRIM_GAP_THRESHOLD)
        ),
        bootstrap_consensus_threshold=float(
            payload.get(
                "bootstrap_consensus_threshold",
                _BOOTSTRAP_CONSENSUS_THRESHOLD,
            )
        ),
        bootstrap_robust_support_threshold=float(
            payload.get(
                "bootstrap_robust_support_threshold",
                _BOOTSTRAP_ROBUST_SUPPORT_THRESHOLD,
            )
        ),
        clade_metadata_columns=tuple(
            payload.get("clade_metadata_columns", list(_CLADE_METADATA_COLUMNS))
        ),
        comparative_formula=payload.get("comparative_formula", _COMPARATIVE_FORMULA),
        comparative_response=payload.get("comparative_response", _COMPARATIVE_RESPONSE),
        comparative_branch_length_floor=float(
            payload.get(
                "comparative_branch_length_floor",
                _COMPARATIVE_BRANCH_LENGTH_FLOOR,
            )
        ),
    )
