from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from bijux_phylogenetics.datasets.rabies_method_sensitivity.slurm import (
    build_rabies_method_sensitivity_slurm_planning_report,
)
from bijux_phylogenetics.datasets.rabies_method_sensitivity.slurm.arrays import (
    build_rabies_method_sensitivity_slurm_array_strategy_report,
    write_rabies_method_sensitivity_slurm_array_members_table,
    write_rabies_method_sensitivity_slurm_array_partition_scripts,
    write_rabies_method_sensitivity_slurm_array_partitions_table,
    write_rabies_method_sensitivity_slurm_array_strategy_json,
)


def _write_bytes(path: Path, *, byte_count: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x" * byte_count)


def _build_planning_report(tmp_path: Path):
    auto_gap_root = tmp_path / "variants" / "auto-gap-threshold"
    auto_gappyout_root = tmp_path / "variants" / "auto-gappyout"
    ginsi_gap_root = tmp_path / "variants" / "ginsi-gap-threshold"
    ginsi_gappyout_root = tmp_path / "variants" / "ginsi-gappyout"
    for root, byte_count in (
        (auto_gap_root, 64_000),
        (auto_gappyout_root, 72_000),
        (ginsi_gap_root, 2_500_000),
        (ginsi_gappyout_root, 2_600_000),
    ):
        _write_bytes(root / "alignment.aln", byte_count=byte_count)
    workflow_report = SimpleNamespace(
        dataset=SimpleNamespace(
            dataset_id="rabies_method_sensitivity_panel",
            workflow_prefix="rabies-method-sensitivity-panel",
            taxon_count=9,
        ),
        task_records=(
            SimpleNamespace(variant_id="auto-gap-threshold", output_root=auto_gap_root),
            SimpleNamespace(variant_id="auto-gappyout", output_root=auto_gappyout_root),
            SimpleNamespace(
                variant_id="ginsi-gap-threshold", output_root=ginsi_gap_root
            ),
            SimpleNamespace(
                variant_id="ginsi-gappyout", output_root=ginsi_gappyout_root
            ),
        ),
        variant_runs=(
            SimpleNamespace(
                config=SimpleNamespace(
                    variant_id="auto-gap-threshold",
                    alignment_mode="auto",
                    trimming_mode="gap-threshold",
                ),
                alignment_length=1353,
                trimmed_alignment_length=1353,
            ),
            SimpleNamespace(
                config=SimpleNamespace(
                    variant_id="auto-gappyout",
                    alignment_mode="auto",
                    trimming_mode="gappyout",
                ),
                alignment_length=1353,
                trimmed_alignment_length=1353,
            ),
            SimpleNamespace(
                config=SimpleNamespace(
                    variant_id="ginsi-gap-threshold",
                    alignment_mode="ginsi",
                    trimming_mode="gap-threshold",
                ),
                alignment_length=1353,
                trimmed_alignment_length=1353,
            ),
            SimpleNamespace(
                config=SimpleNamespace(
                    variant_id="ginsi-gappyout",
                    alignment_mode="ginsi",
                    trimming_mode="gappyout",
                ),
                alignment_length=1353,
                trimmed_alignment_length=1353,
            ),
        ),
        iqtree_threads=1,
        bootstrap_replicates=1000,
    )
    return build_rabies_method_sensitivity_slurm_planning_report(workflow_report)


def test_build_rabies_method_sensitivity_slurm_array_strategy_report_groups_jobs(
    tmp_path: Path,
) -> None:
    planning_report = _build_planning_report(tmp_path)

    strategy = build_rabies_method_sensitivity_slurm_array_strategy_report(
        planning_report
    )

    assert strategy.dataset_id == "rabies_method_sensitivity_panel"
    assert strategy.partition_count == 2
    assert strategy.script_count == 2
    assert strategy.total_job_count == 4
    assert strategy.largest_partition_size == 2
    partition_ids = [row.partition_id for row in strategy.partitions]
    assert partition_ids == [
        "compact-mafft-auto-standard",
        "compact-mafft-ginsi-elevated",
    ]
    auto_partition, ginsi_partition = strategy.partitions
    assert auto_partition.variant_ids == ("auto-gap-threshold", "auto-gappyout")
    assert auto_partition.array_spec == "0-1"
    assert ginsi_partition.variant_ids == ("ginsi-gap-threshold", "ginsi-gappyout")
    assert all(member.script_path.endswith(".sbatch") for member in strategy.members)


def test_write_rabies_method_sensitivity_slurm_array_artifacts(tmp_path: Path) -> None:
    strategy = build_rabies_method_sensitivity_slurm_array_strategy_report(
        _build_planning_report(tmp_path)
    )

    partitions_path = write_rabies_method_sensitivity_slurm_array_partitions_table(
        tmp_path / "slurm-array-partitions.tsv",
        strategy,
    )
    members_path = write_rabies_method_sensitivity_slurm_array_members_table(
        tmp_path / "slurm-array-members.tsv",
        strategy,
    )
    strategy_path = write_rabies_method_sensitivity_slurm_array_strategy_json(
        tmp_path / "slurm-array-strategy.json",
        strategy,
    )
    scripts_root = write_rabies_method_sensitivity_slurm_array_partition_scripts(
        tmp_path / "slurm-arrays",
        strategy,
    )

    partitions_text = partitions_path.read_text(encoding="utf-8")
    assert "suggested_sbatch_command" in partitions_text
    assert "compact-mafft-auto-standard" in partitions_text
    members_text = members_path.read_text(encoding="utf-8")
    assert "array_index" in members_text
    assert "ginsi-gappyout" in members_text
    payload = json.loads(strategy_path.read_text(encoding="utf-8"))
    assert payload["partition_count"] == 2
    assert payload["largest_partition_size"] == 2
    auto_script = scripts_root / "compact-mafft-auto-standard.sbatch"
    assert auto_script.is_file()
    script_text = auto_script.read_text(encoding="utf-8")
    assert "#SBATCH --array=0-1" in script_text
    assert '--variant-id "${variant_id}"' in script_text
