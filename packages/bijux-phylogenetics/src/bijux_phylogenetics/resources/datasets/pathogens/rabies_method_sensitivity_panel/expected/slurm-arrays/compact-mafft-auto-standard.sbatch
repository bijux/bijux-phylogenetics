#!/usr/bin/env bash
set -euo pipefail

# partition_id: compact-mafft-auto-standard
# dataset_size_class: compact
# method_group: mafft-auto
# resource_class: standard
# array_spec: 0-1
# maximum_cpus_per_task: 1
# maximum_memory_mib: 1280
# maximum_wallclock_minutes: 27

#SBATCH --job-name=compact-mafft-auto-standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=1280M
#SBATCH --time=00:27:00
#SBATCH --array=0-1

variant_ids=(
  "auto-gap-threshold"
  "auto-gappyout"
)

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "SLURM_ARRAY_TASK_ID is required" >&2
  exit 2
fi

variant_id="${variant_ids[$SLURM_ARRAY_TASK_ID]}"
array_root="${BIJUX_PHYLOGENETICS_ARRAY_ROOT:-artifacts/rabies-method-sensitivity-arrays}"
output_root="${array_root}/compact-mafft-auto-standard/${variant_id}"
exec bijux-phylogenetics demo rabies-method-sensitivity-panel \
  --out "${output_root}" \
  --variant-id "${variant_id}" \
  --parallel-workers 1 \
  "$@"
