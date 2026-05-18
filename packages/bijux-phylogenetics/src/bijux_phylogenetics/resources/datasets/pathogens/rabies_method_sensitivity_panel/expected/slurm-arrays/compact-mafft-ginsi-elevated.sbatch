#!/usr/bin/env bash
set -euo pipefail

# partition_id: compact-mafft-ginsi-elevated
# dataset_size_class: compact
# method_group: mafft-ginsi
# resource_class: elevated
# array_spec: 0-1
# maximum_cpus_per_task: 2
# maximum_memory_mib: 1536
# maximum_wallclock_minutes: 35

#SBATCH --job-name=compact-mafft-ginsi-elevated
#SBATCH --cpus-per-task=2
#SBATCH --mem=1536M
#SBATCH --time=00:35:00
#SBATCH --array=0-1

variant_ids=(
  "ginsi-gap-threshold"
  "ginsi-gappyout"
)

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "SLURM_ARRAY_TASK_ID is required" >&2
  exit 2
fi

variant_id="${variant_ids[$SLURM_ARRAY_TASK_ID]}"
array_root="${BIJUX_PHYLOGENETICS_ARRAY_ROOT:-artifacts/rabies-method-sensitivity-arrays}"
output_root="${array_root}/compact-mafft-ginsi-elevated/${variant_id}"
exec bijux-phylogenetics demo rabies-method-sensitivity-panel \
  --out "${output_root}" \
  --variant-id "${variant_id}" \
  --parallel-workers 1 \
  "$@"
