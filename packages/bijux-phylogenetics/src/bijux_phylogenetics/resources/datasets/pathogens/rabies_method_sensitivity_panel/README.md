# Rabies Method-Sensitivity Panel

This packaged dataset contains one real rabies virus nucleoprotein panel built
for method-sensitivity workflow review. It starts from the same compact
cross-host panel used by the broader rabies demonstration surfaces, but the
owned workflow here focuses on one narrower question: which biological
conclusions remain stable when alignment, trimming, and inference-engine
choices are varied deliberately.

## Source

Source accessions:

- `MG458305` rabies virus isolate `RV108`
- `MG458304` rabies virus isolate `RV50`
- `PV641713` rabies virus isolate `PKDOW/PDB011`
- `PX845689` rabies virus isolate `4092/26/2025/Novosibirsk`
- `OQ693985` rabies virus isolate `RV_2277140521J_Mazovie_Raccoon_dog_Pol_2021`
- `PX845683` rabies virus isolate `4092/20/2025/Novosibirsk`
- `PX845681` rabies virus isolate `4092/18/2025/Novosibirsk`
- `PX845678` rabies virus isolate `4092/15/2025/Novosibirsk`
- `PX845676` rabies virus isolate `4092/13/2025/Novosibirsk`

All packaged sequences derive from the rabies virus nucleoprotein coding
region. Host and locality labels come from the accession source qualifiers
shipped with those records.

## Metadata Coding

The metadata table carries both raw provenance columns and grouped workflow
states:

- `host_species`: raw host label from the accession source feature
- `host_group`: grouped host state used by the broader rabies host workflow
- `country`: raw country or locality label from the accession source feature
- `region_group`: grouped macroregion used by the broader rabies geography
  workflow

The method-sensitivity workflow does not rerun host-switching or biogeography
analysis directly. It keeps the grouped metadata because those columns make the
reviewed clades interpretable when the final report describes which host- and
region-associated lineages remain stable under method choice.

## Governed Workflow

The packaged `workflow-config.json` file defines one fixed four-way
preprocessing matrix:

- MAFFT `auto` plus trimAl `gap-threshold`
- MAFFT `ginsi` plus trimAl `gap-threshold`
- MAFFT `auto` plus trimAl `gappyout`
- MAFFT `ginsi` plus trimAl `gappyout`

The same config also declares `parallel_workers: 2`, which is the governed
default for this batch workflow. Each variant still runs in its own isolated
output root, so the workflow can execute variants in parallel without sharing
one manifest prefix or one log file.

For each declared variant, the owned workflow reruns:

- MAFFT multiple-sequence alignment
- trimAl alignment trimming
- IQ-TREE model selection plus bootstrap-supported inference
- FastTree approximate inference on the same trimmed alignment
- outgroup rooting on `bat_chile_rv108` (`MG458305`)
- unrooted engine comparison and rooted engine comparison
- cross-variant rooted IQ-TREE comparison
- stable-clade and changed-conclusion ledgers
- final reviewer-facing HTML report

The goal is not to claim that every method is interchangeable. The workflow is
deliberately designed to show where the compact rabies panel is stable and
where the engine choice still changes internal clade conclusions.

## Contents

- `sequences.fasta`: raw rabies nucleoprotein coding sequences
- `metadata.csv`: combined host and geography metadata keyed by `taxon`
- `workflow-config.json`: one durable sensitivity matrix and execution config
- `expected/`: governed stable workflow outputs regenerated from the owned
  runtime surface

The governed expected bundle now also includes:

- `parallel-execution-summary.tsv`: one stable ledger of per-variant execution
  status, execution mode, and task log path
- `rabies-method-sensitivity-panel.run.json`: one workflow execution record
  over selected variants, task status, execution mode, and worker count
- `rabies-method-sensitivity.manifest.json`: one workflow manifest describing
  the parallel batch surface and its reviewer-facing outputs
- `slurm-job-plan.tsv`: one Slurm-ready per-variant job ledger with estimated
  CPUs, memory, wallclock, scratch, output size, and suggested `sbatch`
  options
- `slurm-estimation-assumptions.tsv`: one explicit sizing contract for the
  Slurm plan so reviewers can see why the planner reserved those resources
- `slurm-planning-summary.json`: one machine-readable planning summary over
  job count, total core-hours, maximum memory, wallclock, scratch, and output
  estimates
- `slurm-array-partitions.tsv`: one partition strategy that groups variants
  into real Slurm arrays by dataset-size class, method group, and resource
  class
- `slurm-array-members.tsv`: one task-level ledger mapping each variant to its
  array partition and array index
- `slurm-array-strategy.json`: one machine-readable export of the array
  partitioning contract
- `slurm-storage-categories.tsv`: one retained-storage ledger that separates
  workflow outputs, canonical logs, tree artifacts, posterior samples, and
  reviewer-facing reports
- `slurm-storage-variants.tsv`: one per-variant retained-storage ledger over
  those same categories
- `slurm-storage-report.json`: one machine-readable storage summary over total
  retained bytes, total retained MiB, category totals, and the largest variant
- `slurm-storage-report.html`: one reviewer-facing summary of the retained
  storage estimate before scaling the workflow further
- `slurm-output-explosion-checks.tsv`: one check-level ledger over the global
  and per-variant warnings that guard against retained-output growth
- `slurm-output-explosion-variants.tsv`: one per-variant risk ledger over
  retained output size, tree burden, posterior burden, and dominant output
  share
- `slurm-output-explosion-report.json`: one machine-readable warning summary
  over overall risk status, failed checks, and warning or high-risk variants
- `slurm-output-explosion-report.html`: one reviewer-facing summary of whether
  retained outputs, tree files, posterior samples, or report artifacts are
  starting to scale badly
- `slurm-tree-retention-checks.tsv`: one check-level ledger over the
  consistency checks behind safe tree thinning and compression policy
- `slurm-tree-retention-files.tsv`: one per-file retention ledger recording
  tree count, thinning policy, compression policy, and the retained-tree
  target for each tree-bearing file
- `slurm-tree-retention-policy.json`: one machine-readable summary over
  tree-set counts, thinning requirements, compression requirements, and the
  current overall policy status
- `slurm-tree-retention-policy.html`: one reviewer-facing summary of whether
  any retained tree sets need interval thinning or gzip compression before the
  workflow scales further
- `slurm-job-evidence.tsv`: one workflow-wide index over the per-job
  provenance packages written for each planned Slurm job
- `slurm-job-evidence-summary.json`: one machine-readable summary over the
  complete per-job evidence surface
- `slurm-job-evidence/<variant_id>/`: one self-contained debugging package for
  that job, including a copied task log, copied step manifests, one evidence
  JSON, and one reviewer-facing HTML summary
- `slurm-output-freshness.tsv`: one per-job stale-output ledger proving
  whether each planned job still matches the current packaged inputs and
  output-affecting settings
- `slurm-output-freshness-checks.tsv`: one check-level ledger over the exact
  input checksum and workflow-setting comparisons behind the freshness result
- `slurm-output-freshness.json`: one machine-readable freshness summary over
  fresh jobs, stale jobs, and failed freshness checks
- `slurm-failure-recovery-jobs.tsv`: one per-job recovery ledger over rerun
  eligibility, likely failure cause, and the next concrete recovery action
- `slurm-failure-recovery-partitions.tsv`: one partition-level recovery rollup
  over rerunnable jobs, blocked jobs, and recommended rerun scope
- `slurm-failure-recovery-report.json`: one machine-readable recovery summary
  over rerunnable jobs, blocked jobs, workflow state, and recovery partitions
- `slurm-failure-recovery-report.html`: one reviewer-facing summary of which
  jobs should be rerun, which ones are waiting on a live workflow, and the
  likely failure causes inferred from the governed task logs
- `slurm-merge-checks.tsv`: one batch-level ledger over the checks that decide
  whether the distributed variant outputs are globally merge-ready
- `slurm-merge-variants.tsv`: one per-variant merge ledger recording job
  status, freshness status, evidence presence, and whether that variant was
  included in the global merge
- `slurm-merge-report.json`: one machine-readable summary over merge readiness,
  merged variant count, failed merge checks, and the merged workflow totals
- `slurm-merge-report.html`: one reviewer-facing summary of the global merge
  decision with direct links back to the per-job evidence packages
- `slurm-job-status.tsv`: one per-job ledger classifying planned jobs as
  completed, failed, pending, or stale from the real execution evidence
- `slurm-partition-status.tsv`: one per-partition rollup of job-state counts
  for resumable batch review
- `slurm-workflow-status.json`: one machine-readable workflow-wide summary of
  completed, failed, pending, and stale job counts
- `slurm-arrays/*.sbatch`: one executable per-partition array script that runs
  the governed CLI with one selected variant per Slurm array task
- `reproducibility-checks.tsv`: one bundle-level audit ledger proving the
  current summary artifacts still match their manifests and current inputs
- `reproducibility-variants.tsv`: one per-variant provenance ledger over
  current file inventory, alignment lengths, and output digests
- `reproducibility-audit.json`: one machine-readable audit summary with the
  pass/fail result and all recorded checks
- `parallel-logs/*.log`: one orchestration log per variant, kept separate so
  concurrent batch execution remains inspectable
