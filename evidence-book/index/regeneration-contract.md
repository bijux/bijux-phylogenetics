# Regeneration Contract

This contract records which study files are governed sources, which are
generated durable outputs, and how each study is rerun.

- studies: `2`

## Primate Longevity Signal

- study id: `primate-longevity-signal`
- build script: `None`
- rerun command: `UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python -m bijux_phylogenetics.command_line evidence book build --study-id primate-longevity-signal`
- bundle ids: `evidence-001, evidence-002, evidence-003, evidence-004, evidence-005, evidence-006, evidence-007, evidence-008, evidence-009`
- source path count: `6`
- generated path count: `174`

## Primate PGLS And Signal

- study id: `primate-pgls-and-signal`
- build script: `None`
- rerun command: `UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python -m bijux_phylogenetics.command_line evidence book build --study-id primate-pgls-and-signal`
- bundle ids: `evidence-001, evidence-002, evidence-003, evidence-004, evidence-005, evidence-006, evidence-007, evidence-008, evidence-009, evidence-010`
- source path count: `4`
- generated path count: `144`
