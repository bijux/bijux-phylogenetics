# Regeneration Contract

This contract records which study files are governed sources, which are
generated durable outputs, and how each study is rerun.

- studies: `4`

## Comparative trust boundary evidence study

- study id: `comparative-trust-boundaries`
- build script: `evidence-book/studies/comparative-trust-boundaries/build_evidence.py`
- rerun command: `UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python evidence-book/studies/comparative-trust-boundaries/build_evidence.py`
- bundle ids: `evidence-001, evidence-002, evidence-003`
- source path count: `2`
- generated path count: `19`

## Primate Longevity Signal

- study id: `primate-longevity-signal`
- build script: `evidence-book/studies/primate-longevity-signal/build_evidence.py`
- rerun command: `UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python evidence-book/studies/primate-longevity-signal/build_evidence.py`
- bundle ids: `evidence-001, evidence-002, evidence-003, evidence-004, evidence-005, evidence-006, evidence-007, evidence-008, evidence-009`
- source path count: `3`
- generated path count: `98`

## Primate PGLS and signal evidence study

- study id: `primate-pgls-and-signal`
- build script: `evidence-book/studies/primate-pgls-and-signal/build_evidence.py`
- rerun command: `UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python evidence-book/studies/primate-pgls-and-signal/build_evidence.py`
- bundle ids: `evidence-001, evidence-002, evidence-003, evidence-004, evidence-005, evidence-006`
- source path count: `4`
- generated path count: `35`

## Taxon Trust

- study id: `taxon-trust`
- build script: `evidence-book/studies/taxon-trust/build_evidence.py`
- rerun command: `UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python evidence-book/studies/taxon-trust/build_evidence.py`
- bundle ids: `evidence-001`
- source path count: `1`
- generated path count: `16`
