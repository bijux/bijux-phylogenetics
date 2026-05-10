# Weak signal instability bundle

Governed weak-signal case showing that significance can cross the 0.05 boundary across fixed-seed permutation reruns.

Claim surfaces:

- `weak-signal-instability-visible` — Weak phylogenetic signal instability remains visible (`matched_with_tolerance`)

Sources:

- `permutation-signal-instability` — Permutation-backed weak phylogenetic signal boundary via `bijux_phylogenetics.comparative.signal:compute_phylogenetic_signal_test`

Instability summary:

- alpha threshold: `0.05`
- min p-value: `0.03`
- max p-value: `0.08`
- runs below threshold: `7`
- runs at or above threshold: `13`

Limitations:

- This bundle proves one governed weak-signal boundary case and does not claim to exhaust every instability mode in comparative workflows.
