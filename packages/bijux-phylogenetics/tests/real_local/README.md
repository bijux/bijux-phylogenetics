# Real Engine Validation

The `real_local/` suite is the governed external-engine validation lane.

It is intentionally separate from the fast `engine_contract` tests that use
fixture executables and parser-only contract coverage.

Use this lane when you need evidence that Bijux still works with installed
external engines, not only with fake runners:

- `test_alignment_engine_validation_matrix.py` exercises MAFFT, trimAl,
  IQ-TREE, and FastTree on small real inputs and writes one matrix JSON artifact
  under the pytest temporary directory.
- `test_bayesian_engine_validation_matrix.py` exercises a real MrBayes run and
  records BEAST either from a live executable run when available or from the
  governed checked-in BEAST XML/log/tree corpus when the executable is absent.

Each matrix case records the reviewer-facing engine name, validation name,
validation mode, executable path when one was used, version text, command,
exit code, runtime, output paths, and output hashes.

Run only this lane with:

```bash
./.venv/bin/pytest tests/real_local -m engine_real
```
