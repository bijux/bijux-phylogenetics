# Comparative Trust Boundaries

This study turns comparative failure surfaces and numerically fragile
boundary cases into governed Evidence IDs for `bijux-phylogenetics`.

It exists to keep trust boundaries visible:

- expected comparative input failures are treated as evidence, not as incidental exceptions
- weak-signal cases are rerun across fixed seeds so one lucky significance pass cannot overstate confidence
- OU identifiability warnings stay reviewer-visible on governed reference cases

Current bundles:

- `evidence-001` comparative input rejection
- `evidence-002` weak-signal instability
- `evidence-003` OU identifiability warnings
