---
title: Common Workflows
audience: public
type: how-to
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-12
---

# Common Workflows

Typical public workflows include:

- validate and inspect a tree before downstream use
- trim and inspect an alignment before reporting
- run one command from raw FASTA to a supported inference bundle
- assemble aligned loci into a concatenated supermatrix
- audit a concatenated multi-locus matrix before inference
- run a comparative model and capture JSON plus report artifacts
- generate a review-ready figure package for a tree

The public workflow contract is that important outputs should be inspectable
after the command finishes.

When the goal is to fit a phylogenetic regression rather than only measure
signal, use `comparative pgls`. The command inspects the requested response and
predictors, fits generalized least squares on the phylogenetic covariance, and
reports coefficient estimates with standard errors, Student-t test statistics,
p-values, and explicit 95% confidence intervals.

```bash
bijux-phylogenetics comparative pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --response longevity \
  --predictors social_group_size brain_mass_g \
  --taxon-column species \
  --json
```

When `--lambda-value estimate` is used, Bijux first chooses the covariance
strength that best fits the data and then reports coefficient uncertainty for
that fitted covariance. When an external workflow already fixed one
phylogenetic covariance, pass that exact numeric lambda instead. That keeps
coefficient and p-value review tied to one covariance assumption rather than
silently mixing model-selection differences with coefficient-level inference.

When the response is binary rather than continuous, use `comparative logistic`.
This workflow keeps the comparative formula surface, but fits a binary
working-correlation approximation instead of reusing continuous-trait PGLS
output. It requires the response to be encoded explicitly as `0` and `1`,
supports the same predictor encoding rules as `comparative pgls`, and reports
Wald-normal coefficient uncertainty plus explicit separation-risk warnings.

```bash
bijux-phylogenetics comparative logistic \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --response sociality_present \
  --predictors brain_mass_g habitat \
  --taxon-column species \
  --lambda-value 1.0 \
  --coefficients-out artifacts/primates.logistic-coefficients.tsv \
  --fitted-out artifacts/primates.logistic-fitted.tsv \
  --excluded-taxa-out artifacts/primates.logistic-excluded.tsv \
  --json
```

The coefficient ledger keeps one row per fitted term with the estimate,
standard error, Wald test statistic, p-value, and 95% confidence interval. The
fitted ledger keeps one row per analyzed taxon with the observed binary
response, fitted probability, linear predictor, and raw residual. The
excluded-taxa ledger keeps the same explicit pruning contract as the other
comparative workflows. If the fitted probabilities collapse toward `0` or `1`,
or the working information matrix becomes singular enough to require
stabilization, the JSON result marks that as separation risk instead of
pretending the approximation is as stable as an ordinary continuous-trait fit.

When the goal is to compare competing comparative hypotheses rather than fit
just one formula, use `comparative model-selection`. This workflow keeps one
shared complete-case taxon set across every candidate formula, auto-detects
whether the shared response should use continuous-trait PGLS or the binary
working-correlation logistic surface, and then ranks the candidate formulas by
information criteria on that one fixed analysis set.

```bash
bijux-phylogenetics comparative model-selection \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --formula "sociality_present ~ brain_mass_g" \
  --formula "sociality_present ~ habitat" \
  --formula "sociality_present ~ brain_mass_g + habitat" \
  --taxon-column species \
  --lambda-value 1.0 \
  --ranking-out artifacts/primates.model-ranking.tsv \
  --pairwise-out artifacts/primates.model-pairwise.tsv \
  --excluded-taxa-out artifacts/primates.model-excluded.tsv \
  --json
```

The ranking ledger keeps one row per candidate with log-likelihood, AIC, AICc,
BIC, delta values, Akaike weight, selected-model status, and the encoded model
columns actually used in the fit. The pairwise ledger makes the comparison
contract explicit by marking whether each candidate pair is identical, nested,
or non-nested, and by preserving a likelihood-ratio statistic only where a
nested comparison is real. The excluded-taxa ledger records the shared
complete-case rule directly so reviewers can see which taxa were dropped before
any candidate formula was ranked.

When the goal is to review phylogenetic independent contrasts directly, use
`comparative contrasts`. The base workflow computes one standardized contrast
row per internal node for one numeric trait and can optionally fit one
regression-through-origin over matched predictor and response contrasts.

```bash
bijux-phylogenetics comparative contrasts \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait longevity \
  --predictor-trait brain_mass_g \
  --taxon-column species \
  --contrasts-out artifacts/primates.independent-contrasts.tsv \
  --regression-out artifacts/primates.independent-contrast-regression.tsv \
  --json
```

The contrast ledger keeps one row per internal node with the left and right
descendant taxa, the standardized contrast, the expected variance, the local
ancestral value, and the shared root estimate for the analyzed trait. When
`--predictor-trait` is supplied, the regression ledger keeps one row per
matched node with the predictor contrast, response contrast, fitted
through-origin response contrast, residual, and leverage fraction, while the
JSON metrics report the fitted slope and p-value explicitly.

When the goal is to measure whether one numeric trait shows phylogenetic
structure rather than fit a regression, use `comparative signal`. That
workflow reports Blomberg's K, Pagel's lambda, a permutation p-value for the
observed K value, and a likelihood-ratio-style p-value for the fitted lambda
against the zero-signal lambda boundary.

```bash
bijux-phylogenetics comparative signal \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait longevity \
  --taxon-column species \
  --permutations 199 \
  --summary-out artifacts/primates.signal-summary.tsv \
  --permutations-out artifacts/primates.signal-permutations.tsv \
  --json
```

The summary ledger keeps one row with the fitted K, lambda, log-likelihood
context, both p-values, and the permutation exceedance count. The permutation
ledger keeps one row per shuffled trait realization so reviewers can see the
null K distribution directly instead of only one final exceedance count.

When the goal is to fit a standalone continuous-trait evolution model rather
than a regression, use `comparative brownian`. This workflow keeps the
Brownian motion surface explicit by reporting the fitted root state,
evolutionary rate `sigma²`, log-likelihood, AIC, and AICc, while preserving the
taxa pruned before fitting.

```bash
bijux-phylogenetics comparative brownian \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait longevity \
  --taxon-column species \
  --summary-out artifacts/primates.brownian-summary.tsv \
  --excluded-taxa-out artifacts/primates.brownian-excluded.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, excluded taxon
count, fitted root state, `sigma²`, log-likelihood, AIC, AICc, and residual
diagnostic summary fields. The excluded-taxa ledger keeps one row per taxon
that was absent from the tree, absent from the trait table, missing the target
trait value, or pruned because the trait value was non-numeric. This makes the
Brownian fit auditable as a tree-plus-trait workflow rather than a detached
scalar estimate.

When the goal is to test whether a continuous trait is better explained by
constrained evolution toward an optimum, use `comparative ou`. This workflow
fits the stationary-root Ornstein-Uhlenbeck surface directly and reports the
fitted pull strength `alpha`, optimum `theta`, diffusion rate `sigma²`,
log-likelihood, AIC, and AICc, while preserving explicit pruning and
identifiability warnings.

```bash
bijux-phylogenetics comparative ou \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait longevity \
  --taxon-column species \
  --summary-out artifacts/primates.ou-summary.tsv \
  --excluded-taxa-out artifacts/primates.ou-excluded.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, excluded taxon
count, fitted `alpha`, `theta`, `sigma²`, log-likelihood, AIC, AICc, and
residual diagnostic summary fields. The same row also preserves the
convergence-status label and the count of OU identifiability warnings. The
excluded-taxa ledger keeps one row per taxon that was absent from the tree,
absent from the trait table, missing the target trait value, or pruned because
the trait value was non-numeric. This keeps OU review grounded in both the fit
statistics and the data-retention contract.

When the goal is to test whether trait evolution decelerates through time in an
early-burst or adaptive-radiation style pattern, use `comparative
early-burst`. This workflow fits the bounded rate-change surface directly,
compares the retained fit against Brownian and OU on the same pruned taxon set,
and flags weak identifiability when the likelihood surface stays broad or
slides onto the search boundary.

```bash
bijux-phylogenetics comparative early-burst \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait longevity \
  --taxon-column species \
  --summary-out artifacts/primates.early-burst-summary.tsv \
  --excluded-taxa-out artifacts/primates.early-burst-excluded.tsv \
  --comparison-out artifacts/primates.early-burst-comparison.tsv \
  --profile-out artifacts/primates.early-burst-profile.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, excluded taxon
count, fitted `rate_change`, root state, `sigma²`, log-likelihood, AIC, AICc,
the selected best model among Brownian/OU/early-burst, and the count of
identifiability warnings. The excluded-taxa ledger keeps one row per taxon that
was absent from the tree, absent from the trait table, missing the target trait
value, or pruned because the trait value was non-numeric. The comparison ledger
keeps both model-fit rows and likelihood-ratio rows so reviewers can see
whether early-burst is actually preferred over Brownian or OU. The profile
ledger keeps one fixed `rate_change` row per bounded likelihood evaluation so
weak identifiability can be reviewed directly rather than inferred from one
point estimate.

When the goal is to fit the same comparative regression across several response
traits and then inspect how those fitted traits still co-vary, use
`comparative multivariate`. This workflow keeps one shared complete-case taxon
set across every requested response and predictor, fits one comparative
regression per response on that exact taxon set, and then reports residual
covariance and residual trait-trait association explicitly.

```bash
bijux-phylogenetics comparative multivariate \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --responses longevity range_size \
  --predictors brain_mass_g social_group_size \
  --taxon-column species \
  --covariance-out artifacts/primates.multivariate-covariance.tsv \
  --associations-out artifacts/primates.multivariate-associations.tsv \
  --excluded-taxa-out artifacts/primates.multivariate-excluded.tsv \
  --json
```

The covariance ledger keeps one response-pair row with residual covariance,
residual correlation, pair count, and diagonal status. The association ledger
keeps one unique response-pair row with the same covariance and correlation
plus a correlation test statistic, p-value, and Fisher-style interval. The
excluded-taxa ledger makes the complete-case rule explicit by recording which
taxa were dropped because a required response or predictor column was blank.

When the goal is to detect whether one fitted comparative model is leaving
systematic residual structure inside particular subtrees, use
`comparative clade-residuals`. This workflow keeps the fitted taxon-level
residuals from one comparative model and then aggregates those residuals across
every internal non-root clade in the analyzed tree.

```bash
bijux-phylogenetics comparative clade-residuals \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --formula "longevity ~ brain_mass_g + habitat" \
  --taxon-column species \
  --lambda-value estimate \
  --taxa-out artifacts/primates.residual-taxa.tsv \
  --clades-out artifacts/primates.residual-clades.tsv \
  --json
```

The taxon ledger keeps one row per analyzed taxon with the observed value,
fitted value, raw residual, and standardized residual used for clade
aggregation. The clade ledger keeps one row per internal clade with its member
taxa, residual averages, residual sum-of-squares share, influence score,
residual-heavy flag, and rank. That makes it possible to distinguish one
isolated outlier taxon from a whole subtree that is carrying consistent model
misspecification.

When the question is not residual burden after one fit but whether the fitted
comparative conclusion itself depends on one major subtree, use
`comparative clade-stability`. This workflow fits one comparative model on the
baseline analyzed taxon set, derives major internal non-root clades from that
baseline tree, removes each candidate clade in turn, and refits the same model
on the retained taxa.

```bash
bijux-phylogenetics comparative clade-stability \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --formula "longevity ~ brain_mass_g + habitat" \
  --taxon-column species \
  --lambda-value estimate \
  --clades-out artifacts/primates.clade-stability.tsv \
  --terms-out artifacts/primates.clade-stability-terms.tsv \
  --json
```

The clade ledger keeps one row per candidate removal with the dropped clade,
retained taxon count, fit status, blocked-refit reason when one exists,
coefficient comparison count, sign-change count, significance-change count,
largest coefficient and p-value shifts, delta log-likelihood, and the final
influence rank. The term ledger then keeps one row per comparable coefficient
for each successful clade removal, including baseline and refit estimates,
baseline and refit p-values, and explicit sign-change and
significance-change flags.

This surface matters when a biological interpretation looks fragile. A large
or influential clade can drive coefficient direction, apparent significance,
or fit quality without leaving obviously extreme single-taxon residuals. The
workflow keeps those subtree dependencies explicit and preserves blocked rows
instead of silently skipping removals that would collapse the remaining fit.

When the question is not whether one clade drives the fit but whether topology
uncertainty across a posterior or bootstrap tree set changes the comparative
conclusion, use `comparative posterior-pgls`. This workflow keeps one
continuous-trait PGLS formula fixed, applies it to every retained tree in a
tree set, and then summarizes how the coefficients and support calls vary
across those trees.

```bash
bijux-phylogenetics comparative posterior-pgls \
  artifacts/primates.posterior.trees \
  artifacts/primates.csv \
  --formula "longevity ~ social_group_size" \
  --taxon-column species \
  --lambda-value estimate \
  --burnin-fraction 0.25 \
  --significance-threshold 0.05 \
  --trees-out artifacts/primates.posterior-pgls-trees.tsv \
  --coefficients-out artifacts/primates.posterior-pgls-coefficients.tsv \
  --summary-out artifacts/primates.posterior-pgls-summary.tsv \
  --json
```

The per-tree ledger keeps one row per retained tree with its post-burn-in
position, topology identifiers, fitted lambda, and log-likelihood. The
coefficient ledger then keeps one row per coefficient per retained tree,
including estimate, p-value, direction, and whether the coefficient met the
chosen support threshold on that tree. The summary ledger collapses those
coefficient rows into reviewer-facing distributions with empirical estimate
ranges, direction consistency, support fractions, and a conclusion-stability
classification such as `stable_supported`, `stable_unsupported`,
`mixed_support`, or `direction_conflict`.

This surface matters when one MCC tree or one consensus tree would hide the
fact that coefficient support is topology-sensitive. It lets users propagate
tree uncertainty into the comparative conclusion directly instead of treating
the posterior tree set and the trait model as separate review steps.

When the covariance assumption itself must stay fixed to Brownian shared branch
lengths, use `comparative brownian-pgls`. This keeps the regression surface
separate from Pagel-lambda fitting and writes an explicit covariance ledger
when requested.

```bash
bijux-phylogenetics comparative brownian-pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --response longevity \
  --predictors social_group_size \
  --taxon-column species \
  --covariance-out artifacts/primates.brownian-covariance.tsv \
  --json
```

That workflow fits PGLS under the raw Brownian shared-path covariance implied
by the rooted tree, checks that the unstabilized covariance is positive
definite, and records one pairwise row per taxon pair with shared path length
and root-to-tip depths. Rooted ultrametric and rooted non-ultrametric trees are
both supported as long as the Brownian covariance stays valid. If zero or
negative branch lengths make the covariance invalid, the workflow fails
explicitly instead of silently regularizing the tree into a different
scientific assumption.

When the residual structure is expected to reflect pull toward an optimum rather
than unconstrained Brownian diffusion, use `comparative ou-pgls`. The workflow
supports either a fixed positive `--alpha` or a governed `estimate` mode that
profiles alpha across the bounded search grid already used by the OU trait-model
surface.

```bash
bijux-phylogenetics comparative ou-pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --response longevity \
  --predictors social_group_size \
  --taxon-column species \
  --alpha estimate \
  --covariance-out artifacts/primates.ou-covariance.tsv \
  --alpha-profile-out artifacts/primates.ou-alpha-profile.tsv \
  --json
```

That workflow fits the comparative regression under stationary-root OU
covariance, records the fitted log-likelihood and AIC, and writes two optional
review ledgers. The covariance ledger keeps one pairwise row per taxon pair
with the implied OU covariance, shared path length, and root-depth context. The
alpha profile ledger keeps one row per candidate alpha value, the
log-likelihood drop from the best-supported fit, and whether the row stays
inside the likelihood-ratio-supported 95% interval. Fixed-alpha review and
estimated-alpha review therefore stay explicit instead of being folded into one
opaque regression number.

When the main question is how strongly the regression residuals prefer
phylogenetic covariance, write the governed lambda profile in the same run.

```bash
bijux-phylogenetics comparative pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --response longevity \
  --predictors social_group_size \
  --taxon-column species \
  --lambda-value estimate \
  --lambda-profile-out artifacts/primates.lambda-profile.tsv \
  --json
```

The profile ledger keeps one row per candidate lambda value across the bounded
search surface, records the log-likelihood drop from the best fit, and marks
which rows stay inside the likelihood-ratio-supported 95% confidence interval.

When the main question is whether the requested biological formula expands into
the encoded predictors you actually expect, use formula syntax directly and
write the governed model matrix. Intercept-free formulas use the standard
comparative spellings `0 + ...` or `... - 1`.

```bash
bijux-phylogenetics comparative pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --formula 'longevity ~ 0 + social_group_size * habitat' \
  --taxon-column species \
  --model-matrix-out artifacts/primates.model-matrix.tsv \
  --json
```

The written matrix keeps one row per analyzed taxon, the response value, and
one encoded column per fitted predictor term. That makes it possible to review
continuous predictors, categorical indicator columns, interaction columns, and
intercept inclusion before interpreting the fitted coefficients.

When categorical predictors are present and the question is how each biological
group was encoded and interpreted, write the governed categorical-contrast
ledger in the same run.

```bash
bijux-phylogenetics comparative pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --response longevity \
  --predictors social_group_size habitat \
  --taxon-column species \
  --categorical-contrasts-out artifacts/primates.categorical-contrasts.tsv \
  --json
```

That ledger keeps one row for the baseline group when treatment coding is used
and one row for every estimated non-baseline group coefficient. It also records
missing-category taxa explicitly, so a dropped or blank category value does not
disappear into a generic exclusion count.

When interaction terms are present and the question is how effect modification
was encoded and estimated, write the governed interaction-coefficient ledger in
the same run.

```bash
bijux-phylogenetics comparative pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --formula 'longevity ~ social_group_size * habitat' \
  --taxon-column species \
  --interaction-coefficients-out artifacts/primates.interaction-coefficients.tsv \
  --json
```

That ledger keeps one row per fitted interaction coefficient, records whether
the interaction is continuous-by-continuous, continuous-by-categorical, or
categorical-by-categorical, and preserves any omitted treatment-coded baseline
levels so the effect-modification interpretation stays explicit.

When a tree is already inferred and you need to root it on one known outgroup
taxon or one expected outgroup clade, use `topology root-outgroup`. The
command writes the rooted tree and can also emit a one-row TSV report that
records which requested taxa were found, which requested taxa were absent,
whether the matched outgroup is monophyletic in the input tree, which extra
taxa fall inside the matched outgroup MRCA when it is not monophyletic, and
which taxa end up isolated on the rooted outgroup side.

```bash
bijux-phylogenetics topology root-outgroup \
  artifacts/mammals.unrooted.nwk \
  --taxa Ornithorhynchus_anatinus Tachyglossus_aculeatus \
  --out artifacts/mammals.rooted.nwk \
  --report-out artifacts/mammals.rooting.tsv \
  --json
```

If the requested outgroup taxa are missing, the TSV and JSON report that
explicitly instead of silently dropping them. If the requested outgroup taxa
are not monophyletic in the input tree, the workflow still records the rooted
tree but also reports the MRCA spillover taxa and a warning that the root does
not cleanly isolate the requested outgroup as one coherent clade.

When no explicit outgroup is available, use `topology reroot-midpoint` for an
exploratory rooted tree. The command writes the rerooted tree and can also
emit a one-row TSV report that records the anchor tip pair that defined the
selected midpoint path, the total tip-to-tip path length, the midpoint
distance from the anchor tip, whether the midpoint landed on an original node
or within an original branch, and which taxa ended up on each side of the new
root.

```bash
bijux-phylogenetics topology reroot-midpoint \
  artifacts/mammals.unrooted.nwk \
  --out artifacts/mammals.midpoint-rooted.nwk \
  --report-out artifacts/mammals.midpoint-rooting.tsv \
  --json
```

The midpoint report also records whether the input tree was suitable for
straightforward midpoint interpretation. Trees that are not strictly
bifurcating are still rerooted, but the report marks them as exploratory and
adds an explicit warning so downstream review does not overclaim the result.

When the goal is to inspect every clade directly instead of transforming the
tree, use `topology clades`. That command writes one row per node-derived
clade, including tips, internal clades, and the root. Each row keeps the
member taxa, any parsed support label, the incoming branch length, root depth,
descendant tip depths, and `node_age` when the tree is branch-length complete
and ultrametric.

```bash
bijux-phylogenetics topology clades \
  artifacts/mammals.supported.nwk \
  --metadata artifacts/mammals.metadata.tsv \
  --metadata-column species \
  --metadata-column location \
  --out artifacts/mammals.clades.tsv \
  --json
```

When `--metadata` is supplied, Bijux treats the table as a taxon-keyed
metadata or trait table and flattens the requested columns into per-clade
review fields such as `species_values`, `species_distinct_values`, and
`species_missing_taxa`. That keeps trait inspection tied directly to clade
membership instead of forcing a separate join step.

When the goal is to inspect tree shape directly rather than individual clades,
use `topology shape`. That command writes one review row with balance and shape
metrics including Sackin imbalance, Colless imbalance where defined, cherry
count, tree height in edges, branch-length tree height where available, and the
governed `balanced`, `skewed`, or `ladderized` summary.

```bash
bijux-phylogenetics topology shape \
  artifacts/mammals.supported.nwk \
  --out artifacts/mammals.shape.tsv \
  --json
```

The JSON output also preserves whether the tree is star-like, comb-like, or
unusually imbalanced. That keeps obvious ladderization or star-topology risks
visible without forcing users to infer them from raw node depths alone.

When the goal is to inspect branch-length patterns directly, use
`topology branch-lengths`. That command writes one row per non-root branch and
summarizes minimum, maximum, mean, and median branch length alongside explicit
zero-length, negative-length, and outlier flags.

```bash
bijux-phylogenetics topology branch-lengths \
  artifacts/mammals.supported.nwk \
  --out artifacts/mammals.branch-lengths.tsv \
  --json
```

This surface matters when the question is not just whether branch lengths exist
but whether they contain odd zeros, extreme long branches, or other scale
distortions that can mislead downstream interpretation.

When you already have two inferred trees and need an explicit topology
distance, use `compare`. The command now supports rooted and unrooted
Robinson-Foulds review directly, and it records which taxa were shared versus
present on only one side before the distance is computed.

```bash
bijux-phylogenetics compare \
  artifacts/mammals.iqtree.nwk \
  artifacts/mammals.fasttree.nwk \
  --rf-mode unrooted \
  --json
```

Use rooted RF when root placement is part of the scientific claim. Use
unrooted RF when the question is only whether the same splits were recovered
regardless of root placement. By default the command prunes both trees to
their shared taxa before computing RF distance, which is appropriate for
reviewing partially overlapping outputs. If taxon-set mismatch itself should
fail the comparison, add `--taxon-overlap-policy require-identical`.

```bash
bijux-phylogenetics compare \
  artifacts/mammals.reference.nwk \
  artifacts/mammals.candidate.nwk \
  --rf-mode rooted \
  --taxon-overlap-policy require-identical \
  --json
```

When the first question is whether two trees can be reduced to the same taxon
set safely before any deeper comparison, use `compare prune`. That command
writes the two pruned trees plus a governed pruning review bundle in one output
directory.

```bash
bijux-phylogenetics compare prune \
  artifacts/mammals.reference.nwk \
  artifacts/mammals.candidate.nwk \
  --out artifacts/mammals.shared-taxa \
  --json
```

The pruning bundle keeps the evidence explicit:
- `left-shared.nwk` and `right-shared.nwk` are the retained shared-taxon trees
- `shared-taxa-pruning.tsv` records one row per input tree with retained and
  removed taxa, branch-length audit fields, and information-loss fields
- `shared-taxa-removed.tsv` records one row per removed taxon with the side and
  removal reason
- `shared-taxa-comparison.tsv` compares the retained trees directly so the
  post-pruning topology review is durable instead of implicit

The JSON payload also preserves the full left and right pruning audits plus a
`post_pruning_comparison` report. That makes it possible to inspect branch
length preservation, removed taxa, and retained-tree topology without
reconstructing the review from separate commands.

When branch lengths matter as well as topology, use
`compare branch-lengths`. That surface preserves the per-split length table
for shared rooted clades and also computes Felsenstein branch-score distance on
the union of unrooted splits. Missing splits count as zero-length branches for
the score calculation, which makes topology disagreement contribute directly to
the final distance instead of disappearing from the branch-length review.

```bash
bijux-phylogenetics compare branch-lengths \
  artifacts/mammals.reference.nwk \
  artifacts/mammals.candidate.nwk \
  --json
```

By default the command prunes both trees to their shared taxa before computing
branch-score distance. If taxon-set mismatch itself should fail the review, add
`--taxon-overlap-policy require-identical`. When a matched split is present but
lacks a branch length on one side, the per-split ledger records that missing
length explicitly and the branch-score summary becomes unavailable instead of
silently treating the missing value as zero.

When the main question is whether topology conflicts are serious or only weakly
supported, use `compare support`. That surface combines clade presence with the
support values parsed from each tree and writes one support-aware conflict
ledger with `--out`.

```bash
bijux-phylogenetics compare support \
  artifacts/mammals.reference.nwk \
  artifacts/mammals.candidate.nwk \
  --out artifacts/mammals.support-conflicts.tsv \
  --json
```

The JSON and TSV outputs separate three situations explicitly:
- shared clades, where both trees carry the same clade and Bijux reports the
  normalized support delta
- high-support conflicts, where a clade present in only one tree still carries
  normalized support of at least `0.9`
- low-support disagreements, where a conflicting clade is present with support
  below `0.7`

Conflicting clades between `0.7` and `0.9` are preserved as
`moderate_support_disagreement` so the report does not flatten moderate support
into either strong conflict or weak noise. If the present-side tree did not
carry a parseable support label, the conflict row is marked
`support_unavailable`.

When the question is which shared taxa are driving disagreement, use
`compare influence`. That workflow performs a leave-one-taxon-out comparison
for every shared taxon, recomputes the topology and support conflict surfaces
after each exclusion, and ranks taxa by how much the disagreement surface
changes.

```bash
bijux-phylogenetics compare influence \
  artifacts/mammals.reference.nwk \
  artifacts/mammals.candidate.nwk \
  --out artifacts/mammals.taxon-influence.tsv \
  --json
```

The written ledger keeps one row per excluded taxon and preserves:
- whether the baseline rooted conflict disappeared or persisted after exclusion
- rooted and unrooted Robinson-Foulds deltas after the taxon was removed
- support-disagreement, conflicting-clade, and high-support-conflict deltas
- a transparent influence score and resulting rank

This review surface is useful when one taxon appears unstable, misplaced, or
poorly aligned and the practical question is whether disagreement is global or
concentrated around that single tip. The ranking is deliberately heuristic: it
adds the normalized topology shift to the absolute support-surface count shifts
so reviewers can see why one taxon outranks another instead of treating the
rank as an opaque diagnostic.

When the main question is which clades agree or conflict across several trees,
use `compare clades`. The command accepts two required tree paths plus
additional trees through repeated `--tree` flags, computes clade overlap on the
shared taxon set, preserves support values where the underlying tree labels
carry them, and can write one flat clade-overlap table with `--out`.

```bash
bijux-phylogenetics compare clades \
  artifacts/mammals.reference.nwk \
  artifacts/mammals.bootstrap.nwk \
  --tree artifacts/mammals.fasttree.nwk \
  --out artifacts/mammals.clade-overlap.tsv \
  --json
```

The JSON summary separates clades present in every tree from clades that are
conflicting or tree-specific. The written table then gives one row per
clade-per-tree observation, including whether that clade is present in the
given tree and which support value was observed when the tree format exposed
one.

When the goal is not overlap comparison but a direct ledger of every clade from
every sampled tree, use `tree-set clades`. That command writes one row per
clade observation per tree, preserving the tree index alongside the same clade
membership, support, branch-length, depth, and optional metadata summaries
used by the single-tree surface.

```bash
bijux-phylogenetics tree-set clades \
  artifacts/mammals.posterior.nwk \
  --metadata artifacts/mammals.metadata.tsv \
  --metadata-column species \
  --out artifacts/mammals.posterior-clades.tsv \
  --json
```

This surface is for reviewer-readable extraction, not just summary counts. It
is especially useful before consensus-building because it keeps minority and
tree-specific clades visible instead of collapsing them into one consensus or
frequency table immediately.

When the question is how shape varies across many sampled trees rather than how
individual clades vary, use `tree-set shape`. That command writes one row per
tree with the same balance and height metrics as `topology shape`, and its JSON
aggregate counts how many sampled trees are balanced, skewed, ladderized,
star-like, or comb-like.

```bash
bijux-phylogenetics tree-set shape \
  artifacts/mammals.posterior.nwk \
  --out artifacts/mammals.posterior-shape.tsv \
  --json
```

This surface is for tree-set review before or alongside consensus-building. It
keeps the distribution of imbalance and ladderization explicit instead of
reducing a posterior or bootstrap set to one representative tree immediately.

When the question is how branch-length distributions vary across many sampled
trees, use `tree-set branch-lengths`. That command writes one row per
non-root branch per tree and preserves the source tree index so long branches,
zero branches, and missing lengths can be traced back to the individual sample
that produced them.

```bash
bijux-phylogenetics tree-set branch-lengths \
  artifacts/mammals.posterior.nwk \
  --out artifacts/mammals.posterior-branch-lengths.tsv \
  --json
```

The JSON aggregate keeps the set-wide minimum, maximum, median, zero-length
count, negative-length count, and long-outlier count explicit. That makes it
possible to review whether one sampled tree is distorting the branch-length
distribution instead of assuming the full set is numerically homogeneous.

When the input is a bootstrap replicate tree file and the goal is one governed
review bundle rather than a chain of separate tree-set commands, use
`tree-set bootstrap-summary`. The command reads the bootstrap trees directly,
builds a consensus tree, computes clade frequencies and topology diversity, and
exports a dedicated unstable-branch ledger for consensus branches that are not
yet robust.

```bash
bijux-phylogenetics tree-set bootstrap-summary \
  artifacts/mammals.bootstrap.ufboot \
  --out-dir artifacts/mammals-bootstrap-review \
  --prefix mammals-bootstrap \
  --json
```

The written artifact set includes:
- `mammals-bootstrap.summary.tsv`
- `mammals-bootstrap.consensus.nwk`
- `mammals-bootstrap.clade-frequencies.tsv`
- `mammals-bootstrap.unstable-branches.tsv`
- `mammals-bootstrap.unstable-clades.tsv`
- `mammals-bootstrap.distance-matrix.tsv`
- `mammals-bootstrap.topology-clusters.tsv`

This bootstrap-specific bundle matters because majority-rule consensus alone can
hide whether its retained branches are only weakly recovered or compete with
clear alternative clades across the replicate set. The unstable-branch ledger
keeps that distinction explicit.

When your starting point is one aligned FASTA per locus, run
`alignment concatenate` first. That workflow writes the concatenated alignment,
the remapped partition file, and the taxon-by-locus occupancy matrix in one
step while preserving taxon identifiers and inserting `?` blocks for absent
taxa.

```bash
bijux-phylogenetics alignment concatenate loci/alpha-dna.fasta \
  loci/beta-protein.fasta \
  loci/gamma-dna.fasta \
  --data-type DNA \
  --data-type PROTEIN \
  --data-type DNA \
  --out artifacts/mixed-locus-supermatrix.aln.fasta \
  --partitions-out artifacts/mixed-locus-supermatrix.partitions.txt \
  --matrix-out artifacts/mixed-locus-supermatrix.matrix.tsv \
  --json
```

Use repeated `--data-type` flags when short or ambiguity-rich loci would make a
DNA locus and a protein locus look identical by characters alone. The
concatenation workflow records those explicit datatypes in the partition file
so downstream partitioned inference uses the honest locus contract.

For concatenated phylogenomics inputs, run `alignment occupancy` against an
aligned FASTA plus partition file before tree inference. The command can emit
per-taxon and per-locus TSV tables, a taxon-by-locus occupancy matrix, and a
retained alignment plus remapped partition file after applying explicit
coverage thresholds.

Use `--minimum-locus-occupancy` when partial locus fragments should count as
absent for taxon/locus thresholding. This keeps the retained matrix honest on
datasets where one or two recovered characters should not be treated as
meaningful locus presence. The TSV outputs include `site_coverage_fraction`
columns so the binary coverage calls can be reviewed against overall retained
signal.

```bash
bijux-phylogenetics alignment occupancy \
  artifacts/mixed-locus-supermatrix.aln.fasta \
  artifacts/mixed-locus-supermatrix.partitions.txt \
  --taxon-coverage-threshold 0.6 \
  --locus-coverage-threshold 0.6 \
  --minimum-locus-occupancy 0.75 \
  --taxa-out artifacts/occupancy/taxa.tsv \
  --loci-out artifacts/occupancy/loci.tsv \
  --matrix-out artifacts/occupancy/matrix.tsv \
  --filtered-alignment-out artifacts/occupancy/filtered.fasta \
  --filtered-partitions-out artifacts/occupancy/filtered-partitions.txt \
  --json
```

When the matrix is already aligned and you need to validate the partition file
itself, run `alignment partition-summary` first. That command reports assigned
and unassigned sites, mixed declared datatypes, and one row per locus, and it
can write the review table directly as TSV.

## Partitioned Multi-Locus Inference

Use the concatenation workflow first, then run the partition summary command
before sending the resulting matrix into the adapter inference surface. Pass
the same partition file into the adapter step that needs it.

```bash
bijux-phylogenetics alignment partition-summary \
  artifacts/mixed-locus-supermatrix.aln.fasta \
  artifacts/mixed-locus-supermatrix.partitions.txt \
  --out artifacts/multilocus.partition-summary.tsv \
  --json
bijux-phylogenetics adapter model-select \
  artifacts/mixed-locus-supermatrix.aln.fasta \
  --partitions artifacts/mixed-locus-supermatrix.partitions.txt \
  --out-dir artifacts/multilocus-model \
  --prefix multilocus \
  --json
bijux-phylogenetics adapter infer-ml \
  artifacts/mixed-locus-supermatrix.aln.fasta \
  --partitions artifacts/mixed-locus-supermatrix.partitions.txt \
  --out-dir artifacts/multilocus-ml \
  --model GTR+G \
  --prefix multilocus \
  --json
```

If every partition is DNA or every partition is protein, the adapter passes a
normalized partition scheme to IQ-TREE on the original aligned matrix. If the
partition file mixes DNA and protein loci, the adapter writes one extracted
alignment per partition and a generated NEXUS scheme before invoking IQ-TREE.
For that mixed-datatype path, do not force one fixed single model across every
partition. Use a model-selection keyword such as `MF`, `MFP`, `TEST`, or
`TESTMERGE` so the engine can choose partition-appropriate models honestly.

Use `adapter mrbayes-prepare` when you need a Bayesian NEXUS input file for
one aligned matrix and, optionally, one same-datatype partition file. The
command writes a ready-to-run MrBayes file with the data block, charset
definitions, partition declaration, model commands, and MCMC settings in one
step.

```bash
bijux-phylogenetics adapter mrbayes-prepare \
  artifacts/mixed-locus-supermatrix.aln.fasta \
  --partitions artifacts/mixed-locus-supermatrix.partitions.txt \
  --out artifacts/multilocus-bayesian.nex \
  --model gtr \
  --rates gamma \
  --ngen 200000 \
  --samplefreq 100 \
  --printfreq 100 \
  --json
```

When partitions are supplied, the preparation surface validates that the
coordinates fit the alignment and that any declared partition datatypes still
match the inferred alignment alphabet. The generated MrBayes block uses named
charsets plus one active partition declaration so the resulting file is
accepted by MrBayes as a partitioned Bayesian analysis input instead of as a
flat unpartitioned matrix.

Those direct adapter runs now preserve IQ-TREE's own `.iqtree` and `.log`
artifacts, plus the model-selection sidecar, a generated
`.model-candidates.tsv`, and `.treefile`, `.ufboot`, or `.contree` where the
invoked step produces them. The emitted JSON and manifests also include the
parsed `selected_model`, `selected_criterion`, `candidate_model_count`,
`best_model_aic`, `best_model_aicc`, `best_model_bic`, `log_likelihood`, and
support-value counts so reviewers can verify the ML result against structured
fields instead of manually scraping engine text files.

Use `adapter beast-prepare` when you need a real BEAST2 XML template from one
aligned matrix plus optional dating metadata. The command writes a BEAST2-style
XML file with an explicit alignment block, a starting tree, a strict or
uncorrelated lognormal clock, a Yule or birth-death tree prior, MCMC loggers,
and one MRCA prior per validated calibration target.

```bash
bijux-phylogenetics adapter beast-prepare \
  artifacts/mixed-locus-supermatrix.aln.fasta \
  --tree artifacts/mixed-locus-guide.nwk \
  --calibrations artifacts/fossil-calibrations.tsv \
  --tip-dates artifacts/tip-dates.tsv \
  --out artifacts/multilocus-beast.xml \
  --clock-model relaxed-lognormal \
  --tree-prior birth-death \
  --chain-length 2000000 \
  --log-every 1000 \
  --json
```

When `--calibrations` or `--tip-dates` are supplied, `--tree` is required so
the prepared XML can use the same taxa and named clades that were validated
during preparation. If the alignment is nucleotide or RNA, the template uses
an HKY site model; if the alignment is protein, it uses a JTT site model. The
written XML also names the default BEAST run logs as `STEM.$(seed).log` and
`STEM.$(seed).trees` so later execution produces a predictable artifact set
beside the XML.

Calibration handling is explicit rather than silent. If a calibration table
already provides both lower and upper bounds, the template preserves those hard
bounds as a BEAST uniform prior. If a calibration provides only a lower bound,
the template emits an explicit offset density above that minimum bound and
records a preparation warning in JSON so reviewers can see that the prior shape
was translated for template generation instead of copied as a literal hard
uniform interval.

Tip-dated workflows also surface one tree-prior boundary explicitly. If you
ask for the standard `birth-death` prior together with tip dates, the template
still writes valid XML, but the JSON warnings mark that combination as
exploratory because BEAST's own validator reports that the standard birth-death
prior is not serial-sampling aware.

Use `adapter beast-log` once a BEAST run has produced a posterior log and you
need a governed review surface instead of manual spreadsheet work. The command
parses native BEAST log files, accepts BEAST's own `Sample` state header, and
can discard an explicit burn-in fraction before computing summary statistics
and effective sample sizes.

```bash
bijux-phylogenetics adapter beast-log \
  artifacts/multilocus-beast.1.log \
  --burnin-fraction 0.1 \
  --summary-out artifacts/multilocus-beast.log-summary.tsv \
  --json
```

The JSON output keeps the raw parsed log plus a summary block that separates
posterior, likelihood, prior, clock, and tree-related parameters into explicit
lists. The optional `--summary-out` table writes one row per sampled parameter
with the post-burn-in mean, median, sample standard deviation, 95% HPD
interval, min/max range, first-half and second-half means, standardized drift,
and effective sample size so reviewers can audit BEAST log behaviour without
re-parsing the engine file themselves.

Use `adapter beast-parameters` when the main goal is posterior parameter
diagnostics rather than raw log parsing.

```bash
bijux-phylogenetics adapter beast-parameters \
  artifacts/multilocus-beast.1.log \
  --burnin-fraction 0.1 \
  --summary-out artifacts/multilocus-beast.parameter-summary.tsv \
  --json
```

That command emits one burn-in-aware JSON report and, when requested, one TSV
table covering posterior mean, median, sample standard deviation, 95% HPD
interval, effective sample size, and retained state window for every parameter
that survived the requested burn-in cut.

Use `adapter beast-convergence` when you want the same burn-in handling applied
to convergence warnings directly.

```bash
bijux-phylogenetics adapter beast-convergence \
  artifacts/multilocus-beast.1.log \
  --burnin-fraction 0.1 \
  --ess-threshold 200 \
  --mean-shift-threshold 0.5 \
  --json
```

That command emits structured warnings for low ESS and mean drift after the
requested burn-in has been discarded, so convergence review stays aligned with
the same retained posterior window described in the summary table.

Use `adapter beast-burnin-sensitivity` when you want to compare posterior
parameter estimates and clade probabilities across the governed burn-in review
fractions `5%`, `10%`, `25%`, and `50%`, or across an explicit custom set.

```bash
bijux-phylogenetics adapter beast-burnin-sensitivity \
  artifacts/multilocus-beast.1.trees \
  --log artifacts/multilocus-beast.1.log \
  --slice-out artifacts/multilocus-beast.burnin-slices.tsv \
  --parameter-out artifacts/multilocus-beast.burnin-parameters.tsv \
  --clade-out artifacts/multilocus-beast.burnin-clades.tsv \
  --json
```

That workflow writes one per-fraction slice ledger plus separate parameter and
clade comparison ledgers. Parameter instability is reported when the tested
95% HPD intervals fail to share a common overlap. Clade instability is
reported when a clade posterior probability crosses the majority-rule
threshold across the tested burn-in fractions.

Use `adapter beast-trees` once a BEAST run has produced a posterior tree file
and you need a governed tree-sample surface rather than ad hoc NEXUS handling.
The command parses native `.trees` files, keeps the sampled `STATE_*`
generations, applies an explicit burn-in fraction, extracts clade frequencies,
and can write the retained trees back out as normalized Newick for downstream
MCC and consensus workflows.

```bash
bijux-phylogenetics adapter beast-trees \
  artifacts/multilocus-beast.1.trees \
  --burnin-fraction 0.1 \
  --tree-set-out artifacts/multilocus-beast.postburnin.nwk \
  --json
```

The JSON summary reports the total sampled tree count, the discarded burn-in
count, the retained tree count, rooted-tree count, sampled state count, and
the extracted clade table. The optional `--tree-set-out` file is a normalized
Newick tree set that can be handed directly to the existing tree-set
consensus, topology-diversity, and MCC review surfaces without re-scraping the
original BEAST NEXUS container.

Use `adapter beast-subsample` when the posterior tree file is already too large
for the next review step, but you still want a governed retained subset instead
of ad hoc copying. The command supports either evenly spaced thinning or a
seeded random subset and preserves the native `STATE_*` labels in the retained
sample ledger.

```bash
bijux-phylogenetics adapter beast-subsample \
  artifacts/multilocus-beast.1.trees \
  --method evenly-spaced \
  --burnin-fraction 0.1 \
  --thinning-interval 5 \
  --tree-set-out artifacts/multilocus-beast.subsample.nwk \
  --sample-table-out artifacts/multilocus-beast.subsample.tsv \
  --json
```

For random retained subsets, switch to `--method random`, pass
`--sample-count`, and set `--seed` explicitly when you need the same retained
trees again later. The retained tree-set file stays normalized Newick, while
the TSV ledger records the retained source index, post-burn-in index, tree
name, sampled state, and rooted flag for each kept tree.

Use `adapter beast-consensus` when you want a governed majority-rule summary
tree from BEAST posterior samples instead of chaining generic tree-set commands
by hand.

```bash
bijux-phylogenetics adapter beast-consensus \
  artifacts/multilocus-beast.1.trees \
  --burnin-fraction 0.1 \
  --out artifacts/multilocus-beast.consensus.nwk \
  --tree-set-out artifacts/multilocus-beast.postburnin.nwk \
  --clade-table-out artifacts/multilocus-beast.clades.tsv \
  --json
```

That command applies burn-in filtering once, writes the retained posterior tree
set, computes informative clade frequencies, and writes a consensus tree whose
internal labels are posterior clade probabilities on the `0..1` scale. The
clade-frequency ledger preserves every informative retained clade, not only the
majority clades that appear in the consensus topology, so reviewers can see
which alternative groupings remained credible after burn-in removal.

Use `adapter beast-diversity` when you want a governed uncertainty view over
posterior topology dispersion rather than only a single consensus summary.

```bash
bijux-phylogenetics adapter beast-diversity \
  artifacts/multilocus-beast.1.trees \
  --burnin-fraction 0.1 \
  --tree-set-out artifacts/multilocus-beast.postburnin.nwk \
  --distance-out artifacts/multilocus-beast.distances.tsv \
  --topology-out artifacts/multilocus-beast.topologies.tsv \
  --unstable-clade-out artifacts/multilocus-beast.unstable-clades.tsv \
  --json
```

That workflow keeps the retained posterior tree set, writes the full pairwise
RF distance ledger, clusters retained trees by rooted topology, and exports an
unstable-clade ledger for non-unanimous groupings. The JSON summary then adds
the retained tree count, number of unique rooted topologies, dominant topology
frequency, effective topology count, RF pair count, and unstable-clade count
so reviewers can judge whether the posterior is tightly concentrated or broadly
dispersed before treating one consensus tree as sufficient.

Use `adapter mrbayes-run` after preparation when you want the governed runtime
to execute MrBayes and preserve the native posterior artifacts for later
inspection instead of leaving the engine outputs implicit.

```bash
bijux-phylogenetics adapter mrbayes-run \
  artifacts/multilocus-bayesian.nex \
  --json
```

The run keeps the sampled posterior trees (`.run1.t`), parameter traces
(`.run1.p`), MCMC diagnostics (`.mcmc`), and consensus tree (`.con.tre`) in
the same output directory. Those files are then consumable through the parser
surfaces instead of through manual text scraping:

- `adapter mrbayes-traces` for tabular parameter traces from `.run1.p`
- `adapter mrbayes-parameters` for burn-in-aware posterior mean, median, SD, 95% HPD, and ESS summaries from `.run1.p`
- `adapter mrbayes-trees` for sampled posterior trees and generation tags from `.run1.t`
- `adapter mrbayes-subsample` for governed retained posterior subsets with generation metadata from `.run1.t`
- `adapter mrbayes-mcmc` for acceptance-rate and split-frequency diagnostics from `.mcmc`
- `adapter mrbayes-consensus` for the annotated consensus topology and posterior-probability range from `.con.tre`

That separation matters because the consensus-tree file uses MrBayes-specific
annotation syntax that common generic NEXUS readers may not accept directly.
The governed parser strips the native inline annotations only after recording
their posterior-probability values, so reviewers keep both a clean tree object
and the support summary that MrBayes actually emitted.

When you want posterior parameter summaries directly, use the dedicated
diagnostics surface:

```bash
bijux-phylogenetics adapter mrbayes-parameters \
  artifacts/multilocus-bayesian.nex.run1.p \
  --burnin-fraction 0.25 \
  --summary-out artifacts/multilocus-bayesian.parameter-summary.tsv \
  --json
```

That workflow writes one row per retained parameter with posterior mean,
median, sample standard deviation, 95% HPD interval, effective sample size,
first-half versus second-half mean split, and the retained generation window.

Use `adapter mrbayes-burnin-sensitivity` when you want the same cross-fraction
review for MrBayes posterior trees and traces.

```bash
bijux-phylogenetics adapter mrbayes-burnin-sensitivity \
  artifacts/multilocus-bayesian.nex.run1.t \
  --traces artifacts/multilocus-bayesian.nex.run1.p \
  --slice-out artifacts/multilocus-bayesian.burnin-slices.tsv \
  --parameter-out artifacts/multilocus-bayesian.burnin-parameters.tsv \
  --clade-out artifacts/multilocus-bayesian.burnin-clades.tsv \
  --json
```

This workflow tests the same governed default fractions `5%`, `10%`, `25%`,
and `50%` unless you pass an explicit custom set. It reports parameter
instability from non-overlapping 95% HPD intervals and clade instability from
posterior probabilities that move across the majority-rule threshold.

Use `adapter mrbayes-subsample` when a downstream audit only needs a retained
posterior subset instead of the full `.run1.t` file.

```bash
bijux-phylogenetics adapter mrbayes-subsample \
  artifacts/multilocus-bayesian.nex.run1.t \
  --method random \
  --burnin-fraction 0.25 \
  --sample-count 200 \
  --seed 7 \
  --tree-set-out artifacts/multilocus-bayesian.subsample.nwk \
  --sample-table-out artifacts/multilocus-bayesian.subsample.tsv \
  --json
```

That workflow supports both evenly spaced thinning and seeded random
subsampling. The retained tree-set file is normalized Newick, and the sample
ledger keeps the retained source index, post-burn-in index, tree name,
sampled generation, and rooted flag so the subset remains traceable back to the
native MrBayes posterior output.

For ultrafast bootstrap review specifically, `adapter bootstrap` now writes
three reviewer-facing TSV artifacts alongside the native IQ-TREE files:

- `PREFIX.support.tsv` for every supported internal branch
- `PREFIX.low-support.tsv` for the subset of branches below the governed weak-support threshold
- `PREFIX.support-histogram.tsv` for the support distribution buckets used in reports

Those artifacts map support values back onto explicit descendant-taxon clades,
flag low-support branches directly, and preserve the same `lt50`, `50to69`,
`70to89`, and `ge90` buckets that appear in the manifest and HTML workflow
report.

Use `adapter infer-fast` when you need a rapid approximate tree for one aligned
DNA or protein matrix and you want reviewer-facing local-support evidence
instead of a bare Newick file.

```bash
bijux-phylogenetics adapter infer-fast \
  aligned-matrix.fasta \
  --out artifacts/mammals.fasttree.nwk \
  --sequence-type protein \
  --json
```

That workflow runs FastTree directly, writes the inferred tree to the requested
output path, and emits three sidecar TSV artifacts next to it:

- `mammals.fasttree.support.tsv` for every internal branch with a parsable local-support label
- `mammals.fasttree.low-support.tsv` for the subset of branches below the governed weak-support threshold
- `mammals.fasttree.support-histogram.tsv` for the local-support distribution buckets used in reports

The structured manifest and JSON also expose the FastTree approximation
contract explicitly: the method is approximately maximum-likelihood, the native
support labels are SH-like local-support proportions, and the governed support
scale is `0..1` rather than bootstrap percentages.

Use `adapter infer-large` when the matrix is already aligned and you need a
large-alignment inference path that avoids copying the matrix through multiple
Python-side structures before FastTree runs.

```bash
bijux-phylogenetics adapter infer-large \
  aligned-matrix.fasta \
  --out-dir artifacts/large-alignment \
  --prefix mammals \
  --sequence-type protein \
  --timeout-seconds 600 \
  --resume \
  --json
```

That workflow performs a streamed preflight scan of the aligned FASTA, runs
FastTree in place on the original matrix, and writes these review outputs:

- `mammals.tree`
- `mammals.support.tsv`
- `mammals.low-support.tsv`
- `mammals.support-histogram.tsv`
- `mammals.resources.tsv`
- `mammals.log`
- `mammals.manifest.json`

The streamed preflight records sequence count, alignment width, total site
cells, and inferred sequence type without materializing the full matrix as an
in-memory Python alignment object. The resource ledger separates preflight
allocation observations from sampled FastTree process RSS so large runs expose
both wall time and memory pressure directly. When `--resume` is used, the
workflow reuses a completed manifest only if the input checksum, command, and
recorded outputs still agree; otherwise it reruns the inference step.

Use `adapter compare-engines` when you need a governed side-by-side comparison
between IQ-TREE and FastTree on the same aligned matrix.

```bash
bijux-phylogenetics adapter compare-engines \
  aligned-matrix.fasta \
  --out-dir artifacts/engine-comparison \
  --prefix mammals \
  --sequence-type dna \
  --bootstrap-replicates 1000 \
  --json
```

That workflow runs IQ-TREE model selection, IQ-TREE ultrafast bootstrap
support inference, and FastTree approximate inference on the same alignment.
It then writes these user-facing outputs:

- `mammals.fasttree.nwk`
- `mammals.iqtree-support.nwk`
- `mammals.comparison.html`
- `mammals.comparison.tsv`
- `mammals.shared-clades.tsv`
- `mammals.conflicting-clades.tsv`
- `mammals.manifest.json`

The shared-clade ledger preserves both engines' support values for clades that
appear in both trees. The conflicting-clade ledger separates clades that appear
in only one tree from shared clades whose normalized support fractions diverge
enough to merit review. The normalization rule is explicit and limited:
FastTree SH-like local support and IQ-TREE UFBoot support are shown together as
fractions only for side-by-side review, not as proof that the two support
methods are interchangeable.

Use `adapter reproducibility` when you need to test whether repeated
bootstrap-supported IQ-TREE inference stays deterministic under fixed settings.

```bash
bijux-phylogenetics adapter reproducibility \
  aligned-matrix.fasta \
  --out-dir artifacts/inference-reproducibility \
  --prefix mammals \
  --sequence-type dna \
  --bootstrap-replicates 1000 \
  --repeats 3 \
  --json
```

That workflow runs model selection once to choose a fixed model, then reruns
the same supported IQ-TREE inference settings multiple times on the same
alignment. It writes these review outputs:

- `mammals.runs.tsv`
- `mammals.comparisons.tsv`
- `mammals.support-deltas.tsv`
- `mammals.manifest.json`

The run ledger records one manifest and supported tree per rerun. The
comparison ledger classifies each rerun relative to the baseline as
`deterministic`, `equivalent`, or `unstable` after checking topology,
log-likelihood, and clade support values. The support-delta ledger preserves
the per-clade support shifts as normalized fractions so support drift can be
reviewed directly instead of being inferred from summary text alone.

Use `adapter sh-alrt` when you need SH-aLRT support alongside ultrafast
bootstrap support on the same supported tree.

```bash
bijux-phylogenetics adapter sh-alrt \
  aligned-matrix.fasta \
  --out-dir artifacts/sh-alrt-support \
  --model GTR+G \
  --alrt-replicates 1000 \
  --bootstrap-replicates 1000 \
  --prefix mammals \
  --json
```

That workflow runs IQ-TREE with both `-alrt` and `-bb`, retains the native
`.treefile`, `.ufboot`, `.iqtree`, and `.log` outputs, writes
`mammals.support.tsv` with both support measures on each branch, and writes
`mammals.conflicting-support.tsv` for branches where SH-aLRT and UFBoot imply
different confidence postures under the governed thresholds. The structured
manifest and JSON also expose parsed SH-aLRT minima, maxima, annotated-branch
counts, and conflicting-signal counts directly.

## Coding DNA Alignment

Use `adapter align --codon-aware` when you need a nucleotide alignment that
preserves codon triplets for downstream phylogenetic inference.

```bash
bijux-phylogenetics adapter align coding-cds.fasta \
  --out artifacts/coding.aligned.fasta \
  --mode linsi \
  --codon-aware \
  --json
bijux-phylogenetics adapter model-select artifacts/coding.aligned.fasta \
  --out-dir artifacts/coding-model \
  --prefix coding \
  --sequence-type dna \
  --json
bijux-phylogenetics adapter infer-ml artifacts/coding.aligned.fasta \
  --out-dir artifacts/coding-ml \
  --model GTR+G \
  --prefix coding \
  --sequence-type dna \
  --json
```

The codon-aware alignment workflow excludes frame-broken sequences and
sequences with premature stop codons before MAFFT runs. It then aligns a
translated amino-acid guide and projects guide gaps back as nucleotide triplets,
so the final alignment length stays divisible by three and codon boundaries are
retained. The workflow also writes the translated guide input, the aligned
guide, and a TSV ledger of excluded sequences for review.

## Raw FASTA To Tree

Use `adapter fasta-to-tree` when you need one governed command from unaligned
FASTA to a reviewable inference bundle.

```bash
bijux-phylogenetics alignment sequence-type raw-sequences.fasta --json
bijux-phylogenetics alignment validate-input raw-sequences.fasta --json
bijux-phylogenetics alignment repair-input raw-sequences.fasta \
  --out artifacts/raw-sequences.repaired.fasta \
  --normalize-identifiers \
  --remove-invalid-records \
  --json
bijux-phylogenetics adapter fasta-to-tree raw-sequences.fasta \
  --out-dir artifacts/fasta-to-tree \
  --prefix mammals \
  --iqtree-seed 1 \
  --iqtree-threads 1 \
  --normalize-identifiers \
  --remove-invalid-records \
  --bootstrap-replicates 1000 \
  --json
```

The workflow accepts DNA and protein FASTA inputs. It writes these durable
user-facing outputs:

- `mammals.aln`
- `mammals.trimmed.aln`
- `mammals.tree`
- `mammals.log`
- `mammals.model.tsv`
- `mammals.support.tsv`
- `mammals.manifest.json`

It also retains step-specific engine artifacts under
`artifacts/fasta-to-tree/engine-artifacts/mammals/` for auditability.

Within that engine-artifact directory, the IQ-TREE stages preserve the native
`.iqtree`, `.log`, model-selection sidecar, `.model-candidates.tsv`,
`.treefile`, `.ufboot`, and `.contree` files where each stage produces them.
The bootstrap-support stage also exports `.support.tsv`, `.low-support.tsv`,
and `.support-histogram.tsv` so the supported tree can be reviewed without
re-parsing Newick labels manually. The corresponding manifest entries expose
the parsed selected model, selected criterion, AIC/AICc/BIC winners,
candidate-model counts, log-likelihood, support summary, and weak-backbone
summary directly.

The dedicated SH-aLRT support workflow follows the same pattern: it retains the
native IQ-TREE files, exports `.support.tsv` and
`.conflicting-support.tsv`, and records the parsed combined SH-aLRT/UFBoot
branch summary in the manifest.

The IQ-TREE part of the workflow now defaults to deterministic execution with
`--iqtree-seed 1` and `--iqtree-threads 1`. Ultrafast bootstrap support is the
governed support workflow here, so `--bootstrap-replicates` must be at least
`1000`.

Use `alignment sequence-type` when you need the raw FASTA type decision before
any engine is invoked. It reports compatible types, the selected default, the
confidence level, and mixed or invalid blocking signals.

Use `alignment validate-input` when you need one broader report of duplicate
identifiers, illegal sequence characters, empty records, raw-sequence length
outliers, and sequence-type detection before any engine is invoked.

Use `alignment repair-input` or the matching `adapter fasta-to-tree`
`--normalize-identifiers` and `--remove-invalid-records` controls when you want
the runtime to prepare a repaired FASTA explicitly rather than proceeding on
silent assumptions. Without those repair flags, `adapter fasta-to-tree` now
fails fast on duplicate identifiers, empty sequences, or illegal characters.
Mixed raw inputs also fail fast unless you declare a compatible
`--sequence-type` explicitly and remove incompatible records.

The checked real-dataset workflow corpus now lives under:

- `packages/bijux-phylogenetics/tests/fixtures/fasta_to_tree/real/`
- `packages/bijux-phylogenetics/tests/fixtures/expected/fasta_to_tree/`

Those checks pin reviewer-facing output bundles for:

- `gnathostome-ortholog-proteins`
- `gnathostome-ortholog-coding-sequences`
- `strnog-enog411bqtj-proteins`

They verify that the workflow emits the aligned matrix, trimmed matrix, tree,
log, model table, and support table with stable names on real DNA and protein
inputs.

For coding DNA, prefer the codon-aware alignment workflow above and then run
the downstream adapter steps explicitly. The current `adapter fasta-to-tree`
path is still a generic trimAl-based pipeline rather than the codon-preserving
entrypoint.
